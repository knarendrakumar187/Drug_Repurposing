"""Reconstruct the complete Paper 2 graph from the public Paper 1 graph.

Read-only plan by default. --apply requires an explicit URI, database and NEW
backup file. Run against a separate Paper 2 database; never the Paper 1 graph.
KEGG inputs remain local. Saved drug evidence is versioned in this repository.
Every graph change and its verification occur in one transaction.
"""
import argparse
from collections import Counter
from datetime import date
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys

import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'scripts'))
import restore_kegg
from scripts.unify_drug_schema import snapshot,validate,needs_migration
spec=importlib.util.spec_from_file_location('drug_loader',ROOT/'Paper2_Results/code/write_new_links_to_graph.py')
drug_loader=importlib.util.module_from_spec(spec);spec.loader.exec_module(drug_loader)
BASE=(126294,35140166);FINAL=(126999,35160254)
DIRECTION_FIELDS=['external_action_direction','external_action_label','external_action_sources','external_action_release','external_action_added']
INDICATION_LAYER='ChEMBL indication (added drugs)'
DIRECTION_RELEASE='Guide to Pharmacology 2026.3.1; DGIdb 5 (API, Oct 2026); ChEMBL REST API (Oct 2026)'


def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()


def counts(tx):
    return (tx.run('MATCH(n) RETURN count(n)').single()[0],tx.run('MATCH()-[r]->() RETURN count(r)').single()[0])


def directions(path):
    d=pd.read_csv(path)
    if len(d)!=297 or d.duplicated(['Drug_ID','Gene_ID']).any():raise ValueError('Expected 297 unique direction records.')
    return [{'drug':r.Drug_ID,'gene':r.Gene_ID,'props':{
        'external_action_direction':r.push,'external_action_label':str(r.action_label)[:300],
        'external_action_sources':r.direction_source.replace('external: ',''),'external_action_release':DIRECTION_RELEASE}} for r in d.itertuples()]


def indications(path):
    p=pd.read_csv(path);p=p[p.status.eq('write')]
    if len(p)!=2087 or p.duplicated(['Drug_ID','efo_id_graph']).any():raise ValueError('Expected 2,087 unique saved indication pairs.')
    def clean(v):return '' if pd.isna(v) else str(v)
    return [{'drug':r.Drug_ID,'disease':r.efo_id_graph,'props':{
        'source':'ChEMBL indication','layer':INDICATION_LAYER,'max_phase_for_ind':float(r.max_phase_for_ind),
        'efo_terms':clean(r.efo_terms),'mesh_headings':clean(r.mesh_headings),'drugind_ids':clean(r.drugind_ids),
        'chembl_release':str(r.chembl_release),'chembl_release_date':str(r.chembl_release_date),'retrieved':str(r.retrieved)}} for r in p.itertuples()]


def gap_properties(path):
    d=pd.read_csv(path,keep_default_na=False)
    if len(d)!=28 or d.Drug_ID.duplicated().any():raise ValueError('Expected 28 saved gap-fill drug records.')
    rows=[]
    for r in d.to_dict('records'):
        props={k:v for k,v in r.items() if k!='Drug_ID' and v!=''}
        if 'max_phase' in props:props['max_phase']=float(props['max_phase'])
        rows.append({'drug':r['Drug_ID'],'props':props})
    return rows


def make_plan(tx,gmt,ctd):
    n=counts(tx)
    if n not in (BASE,FINAL):raise ValueError(f'Expected a clean public Paper 1 import {BASE} or complete Paper 2 graph {FINAL}; found {n}.')
    before=snapshot(tx);validate(before)
    if needs_migration(before):raise ValueError('Legacy schema found. Use the separate migration tool first.')
    res=ROOT/'Paper2_Results'
    files={'kegg_gmt':gmt,'ctd':ctd,'enrichment':ROOT/'Paper1_Results/Task4/task4_all_significant.csv',
        'drug_links':res/'New_drug_gene_links/new_links_medicines_curated.csv',
        'drug_aliases':res/'New_drug_gene_links/raw_rows_all_sources.csv.gz',
        'drug_releases':res/'New_drug_gene_links/numbers.json',
        'directions':res/'Task7_external_layer/external_layer_links.csv',
        'indications':res/'Task4/graph_write/indication_pairs_plan.csv',
        'metadata':res/'Recovery2_ChEMBL/graph_gap_properties.csv'}
    for p in files.values():
        if not p.is_file():raise ValueError(f'Missing required input: {p}')
    kegg=restore_kegg.plan(tx,gmt,ctd,files['enrichment'])
    drugs=drug_loader.load_plan(tx,res/'New_drug_gene_links')
    if len(drugs['links'])!=5216 or (n==BASE and len(drugs['nodes'])!=519):raise ValueError('Saved drug records do not reproduce expected additions.')
    data={'before_counts':n,'before':before,'kegg':kegg,'drugs':drugs,'directions':directions(files['directions']),
          'indications':indications(files['indications']),'metadata':gap_properties(files['metadata']),
          'inputs':{k:{'file':p.name,'sha256':digest(p)} for k,p in files.items()}}
    # All indication diseases must already be in the imported baseline.
    ids={r['disease'] for r in data['indications']}
    found={r['id'] for r in tx.run('UNWIND $ids AS id MATCH(d:Disease {node_id:id}) RETURN d.node_id AS id',ids=sorted(ids))}
    if ids!=found:raise ValueError('Some indication diseases are absent; no changes made.')
    if n==BASE and drug_loader.existing_layer(tx):raise ValueError('Unexpected existing drug layer.')
    return data


def apply(tx,data):
    if counts(tx)!=tuple(data['before_counts']) or snapshot(tx)!=data['before']:raise ValueError('Graph changed since backup.')
    today=date.today().isoformat()
    existing=drug_loader.existing_layer(tx)
    if existing:
        if not drug_loader.layer_matches(existing,data['drugs']['links']):raise ValueError('Existing drug layer differs from saved evidence.')
    else:drug_loader.add_links(tx,data['before'],data['drugs'],today)
    # Fill only the historical ChEMBL gaps; retain all OptimusKG values.
    for row in data['metadata']:
        found=tx.run('MATCH(d:Drug {node_id:$id}) RETURN properties(d) AS p',id=row['drug']).data()
        if len(found)!=1:raise ValueError('Gap-fill drug is missing or ambiguous.')
        props=found[0]['p']
        if any(props.get(k) not in (None,'',v) for k,v in row['props'].items()):raise ValueError('Metadata gap fill would overwrite existing evidence.')
        tx.run('MATCH(d:Drug {node_id:$id}) SET d += $props',id=row['drug'],props=row['props']).consume()
    for row in data['directions']:
        found=tx.run('MATCH(d:Drug {node_id:$drug})-[r:ACTS_ON]->(g:Gene {ensembl_id:$gene}) RETURN properties(r) AS p',drug=row['drug'],gene=row['gene']).data()
        if len(found)!=1:raise ValueError('Direction record must match exactly one original drug link.')
        if any(found[0]['p'].get(k) not in (None,v) for k,v in row['props'].items()):raise ValueError('Conflicting direction evidence.')
        tx.run('MATCH(d:Drug {node_id:$drug})-[r:ACTS_ON]->(g:Gene {ensembl_id:$gene}) SET r += $props,r.external_action_added=coalesce(r.external_action_added,$today)',drug=row['drug'],gene=row['gene'],props=row['props'],today=today).consume()
    for row in data['indications']:
        found=tx.run('MATCH(d:Drug {node_id:$drug})-[r:TREATS]->(x:Disease {node_id:$disease}) RETURN properties(r) AS p',drug=row['drug'],disease=row['disease']).data()
        if found:
            if len(found)!=1 or any(found[0]['p'].get(k)!=v for k,v in row['props'].items()):raise ValueError('Existing indication evidence differs from the saved record.')
        else:
            n=tx.run('MATCH(d:Drug {node_id:$drug}) MATCH(x:Disease {node_id:$disease}) CREATE(d)-[r:TREATS]->(x) SET r=$props,r.added=$today RETURN count(r)',drug=row['drug'],disease=row['disease'],props=row['props'],today=today).single()[0]
            if n!=1:raise ValueError('Indication endpoint missing or ambiguous.')
    restore_kegg.apply(tx,data['kegg'])
    return verify(tx,data)


def verify(tx,data):
    n=counts(tx)
    if n!=FINAL:raise ValueError(f'Final counts mismatch: {n}, expected {FINAL}.')
    restore_kegg.verify(tx,data['kegg'])
    if not drug_loader.layer_matches(drug_loader.existing_layer(tx),data['drugs']['links']):raise ValueError('Added drug evidence mismatch.')
    after=snapshot(tx);validate(after)
    if needs_migration(after):raise ValueError('Legacy schema remains.')
    actual_nodes={n['props']['node_id']:n['props'] for n in after['nodes']}
    gaps={r['drug']:r['props'] for r in data['metadata']}
    for old in data['before']['nodes']:
        p=old['props'];wanted=dict(p);wanted.update(gaps.get(p['node_id'],{}))
        if actual_nodes.get(p['node_id'])!=wanted:raise ValueError('Original drug metadata was unexpectedly changed.')
    before_links={r['id']:r for r in data['before']['links']};after_links={r['id']:r for r in after['links']}
    action_expected={(r['drug'],r['gene']):r['props'] for r in data['directions']}
    action_actual=tx.run('MATCH(d:Drug)-[r:ACTS_ON]->(g:Gene) WHERE r.external_action_direction IS NOT NULL RETURN d.node_id AS drug,g.ensembl_id AS gene,properties(r) AS p,elementId(r) AS id').data()
    if len(action_actual)!=297:raise ValueError('Direction count mismatch.')
    allowed={}
    for r in action_actual:
        expected=action_expected.get((r['drug'],r['gene']))
        if expected is None or any(r['p'].get(k)!=v for k,v in expected.items()):raise ValueError('Direction property mismatch.')
        allowed[r['id']]=dict(expected,external_action_added=r['p']['external_action_added'])
    for id,old in before_links.items():
        expected=dict(old['props']);expected.update(allowed.get(id,{}))
        new=after_links.get(id)
        if not new or new['props']!=expected or new['start']!=old['start'] or new['end']!=old['end']:raise ValueError('Original drug-link evidence changed unexpectedly.')
    actual=tx.run('MATCH(d:Drug)-[r:TREATS {layer:$layer}]->(x:Disease) RETURN d.node_id AS drug,x.node_id AS disease,properties(r) AS props',layer=INDICATION_LAYER).data()
    def sig(rows):return Counter(json.dumps({'drug':r['drug'],'disease':r['disease'],'props':{k:v for k,v in r['props'].items() if k!='added'}},sort_keys=True) for r in rows)
    if sig(actual)!=sig(data['indications']):raise ValueError('Indication evidence mismatch.')
    return {'verified':True,'nodes':n[0],'relationships':n[1],'kegg_nodes':186,'kegg_memberships':12724,'kegg_disease_links':61,
            'drug_nodes':12544,'acts_on':25890,'added_drug_links':5216,'action_annotations':297,'added_indications':2087,'metadata_gap_drugs':28}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--kegg-gmt',type=Path,required=True);ap.add_argument('--ctd',type=Path,required=True)
    ap.add_argument('--uri',required=True);ap.add_argument('--database',required=True)
    ap.add_argument('--user',default=os.environ.get('NEO4J_USERNAME',os.environ.get('NEO4J_USER','neo4j')))
    ap.add_argument('--apply',action='store_true');ap.add_argument('--backup',type=Path)
    args=ap.parse_args()
    if args.apply and not args.backup:ap.error('--apply requires a NEW --backup file.')
    if not os.environ.get('NEO4J_PASSWORD'):ap.error('Set NEO4J_PASSWORD for the Paper 2 database.')
    from neo4j import GraphDatabase
    with GraphDatabase.driver(args.uri,auth=(args.user,os.environ['NEO4J_PASSWORD'])) as d,d.session(database=args.database) as s:
        data=s.execute_read(make_plan,args.kegg_gmt,args.ctd)
        if tuple(data['before_counts'])==FINAL:
            print(json.dumps(s.execute_read(verify,data),indent=2));print('Already complete; no changes made.');return
        print('Plan: restore KEGG; add 519 drugs, 5,216 drug links, 297 directions, 2,087 indications and the saved ChEMBL metadata gaps.')
        if not args.apply:print('Dry run: no graph changes.');return
        args.backup.parent.mkdir(parents=True,exist_ok=True)
        backup={'format':'paper2-additions-v1','before_counts':data['before_counts'],'drug_layer_before':data['before'],'inputs':data['inputs'],
                'scope':'Additive rebuild checkpoint, not a full database backup. Keep the imported dump for full restoration.'}
        with args.backup.open('x',encoding='utf-8') as f:json.dump(backup,f,indent=2)
        result=s.execute_write(apply,data)
        print(json.dumps(result,indent=2))


if __name__=='__main__':main()
