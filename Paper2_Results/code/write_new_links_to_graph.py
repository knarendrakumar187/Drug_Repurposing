"""Load Paper 2 drug evidence as Drug -ACTS_ON-> Gene, without merging evidence.

Dry run by default. --apply and --undo require --backup to a NEW JSON file.
All writes and preservation checks use one transaction. Stop other graph
writers during this operation. Legacy graphs must first use
scripts/unify_drug_schema.py. Credentials are never saved in the backup.
"""
import argparse
from collections import Counter
from datetime import date
import json
import os
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.unify_drug_schema import snapshot, validate, needs_migration

LAYER = 'external (Idea 4)'  # Retained provenance value; not a label or relationship type.


def signature(rows):
    return Counter(json.dumps(r, sort_keys=True) for r in rows)


def load_plan(tx, inputs):
    m = pd.read_csv(inputs / 'new_links_medicines_curated.csv')
    if m.empty or not m.evidence.eq('curated drug-target').all() or m.duplicated(['Drug_ID', 'Gene_ID']).any():
        raise ValueError('Expected unique curated drug-target pairs.')
    raw = pd.read_csv(inputs / 'raw_rows_all_sources.csv.gz', low_memory=False, usecols=['Drug_key', 'chembl'])
    aliases = raw.dropna().groupby('Drug_key').chembl.agg(lambda s: sorted(set(s))).to_dict()
    releases = json.loads((inputs / 'numbers.json').read_text())['release']
    release = '; '.join(f'{k} {v}' for k, v in releases.items())
    candidates = [{'key': d, 'ids': [d] + [a for a in aliases.get(d, []) if a != d]} for d in sorted(m.Drug_ID.unique())]
    matches = tx.run('''UNWIND $rows AS row
        UNWIND range(0, size(row.ids)-1) AS i
        MATCH (d:Drug {node_id: row.ids[i]})
        WITH row.key AS key, i, d.node_id AS id ORDER BY i
        RETURN key, collect(id)[0] AS id''', rows=candidates).data()
    node_of = {r['key']: r['id'] for r in matches}
    genes = {r['id'] for r in tx.run('UNWIND $ids AS id MATCH(g:Gene {ensembl_id:id}) RETURN g.ensembl_id AS id', ids=sorted(m.Gene_ID.unique()))}
    if genes != set(m.Gene_ID):
        raise ValueError('Some input genes are absent from the graph.')
    def text(v):
        return '' if pd.isna(v) else str(v)
    new_nodes, rows = {}, []
    for r in m.itertuples():
        drug = node_of.get(r.Drug_ID, r.Drug_ID)
        if r.Drug_ID not in node_of:
            props = {'node_id': drug, 'name': text(r.Drug), 'drug_stage': text(r.drug_stage), 'source': LAYER, 'release': release}
            if pd.notna(r.max_phase): props['max_phase'] = float(r.max_phase)
            if pd.notna(r.molecule_type): props['molecule_type'] = str(r.molecule_type)
            # Match the historical writer's drop_duplicates('Drug_ID') rule:
            # first input row supplies node metadata; every edge keeps its own evidence.
            new_nodes.setdefault(drug, props)
        rows.append({'drug': drug, 'gene': r.Gene_ID, 'props': {
            'action_labels': text(r.action_labels), 'push': text(r.push),
            'sources': text(r.sources), 'source_detail': text(r.source_detail),
            'evidence': text(r.evidence), 'n_sources': int(r.n_sources),
            'drug_stage': text(r.drug_stage), 'mechanism': text(r.mechanism_text),
            'layer': LAYER, 'release': release}})
    if len({(r['drug'], r['gene']) for r in rows}) != len(rows):
        raise ValueError('Drug alias mapping collapses distinct input pairs; resolve before loading.')
    return {'nodes': list(new_nodes.values()), 'links': rows, 'release': release}


def existing_layer(tx):
    return tx.run('''MATCH(d:Drug)-[r:ACTS_ON {layer:$layer}]->(g:Gene)
        RETURN d.node_id AS drug, g.ensembl_id AS gene, properties(r) AS props''', layer=LAYER).data()


def layer_matches(actual, wanted):
    def norm(row):
        props = {k:v for k,v in row['props'].items() if k not in ('added', 'original_relationship_type', 'migration_record_id')}
        return {'drug': row['drug'], 'gene': row['gene'], 'props': props}
    return signature(map(norm, actual)) == signature(map(norm, wanted))


def add_links(tx, before, plan, added):
    if snapshot(tx) != before:
        raise ValueError('Graph changed since backup; retry with a fresh backup.')
    if needs_migration(before):
        raise ValueError('Convert legacy drug labels and relationships first.')
    if existing_layer(tx):
        raise ValueError('The Paper 2 layer already exists; verify it before loading again.')
    # Evidence records use CREATE, never MERGE by the drug/gene endpoints.
    for offset in range(0, len(plan['nodes']), 1000):
        nodes = [dict(p, added=added) for p in plan['nodes'][offset:offset+1000]]
        tx.run('UNWIND $rows AS row CREATE(d:Drug) SET d = row', rows=nodes).consume()
    for offset in range(0, len(plan['links']), 1000):
        rows = [dict(r, props=dict(r['props'], added=added)) for r in plan['links'][offset:offset+1000]]
        tx.run('''UNWIND $rows AS row MATCH(d:Drug {node_id:row.drug})
            MATCH(g:Gene {ensembl_id:row.gene}) CREATE(d)-[r:ACTS_ON]->(g)
            SET r = row.props''', rows=rows).consume()
    after = snapshot(tx)
    validate(after)
    if needs_migration(after) or not layer_matches(existing_layer(tx), plan['links']):
        raise ValueError('Added endpoint/property verification failed; rolling back.')
    nodes = {r['id']:r for r in after['nodes']}; links = {r['id']:r for r in after['links']}
    if any(nodes.get(r['id']) != r for r in before['nodes']) or any(links.get(r['id']) != r for r in before['links']):
        raise ValueError('Pre-existing drug evidence changed; rolling back.')
    if len(nodes) != len(before['nodes']) + len(plan['nodes']) or len(links) != len(before['links']) + len(plan['links']):
        raise ValueError('Added counts do not match the plan; rolling back.')
    return {'drug_nodes_created':len(plan['nodes']), 'links_written':len(plan['links']),
            'relationship_type':'ACTS_ON', 'new_node_label':'Drug', 'existing_evidence_preserved':True,
            'links_match_file':True, 'verified':True}


def undo(tx, before):
    if snapshot(tx) != before:
        raise ValueError('Graph changed since backup.')
    doomed = tx.run('''MATCH()-[r:ACTS_ON|EXTERNAL_ACTS_ON {layer:$layer}]->()
        RETURN elementId(r) AS id''', layer=LAYER).value()
    tx.run('''MATCH()-[r:ACTS_ON|EXTERNAL_ACTS_ON {layer:$layer}]->() DELETE r''', layer=LAYER).consume()
    removed_nodes = tx.run('''MATCH(d) WHERE (d:Drug OR d:ExternalDrug) AND d.source=$layer AND NOT (d)--()
        WITH d, elementId(d) AS id DELETE d RETURN id''', layer=LAYER).value()
    after = snapshot(tx)
    if after['links'] != [r for r in before['links'] if r['id'] not in set(doomed)] or after['nodes'] != [r for r in before['nodes'] if r['id'] not in set(removed_nodes)]:
        raise ValueError('Undo changed unrelated evidence; rolling back.')
    return {'links_removed':len(doomed), 'orphan_drugs_removed':len(removed_nodes), 'verified':True}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    modes = ap.add_mutually_exclusive_group()
    modes.add_argument('--apply', action='store_true'); modes.add_argument('--undo', action='store_true')
    ap.add_argument('--backup', type=Path)
    ap.add_argument('--inputs', type=Path, default=ROOT/'Paper2_Results/New_drug_gene_links')
    ap.add_argument('--out', type=Path, default=ROOT/'Paper2_Results/New_drug_gene_links/graph_write')
    ap.add_argument('--uri', default=os.environ.get('NEO4J_URI', 'bolt://localhost:7687'))
    ap.add_argument('--database', default=os.environ.get('NEO4J_DATABASE', 'neo4j'))
    ap.add_argument('--user', default=os.environ.get('NEO4J_USERNAME', os.environ.get('NEO4J_USER', 'neo4j')))
    args = ap.parse_args()
    if (args.apply or args.undo) and not args.backup: ap.error('Graph writes require --backup to a new file.')
    if not os.environ.get('NEO4J_PASSWORD'): ap.error('Set NEO4J_PASSWORD in the environment.')
    from neo4j import GraphDatabase
    with GraphDatabase.driver(args.uri, auth=(args.user, os.environ['NEO4J_PASSWORD'])) as drv, drv.session(database=args.database) as s:
        before = s.execute_read(snapshot); validate(before)
        if not args.undo and needs_migration(before): ap.error('Run scripts/unify_drug_schema.py first; do not rerun the loader on legacy data.')
        plan = None if args.undo else s.execute_read(load_plan, args.inputs)
        if plan:
            actual = s.execute_read(existing_layer)
            if actual:
                if not layer_matches(actual, plan['links']): raise ValueError('Existing layer differs from the input; no changes made.')
                print('Already loaded as Drug / ACTS_ON; input evidence matches. No changes made.'); return
            print(json.dumps({'mode':'apply' if args.apply else 'dry run', 'new_drugs':len(plan['nodes']), 'links':len(plan['links'])}))
        if not (args.apply or args.undo): return
        args.backup.parent.mkdir(parents=True, exist_ok=True)
        with args.backup.open('x', encoding='utf-8') as f: json.dump(before, f, indent=2)
        result = s.execute_write(undo, before) if args.undo else s.execute_write(add_links, before, plan, date.today().isoformat())
        args.out.mkdir(parents=True, exist_ok=True)
        filename = 'undo_log.json' if args.undo else 'write_log.json'
        (args.out/filename).write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
