"""Restore only the KEGG data omitted from the public Paper 1 release.

Source: authorised MSigDB 2024.1 KEGG legacy GMT and CTD disease-pathway CSV.
No KEGG records are bundled with this repository. Functions are called by
rebuild_graph.py; matching/enrichment/CTD rules follow the original Paper 1 code.
"""
import re
import pandas as pd

ACRONYMS = {"ABC", "TGF", "ECM", "DNA", "RNA", "MAPK", "JAK", "STAT", "NOD", "RIG",
            "PPAR", "MTOR", "ERBB", "VEGF", "GNRH", "WNT", "NK", "T", "B", "II",
            "UV", "KRAS", "TNFA", "NFKB", "MYC", "E2F", "G2M", "UPR", "EMT",
            "SNARE", "GPI", "TCA", "FC", "HIF", "HCM", "ARVC", "SHH"}


def readable(term, prefix):
    """KEGG_TGF_BETA_SIGNALING_PATHWAY -> 'TGF beta signaling pathway'."""
    words = term[len(prefix):].split("_")
    # keep acronyms and gene-like words (with a digit) upper case
    out = [w if w in ACRONYMS or any(c.isdigit() for c in w) else w.lower()
           for w in words]
    s = " ".join(out)
    return s[0].upper() + s[1:]


def norm(value):
    return re.sub(r"[^A-Z0-9]+", "_", str(value).upper()).strip("_")


def plan(tx, gmt, ctd, enrichment):
    symbols = {r['symbol'] for r in tx.run('MATCH(g:Gene) WHERE g.symbol IS NOT NULL RETURN g.symbol AS symbol')}
    nodes, members = [], []
    with open(gmt, encoding='utf-8') as f:
        for line in f:
            parts=line.rstrip('\n').split('\t');term=parts[0];genes=sorted(set(parts[2:]))
            if len(parts)<3 or not term.startswith('KEGG_'):raise ValueError('Expected human KEGG legacy GMT.')
            node='KEGG:'+term
            nodes.append({'node_id':node,'name':readable(term,'KEGG_'),'msigdb_id':term,'source':'KEGG_legacy','msigdb_version':'2024.1','n_genes':len(genes)})
            members.extend({'symbol':g,'pathway':node} for g in genes if g in symbols)
    if len(nodes)!=186 or len(members)!=12724:raise ValueError('KEGG source does not reproduce the expected 186 pathways / 12,724 memberships.')
    e=pd.read_csv(enrichment)
    e=e[(~e.group.str.startswith('inflammation')) & e.database.eq('KEGG') & (e.q_value<.05)]
    c=pd.read_csv(ctd,comment='#',header=None,names=['DiseaseName','DiseaseID','PathwayName','PathwayID','Genes'])
    c=c[c.DiseaseID.eq('MESH:D065626') & c.PathwayID.str.startswith('KEGG:')]
    support={}
    for (pid,name),g in c.groupby(['PathwayID','PathwayName']):
        genes={v for text in g.Genes.dropna() for v in re.split(r'[|;]',text)}-{''}
        key=norm(name)
        if key in support:raise ValueError('Ambiguous CTD pathway-name mapping.')
        support[key]=(pid,len(genes))
    available={n['node_id'] for n in nodes};links=[]
    patterns={'C1':'early increase','C2':'late increase','C3':'early rise then plateau','C4':'mid decrease','C5':'late decrease'}
    for term,g in e.groupby('term'):
        node='KEGG:'+term;ctd_row=support.get(norm(term[len('KEGG_'):]))
        if node not in available or not ctd_row or ctd_row[1]<2:continue
        best=g.sort_values('q_value').iloc[0]
        clusters=sorted({x.split('_')[0] for x in g.group if x[:2] in patterns})
        links.append({'pathway':node,'disease':'EFO_0003095','props':{
            'source':'task4_enrichment+CTD','database':'KEGG','task4_best_q':float(best.q_value),
            'task4_best_fold':round(float(best.fold_enrichment),3),'task4_best_group':best.group,
            'task4_groups':';'.join(sorted(g.group.unique())),'clusters':';'.join(clusters),
            'cluster_patterns':';'.join(patterns[x] for x in clusters),'ctd_pathway_id':ctd_row[0],'ctd_n_genes':ctd_row[1]}})
    if len(links)!=61:raise ValueError(f'CTD/enrichment inputs reproduce {len(links)} KEGG disease links; expected 61. Use the matching source snapshot.')
    return {'nodes':nodes,'members':members,'links':links}


def apply(tx, data):
    tx.run('UNWIND $rows AS row MERGE(p:Pathway {node_id:row.node_id}) SET p += row',rows=data['nodes']).consume()
    for i in range(0,len(data['members']),2000):
        tx.run('UNWIND $rows AS row MATCH(g:Gene {symbol:row.symbol}) MATCH(p:Pathway {node_id:row.pathway}) MERGE(g)-[r:IN_PATHWAY]->(p) SET r.source="KEGG_legacy"',rows=data['members'][i:i+2000]).consume()
    tx.run('UNWIND $rows AS row MATCH(p:Pathway {node_id:row.pathway}) MATCH(d:Disease {node_id:row.disease}) MERGE(p)-[r:IMPLICATED_IN {source:"task4_enrichment+CTD"}]->(d) SET r += row.props',rows=data['links']).consume()


def verify(tx, data):
    expected={n['node_id']:n for n in data['nodes']}
    actual={r['id']:r['props'] for r in tx.run('MATCH(p:Pathway {source:"KEGG_legacy"}) RETURN p.node_id AS id,properties(p) AS props')}
    if actual!=expected:raise ValueError('KEGG node/property verification failed.')
    wanted={(r['symbol'],r['pathway']) for r in data['members']}
    found=[(r['symbol'],r['pathway']) for r in tx.run('MATCH(g:Gene)-[r:IN_PATHWAY]->(p:Pathway {source:"KEGG_legacy"}) RETURN g.symbol AS symbol,p.node_id AS pathway,r.source AS source') if r['source']=='KEGG_legacy']
    if len(found)!=len(wanted) or set(found)!=wanted:raise ValueError('KEGG membership verification failed.')
    actual=[{'pathway':r['pathway'],'disease':r['disease'],'props':r['props']} for r in tx.run('MATCH(p:Pathway {source:"KEGG_legacy"})-[r:IMPLICATED_IN]->(d:Disease) RETURN p.node_id AS pathway,d.node_id AS disease,properties(r) AS props')]
    if sorted(actual,key=lambda r:r['pathway'])!=sorted(data['links'],key=lambda r:r['pathway']):raise ValueError('KEGG disease-link/property verification failed.')
