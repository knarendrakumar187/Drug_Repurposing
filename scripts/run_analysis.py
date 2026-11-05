"""Run Paper 2 analyses against an explicitly selected, fully rebuilt graph.

Writes outputs in Paper2_Results using the existing analysis scripts. Back up
existing results or run in a separate checkout before deliberately recomputing.
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--uri',required=True);ap.add_argument('--http',required=True);ap.add_argument('--database',required=True)
    ap.add_argument('--stage',choices=['core','metapaths','all'],default='core')
    args=ap.parse_args()
    if not os.environ.get('NEO4J_PASSWORD'):ap.error('Set NEO4J_PASSWORD for the Paper 2 database.')
    env=dict(os.environ,NEO4J_URI=args.uri,NEO4J_HTTP=args.http,NEO4J_DATABASE=args.database,PYTHONIOENCODING='utf-8')
    from neo4j import GraphDatabase
    with GraphDatabase.driver(args.uri,auth=(env.get('NEO4J_USERNAME','neo4j'),env['NEO4J_PASSWORD'])) as d,d.session(database=args.database) as s:
        n=s.run('MATCH(n) RETURN count(n)').single()[0];r=s.run('MATCH()-[r]->() RETURN count(r)').single()[0]
        if (n,r)!=(126999,35160254):raise ValueError('Rebuild and verify the complete Paper 2 graph before running analyses.')
        if s.run('MATCH()-[r:EXTERNAL_ACTS_ON]->() RETURN count(r)').single()[0]:raise ValueError('Convert the legacy drug schema first.')
    files=[]
    if args.stage in ('core','all'):
        files += ['task1_master_gene_table.py','task2_drug_gene_edges.py','task3_druggability.py','task4_drug_coverage.py','direction_consistent_drugs.py']
    if args.stage in ('metapaths','all'):files += ['metapath_pipeline.py','compare_top25.py']
    for name in files:
        print('Running',name,flush=True)
        subprocess.run([sys.executable,str(ROOT/'Paper2_Results/code'/name)],cwd=ROOT,env=env,check=True)


if __name__=='__main__':main()
