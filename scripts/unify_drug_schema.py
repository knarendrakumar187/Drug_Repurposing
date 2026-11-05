"""Unify legacy drug labels/links without merging evidence. Dry run by default.

--apply requires a new JSON backup; all changes and preservation checks run
in one transaction. Stop other graph writers during migration. A repeated
run on an already unified graph makes no changes. Credentials are not saved.
"""
import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path


def snapshot(tx):
    links = tx.run('''MATCH (a)-[r:ACTS_ON|EXTERNAL_ACTS_ON]->(b)
        RETURN elementId(r) AS id, type(r) AS type, elementId(a) AS start,
        elementId(b) AS end, labels(a) AS start_labels, labels(b) AS end_labels,
        properties(r) AS props ORDER BY id''').data()
    nodes = tx.run('''MATCH (d) WHERE d:Drug OR d:ExternalDrug
        RETURN elementId(d) AS id, labels(d) AS labels, properties(d) AS props ORDER BY id''').data()
    return {'links': links, 'nodes': nodes}


def validate(before):
    ids = [n['props'].get('node_id') for n in before['nodes']]
    if any(not v for v in ids) or len(set(ids)) != len(ids):
        raise ValueError('Drug identifiers must be present and unique; resolve aliases first.')
    for r in before['links']:
        if not ({'Drug', 'ExternalDrug'} & set(r['start_labels'])) or 'Gene' not in r['end_labels']:
            raise ValueError('Unexpected drug-link endpoint labels; no changes made.')
        if r['type'] == 'EXTERNAL_ACTS_ON' and r['props'].get('original_relationship_type') not in (None, 'EXTERNAL_ACTS_ON'):
            raise ValueError('Conflicting migration provenance; no changes made.')


def needs_migration(before):
    return any(r['type'] == 'EXTERNAL_ACTS_ON' for r in before['links']) or any(
        'ExternalDrug' in n['labels'] for n in before['nodes'])


def migrate(tx, before, record):
    if snapshot(tx) != before:
        raise ValueError('Graph changed since backup; retry with a fresh backup.')
    legacy = [r for r in before['links'] if r['type'] == 'EXTERNAL_ACTS_ON']
    result = tx.run('''UNWIND $rows AS row
        MATCH (a)-[old:EXTERNAL_ACTS_ON]->(b) WHERE elementId(old) = row.id
        CREATE (a)-[new:ACTS_ON]->(b)
        SET new = properties(old)
        SET new.original_relationship_type = 'EXTERNAL_ACTS_ON',
            new.migration_record_id = coalesce(new.migration_record_id, $record)
        DELETE old RETURN row.id AS old_id, elementId(new) AS new_id''',
        rows=legacy, record=record).data()
    if len(result) != len(legacy):
        raise ValueError('Migration count mismatch; transaction rolled back.')
    tx.run('MATCH (d:ExternalDrug) SET d:Drug REMOVE d:ExternalDrug').consume()
    after = snapshot(tx)
    mapping = {r['old_id']: r['new_id'] for r in result}
    actual = {r['id']: r for r in after['links']}
    if len(actual) != len(before['links']):
        raise ValueError('Relationship count changed; transaction rolled back.')
    for old in before['links']:
        new = actual[mapping.get(old['id'], old['id'])]
        props = dict(old['props'])
        if old['type'] == 'EXTERNAL_ACTS_ON':
            props['original_relationship_type'] = 'EXTERNAL_ACTS_ON'
            props.setdefault('migration_record_id', record)
        if (new['type'], new['start'], new['end'], new['props']) != ('ACTS_ON', old['start'], old['end'], props):
            raise ValueError('Endpoint/property preservation failed; transaction rolled back.')
    nodes = {n['id']: n for n in after['nodes']}
    if len(nodes) != len(before['nodes']):
        raise ValueError('Drug-node count changed; transaction rolled back.')
    for old in before['nodes']:
        new = nodes[old['id']]
        labels = set(old['labels'])
        if 'ExternalDrug' in labels:
            labels = (labels - {'ExternalDrug'}) | {'Drug'}
        if new['props'] != old['props'] or set(new['labels']) != labels:
            raise ValueError('Drug-node preservation failed; transaction rolled back.')
    return {'relationships_converted': len(legacy), 'nodes_relabelled': sum(
        'ExternalDrug' in n['labels'] for n in before['nodes']), 'verified': True}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--apply', action='store_true')
    ap.add_argument('--backup', type=Path, help='new backup path, required with --apply')
    ap.add_argument('--uri', default=os.environ.get('NEO4J_URI', 'bolt://localhost:7687'))
    ap.add_argument('--database', default=os.environ.get('NEO4J_DATABASE', 'neo4j'))
    ap.add_argument('--user', default=os.environ.get('NEO4J_USERNAME', os.environ.get('NEO4J_USER', 'neo4j')))
    args = ap.parse_args()
    if args.apply and not args.backup:
        ap.error('--apply requires --backup to a new file')
    if not os.environ.get('NEO4J_PASSWORD'):
        ap.error('set NEO4J_PASSWORD in the environment')
    from neo4j import GraphDatabase
    with GraphDatabase.driver(args.uri, auth=(args.user, os.environ['NEO4J_PASSWORD'])) as drv:
        with drv.session(database=args.database) as s:
            before = s.execute_read(snapshot)
            validate(before)
            if not needs_migration(before):
                print('Already unified; no changes needed.')
                return
            print(json.dumps({'mode': 'apply' if args.apply else 'dry run',
                'legacy_links': sum(r['type'] == 'EXTERNAL_ACTS_ON' for r in before['links']),
                'legacy_labels': sum('ExternalDrug' in n['labels'] for n in before['nodes'])}))
            if args.apply:
                args.backup.parent.mkdir(parents=True, exist_ok=True)
                with args.backup.open('x', encoding='utf-8') as f:
                    json.dump(before, f, indent=2)
                record = 'drug-schema-' + datetime.now(timezone.utc).isoformat()
                print(json.dumps(s.execute_write(migrate, before, record)))


if __name__ == '__main__':
    main()
