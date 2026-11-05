"""Opt-in checks against a disposable import of the cleaned public Paper 1 graph."""
import importlib.util
import os
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('drug_writer', ROOT/'Paper2_Results/code/write_new_links_to_graph.py')
writer = importlib.util.module_from_spec(spec); spec.loader.exec_module(writer)
from scripts.unify_drug_schema import snapshot, validate, migrate, needs_migration


class EvidenceComparisonTests(unittest.TestCase):
    def test_migration_provenance_does_not_hide_evidence_differences(self):
        wanted = [{'drug':'d','gene':'g','props':{'sources':'A;B','push':'down','n_sources':2,'layer':writer.LAYER}}]
        actual = [{'drug':'d','gene':'g','props':dict(wanted[0]['props'], added='2026-11-05', original_relationship_type='EXTERNAL_ACTS_ON')}]
        self.assertTrue(writer.layer_matches(actual, wanted))
        actual[0]['props']['sources'] = 'A'
        self.assertFalse(writer.layer_matches(actual, wanted))
        self.assertFalse(writer.layer_matches(wanted+wanted, wanted))

    def test_ambiguous_drug_ids_rejected(self):
        with self.assertRaises(ValueError):
            validate({'nodes':[{'props':{'node_id':'d'}},{'props':{'node_id':'d'}}], 'links':[]})


@unittest.skipUnless(os.environ.get('TEST_NEO4J_URI') == 'bolt://127.0.0.1:17689', 'isolated public-release test database only')
class RealGraphTests(unittest.TestCase):
    def setUp(self):
        from neo4j import GraphDatabase
        self.driver = GraphDatabase.driver(os.environ['TEST_NEO4J_URI'], auth=None)
        self.addCleanup(self.driver.close)
        self.s = self.driver.session(database='neo4j'); self.addCleanup(self.s.close)
        n = self.s.run('MATCH(n) RETURN count(n)').single()[0]
        r = self.s.run('MATCH()-[r]->() RETURN count(r)').single()[0]
        self.assertEqual((n,r), (126294,35140166), 'Requires the disposable, clean public graph import')

    def test_all_curated_links_preserve_baseline_and_undo(self):
        before = self.s.execute_read(snapshot); validate(before)
        plan = self.s.execute_read(writer.load_plan, ROOT/'Paper2_Results/New_drug_gene_links')
        self.assertEqual((len(plan['nodes']),len(plan['links'])),(519,5216))
        result = self.s.execute_write(writer.add_links, before, plan, 'test-only')
        try:
            self.assertTrue(result['verified'])
            after = self.s.execute_read(snapshot)
            self.assertEqual(len(after['links']),20674+5216)
            self.assertEqual(len(after['nodes']),12025+519)
            self.assertFalse(needs_migration(after))
            self.assertTrue(writer.layer_matches(self.s.execute_read(writer.existing_layer),plan['links']))
            rerun = self.s.execute_read(writer.load_plan, ROOT/'Paper2_Results/New_drug_gene_links')
            self.assertEqual(len(rerun['nodes']),0)
            self.assertTrue(writer.layer_matches(self.s.execute_read(writer.existing_layer),rerun['links']))
            with self.assertRaises(ValueError): self.s.execute_write(writer.add_links,after,rerun,'duplicate-test')
            self.assertEqual(self.s.execute_read(snapshot),after)
            # Recreate the legacy schema for ALL saved records, then convert it.
            def legacy_schema(tx):
                tx.run('''MATCH(d:Drug)-[r:ACTS_ON {layer:$layer}]->(g:Gene)
                    CREATE(d)-[old:EXTERNAL_ACTS_ON]->(g) SET old=properties(r) DELETE r''',layer=writer.LAYER).consume()
                tx.run('MATCH(d:Drug {source:$layer}) SET d:ExternalDrug',layer=writer.LAYER).consume()
            self.s.execute_write(legacy_schema)
            legacy=self.s.execute_read(snapshot);validate(legacy)
            converted=self.s.execute_write(migrate,legacy,'full-layer-test')
            self.assertEqual((converted['relationships_converted'],converted['nodes_relabelled']),(5216,519))
            self.assertFalse(needs_migration(self.s.execute_read(snapshot)))
            self.assertTrue(writer.layer_matches(self.s.execute_read(writer.existing_layer),plan['links']))
        finally:
            self.s.execute_write(writer.undo,self.s.execute_read(snapshot))
        self.assertEqual(self.s.execute_read(snapshot),before)

    def test_legacy_conversion_preserves_parallel_evidence_and_repeats(self):
        self.s.run('''CREATE(d:Drug:ExternalDrug:SchemaFixture {node_id:'schema-fixture', flags:['A','B']}),
            (g:Gene:SchemaFixture {ensembl_id:'schema-fixture'})
            CREATE(d)-[:EXTERNAL_ACTS_ON {sources:['A'],action:'inhibitor',value:0.0,flag:false}]->(g),
            (d)-[:EXTERNAL_ACTS_ON {sources:['B'],action:'inhibitor',text:''}]->(g),
            (d)-[:ACTS_ON {sources:['C'],action:'agonist'}]->(g)''').consume()
        try:
            before = self.s.execute_read(snapshot);validate(before)
            result=self.s.execute_write(migrate,before,'test-migration')
            self.assertEqual(result['relationships_converted'],2)
            after=self.s.execute_read(snapshot);self.assertFalse(needs_migration(after))
            self.assertEqual(len(after['links']),len(before['links']))
            self.assertEqual(self.s.execute_write(migrate,after,'repeat')['relationships_converted'],0)
            self.assertEqual(self.s.execute_read(snapshot),after)
            def fail(tx):
                tx.run("MATCH(d:SchemaFixture:Drug) SET d.name='must roll back'").consume()
                raise ValueError('rollback fixture')
            with self.assertRaises(ValueError):self.s.execute_write(fail)
            self.assertEqual(self.s.execute_read(snapshot),after)
        finally:
            self.s.run('MATCH(n:SchemaFixture) DETACH DELETE n').consume()
