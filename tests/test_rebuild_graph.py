"""Rebuild checks; integration requires disposable local graph and private source files."""
import copy
import importlib.util
import os
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
def module(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/f'{name}.py')
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
rebuild=module('rebuild_graph');prepare=module('prepare_graph_release')


class InputTests(unittest.TestCase):
    def test_wrong_release_archive_is_rejected_without_output(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);(p/'bad.zip').write_bytes(b'wrong release')
            with self.assertRaises(ValueError):prepare.extract(p/'bad.zip',p/'out')
            self.assertFalse((p/'out').exists())

    def test_saved_drug_inputs_are_complete(self):
        rows=rebuild.gap_properties(ROOT/'Paper2_Results/Recovery2_ChEMBL/graph_gap_properties.csv')
        self.assertEqual(len(rows),28)
        self.assertEqual(sum('recovered_name' in r['props'] for r in rows),16)
        self.assertEqual(len(rebuild.directions(ROOT/'Paper2_Results/Task7_external_layer/external_layer_links.csv')),297)
        self.assertEqual(len(rebuild.indications(ROOT/'Paper2_Results/Task4/graph_write/indication_pairs_plan.csv')),2087)


@unittest.skipUnless(os.environ.get('TEST_NEO4J_URI')=='bolt://127.0.0.1:17689' and os.environ.get('TEST_KEGG_GMT') and os.environ.get('TEST_CTD'), 'disposable imported graph and local licensed inputs only')
class FullRebuildTests(unittest.TestCase):
    def test_complete_rebuild_rollback_and_repeat(self):
        from neo4j import GraphDatabase
        with GraphDatabase.driver(os.environ['TEST_NEO4J_URI'],auth=None) as d,d.session(database='neo4j') as s:
            self.assertEqual(s.execute_read(rebuild.counts),rebuild.BASE)
            data=s.execute_read(rebuild.make_plan,Path(os.environ['TEST_KEGG_GMT']),Path(os.environ['TEST_CTD']))
            bad=copy.deepcopy(data);bad['kegg']['members'].pop()
            with self.assertRaises(ValueError):s.execute_write(rebuild.apply,bad)
            self.assertEqual(s.execute_read(rebuild.counts),rebuild.BASE)
            self.assertEqual(s.execute_read(rebuild.snapshot),data['before'])
            result=s.execute_write(rebuild.apply,data)
            self.assertTrue(result['verified']);self.assertEqual(s.execute_read(rebuild.counts),rebuild.FINAL)
            after=s.execute_read(rebuild.snapshot)
            repeat=s.execute_read(rebuild.make_plan,Path(os.environ['TEST_KEGG_GMT']),Path(os.environ['TEST_CTD']))
            self.assertTrue(s.execute_read(rebuild.verify,repeat)['verified'])
            self.assertEqual(s.execute_read(rebuild.snapshot),after)
