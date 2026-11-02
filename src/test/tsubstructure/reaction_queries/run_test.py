#!/usr/bin/env python3
"""Exercise reaction query extraction through the shared -q parser."""
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

EXE = str(Path(sys.argv[1]).resolve())


class ReactionQueries(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        self.reactions = self.directory / 'reactions'
        self.reactions.mkdir()
        self.targets = self.directory / 'targets.smi'
        self.targets.write_text('CCO alcohol\nCN amine\nCF fluoro\nCCCC alkane\n')
        self.reaction = self.reactions / 'reaction.rxn'

    def run_tool(self, query, success=True, **kwargs):
        result = subprocess.run([EXE, '-q', query, '-m', '-', '-m', 'QDT', str(self.targets)],
                                cwd=self.directory, capture_output=True, text=True,
                                timeout=30, **kwargs)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('AddressSanitizer', result.stderr)
        self.assertNotIn('runtime error:', result.stderr)
        return result

    def selected(self, selector=''):
        return self.run_tool(f'RXN:{self.reaction}{selector}').stdout

    def ids(self, text):
        return {line.split()[1] for line in text.splitlines()}

    def count(self, text):
        counts = re.findall(r'\(([0-9]+) matches to ', text)
        self.assertEqual(len(counts), 1, text)
        return int(counts[0])

    def test_scaffold_and_zero_based_sidechains(self):
        self.reaction.write_text('''# Leading comments are supported.
name: "example"
scaffold { query { name: "oxygen" query { smarts: "O" } } }
sidechain { id: 1 smarts: "N" }
sidechain { id: 2 smarts: "F" reagent: "missing_reagent.smi" }
''')
        self.assertEqual(self.ids(self.selected()), {'alcohol'})
        self.assertIn("'oxygen'", self.selected())
        self.assertEqual(self.ids(self.selected(',sidechain=0')), {'amine'})
        self.assertEqual(self.ids(self.selected(',sidechain=1')), {'fluoro'})

    def test_all_site_query_sources(self):
        (self.reactions / 'fluoro.textproto').write_text('# Query file\nquery { smarts: "F" }\n')
        self.reaction.write_text('''scaffold {
  query { query { smarts: "O" } }
  query { query { smarts: "N" } }
  query_file: "fluoro.textproto"
  smarts: "Cl"
}
''')
        self.assertEqual(self.ids(self.selected()), {'alcohol', 'amine', 'fluoro'})

    def test_composite_query_stays_composite(self):
        self.reaction.write_text('''scaffold { query {
  query { smarts: "C" }
  query { smarts: "O" }
  logexp: SS_AND
} }
''')
        self.assertEqual(self.ids(self.selected()), {'alcohol'})

    def test_relative_resources_and_environment_expansion(self):
        (self.reactions / 'groups.textproto').write_text('smarts: "O"\n')
        component = '''query { smarts: "[CD1]" environment {
 attachment { attachment_point: 0 btype: SS_SINGLE_BOND }
 smarts_file: "groups.textproto"
} }'''
        self.targets.write_text('CO alcohol\nCN amine\n')
        self.reaction.write_text(f'scaffold {{ query {{ {component} }} }}')
        self.assertEqual(self.ids(self.selected()), {'alcohol'})
        subdir = self.reactions / 'queries'
        subdir.mkdir()
        (subdir / 'groups.textproto').write_text('smarts: "N"\n')
        (subdir / 'query.textproto').write_text(component)
        self.reaction.write_text('scaffold { query_file: "${RXN_QUERY_DIR}/query.textproto" }')
        env = dict(os.environ, RXN_QUERY_DIR='queries', RXN_REACTION=str(self.reaction))
        result = self.run_tool('RXN:${RXN_REACTION}', env=env)
        self.assertEqual(self.ids(result.stdout), {'amine'})
        # Explicit current-directory lookup takes precedence when requested.
        (self.directory / 'current.textproto').write_text('query { smarts: "O" }')
        self.reaction.write_text('query_files_in_current_directory: true\nscaffold { query_file: "current.textproto" }')
        self.assertEqual(self.ids(self.selected()), {'alcohol'})

    def test_explicit_match_conditions_override_queries(self):
        self.targets.write_text('CCCC alkane\n')
        # Reaction field, corresponding query field, restricted and unrestricted values.
        cases = [('find_unique_embeddings', 'unique_embeddings_only', 'true', 'false'),
                 ('one_embedding_per_start_atom', 'one_embedding_per_start_atom', 'true', 'false'),
                 ('ignore_symmetry_related_matches', 'perceive_symmetric_equivalents', 'true', 'false'),
                 ('embeddings_can_overlap', 'embeddings_do_not_overlap', 'false', 'true'),
                 ('max_matches_to_find', 'max_matches_to_find', '1', '2')]
        for reaction_field, query_field, restricted, unrestricted in cases:
            for desired, original in ((restricted, unrestricted), (unrestricted, restricted)):
                query_desired, query_original = desired, original
                if reaction_field in ('ignore_symmetry_related_matches', 'embeddings_can_overlap'):
                    query_desired = 'false' if desired == 'true' else 'true'
                    query_original = 'false' if original == 'true' else 'true'
                expected_file = self.directory / 'expected.textproto'
                expected_file.write_text(f'query {{ smarts: "CC" {query_field}: {query_desired} }}')
                expected = self.count(self.run_tool(f'PROTO:{expected_file}').stdout)
                original_file = self.directory / 'original.textproto'
                original_file.write_text(f'query {{ smarts: "CC" {query_field}: {query_original} }}')
                self.assertNotEqual(expected, self.count(self.run_tool(f'PROTO:{original_file}').stdout))
                for site, selector in (('scaffold', ''), ('sidechain', ',sidechain=0')):
                    with self.subTest(field=reaction_field, value=desired, site=site):
                        self.reaction.write_text(f'''{site} {{
 query {{ query {{ smarts: "CC" {query_field}: {query_original} }} }}
 match_conditions {{ {reaction_field}: {desired} }}
}}''')
                        self.assertEqual(self.count(self.selected(selector)), expected)
                        # An empty message must preserve the embedded setting.
                        self.reaction.write_text(f'{site} {{ query {{ query {{ smarts: "CC" {query_field}: {query_desired} }} }} match_conditions {{}} }}')
                        self.assertEqual(self.count(self.selected(selector)), expected)

    def test_zero_match_limit_preserves_embedded_limit(self):
        self.targets.write_text('CCCC alkane\n')
        self.reaction.write_text('''scaffold {
 query { query { smarts: "C" max_matches_to_find: 2 } }
 match_conditions { max_matches_to_find: 0 }
}''')
        self.assertEqual(self.count(self.selected()), 2)

    def test_reaction_only_conditions_do_not_change_search(self):
        self.targets.write_text('CCCC alkane\n')
        self.reaction.write_text('scaffold { smarts: "CC" }')
        expected = self.count(self.selected())
        self.reaction.write_text('''scaffold { smarts: "CC" match_conditions {
 process_hit_number: 100
 ignore_not_reacting: true
 suppress_if_more_than_this_many_substructure_search_hits: 1
 enumerate_scaffold_hits_individually: true
} }''')
        self.assertEqual(self.count(self.selected()), expected)

    def test_invalid_sources_and_selectors(self):
        self.reaction.write_text('scaffold { smarts: "O" } sidechain { smarts: "N" }')
        for selector in (',sidechain=-1', ',sidechain=x', ',sidechain=1', ',sidechain=', ',unknown=0'):
            with self.subTest(selector=selector):
                self.run_tool(f'RXN:{self.reaction}{selector}', success=False)
        self.run_tool('RXN:', success=False)
        self.run_tool(f'RXN:{self.reactions / "missing.rxn"}', success=False)
        for text in ('', 'scaffold {}', 'sidechain { smarts: "N" }',
                     'scaffold { query_file: "missing.textproto" }',
                     'scaffold { smarts: "[invalid" }', '$RXN\nlegacy MDL reaction',
                     'scaffold { smarts: "C" match_conditions { max_matches_to_find: 4294967295 } }'):
            with self.subTest(text=text):
                self.reaction.write_text(text)
                self.run_tool(f'RXN:{self.reaction}', success=False)
        self.reaction.write_text('sidechain { smarts: "N" match_conditions { max_matches_to_find: -1 } }')
        self.run_tool(f'RXN:{self.reaction},sidechain=0', success=False)


if __name__ == '__main__':
    result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(ReactionQueries))
    if not result.wasSuccessful():
        sys.exit(1)
    print('Reaction query loading checks passed')
