"""Shared report tests, including optional real RDKit depiction checks."""
import base64
import copy
import html
from html.parser import HTMLParser
import importlib.util
import json
import re
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
import xml.etree.ElementTree as ET

BIN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BIN))
from activity_report import (QUERY_PREVIEW_LENGTH, build_html, depiction_svg,
                             examples_html, query_source_html, ranked_rows)


def row(name, effect, support=2, absent=2):
    return dict(name=name, source=name, delta_median=effect, delta_mean=effect,
                matched=dict(summary=dict(n=support, median=effect, mean=effect), histogram=[0, support]),
                unmatched=dict(summary=dict(n=absent, median=0, mean=0), histogram=[absent, 0]),
                examples=[dict(identifier='<parent>', smiles='C&N', activity=1)])


def report():
    return dict(metadata=dict(feature_kind='fragment', queries=4, activity='activity.txt',
                              molecules_read_with_activity=4),
                **{'global': dict(summary=dict(median=0, mean=0), histogram=[2, 2])},
                queries=[row('positive', 4), row('negative', -8),
                         row('sparse', 10, support=1), row('universal', 99, absent=0)])


class ActivityReportTest(unittest.TestCase):
    def test_ranking_support_and_signed_directions(self):
        data = report()
        self.assertEqual([r['name'] for r in ranked_rows(data, top_n=2)], ['sparse', 'negative'])
        self.assertEqual([r['name'] for r in ranked_rows(data, min_matches=2)], ['negative', 'positive'])
        self.assertEqual([r['name'] for r in ranked_rows(data, sort_by='median', top_n=1)], ['sparse'])
        self.assertEqual([r['name'] for r in ranked_rows(data, sort_by='median', ascending=True, top_n=1)], ['negative'])
        self.assertEqual(len(ranked_rows(data, top_n=0, min_unmatched=0)), 4)
        self.assertEqual(len(data['queries']), 4)

    def test_ties_and_missing_deltas(self):
        data = report()
        data['queries'] = [row('Z', -5), row('A', 5), row('missing', 0)]
        data['queries'][-1]['delta_median'] = None
        self.assertEqual([r['name'] for r in ranked_rows(data)], ['A', 'Z', 'missing'])
        with self.assertRaises(ValueError):
            ranked_rows(data, top_n=-1)

    def test_fragment_and_dictionary_rendering(self):
        data = report()
        html = build_html(data, '<title>', rows=ranked_rows(data, top_n=1))
        self.assertIn('1 of 4 fragments', html)
        self.assertIn('<th>Fragment</th>', html)
        self.assertNotIn('>Multiple</th>', html)
        self.assertIn('&lt;parent&gt;', html)
        self.assertIn('C&amp;N', html)
        self.assertIn('&lt;title&gt;', html)
        dictionary = copy.deepcopy(data)
        dictionary['metadata'].pop('feature_kind')
        for item in dictionary['queries']:
            item['multiplicity'] = dict(molecules_with_multiple_matches=0, max_matches=1)
        html = build_html(dictionary, 'Queries')
        self.assertIn('<th>Query</th>', html)
        self.assertIn('>Multiple</th>', html)

    def test_median_suppression_and_mean_profiles(self):
        data = report()
        data['metadata'].pop('feature_kind')
        data['metadata']['multiplicity_mode'] = 'count'
        for item in data['queries']:
            item['multiplicity'] = dict(molecules_with_multiple_matches=1, max_matches=2,
                profiles={'2': dict(summary=dict(n=2, median=7, mean=8))})
        before = copy.deepcopy(data)
        standard = build_html(data, 'Activity')
        self.assertIn('Delta Median', standard)
        self.assertIn('Global activity median', standard)
        rendered = build_html(data, 'Activity', show_median=False)
        self.assertNotIn('median', rendered.lower())
        self.assertNotIn(' med ', rendered)
        self.assertIn('Global activity mean', rendered)
        self.assertIn('Sort delta mean', rendered)
        self.assertIn('Delta Mean', rendered)
        self.assertIn('n 2 mean 8', rendered)
        self.assertEqual(data, before)

    def test_cli_with_no_site_packages_stdin_and_top_n(self):
        with tempfile.TemporaryDirectory(prefix='activity report ') as tmp:
            root = Path(tmp)
            data = report()
            data['queries'] = [row(f'fragment{i}', i) for i in range(101)]
            data['metadata']['queries'] = 101
            source = root / 'report.json'
            source.write_text(json.dumps(data))
            command = [sys.executable, '-S', str(BIN / 'dicer_fragment_activity_html.py')]
            result = subprocess.run(command + [str(source), '-o', '-'], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('100 of 101 fragments', result.stdout)
            self.assertNotIn('<strong>fragment0</strong>', result.stdout)
            result = subprocess.run(command + ['-', '-o', '-', '--top_n', '0'],
                                    input=source.read_text(), capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('101 of 101 fragments', result.stdout)
            before = source.read_bytes()
            result = subprocess.run(command + [str(source), '-o', str(source)], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(source.read_bytes(), before)
            alias = root / 'hardlink.json'
            alias.hardlink_to(source)
            result = subprocess.run(command + [str(source), '-o', str(alias)], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(source.read_bytes(), before)


class ReportElements(HTMLParser):
    def __init__(self):
        super().__init__()
        self.images = []
        self.searches = []
        self.rows = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == 'img':
            self.images.append(attrs)
        if tag == 'tr' and 'data-search' in attrs:
            self.searches.append(attrs['data-search'])
            self.rows.append({key[5:]: value for key, value in attrs.items() if key.startswith('data-')})


class CompactQueryTest(unittest.TestCase):
    def test_short_source_is_not_truncated(self):
        source = '[N;H2]-[C] & <source>'
        self.assertEqual(query_source_html(source), '<code>' + html.escape(source) + '</code>')
        boundary = 'C' * QUERY_PREVIEW_LENGTH
        self.assertNotIn('<details', query_source_html(boundary))

    def test_long_source_has_escaped_preview_and_full_definition(self):
        source = 'C' * QUERY_PREVIEW_LENGTH + '<script>alert("x")</script>&tail'
        rendered = query_source_html(source)
        self.assertIn('C' * QUERY_PREVIEW_LENGTH + '\u2026</code>', rendered)
        self.assertIn('<summary>Full query</summary>', rendered)
        self.assertIn('<pre>' + html.escape(source) + '</pre>', rendered)
        self.assertNotIn('<script>', rendered)

    def test_full_source_remains_searchable_and_input_is_unchanged(self):
        data = report()
        data['queries'][0]['source'] = 'C' * 120 + ' MixedCaseUniqueTail <&>'
        original = copy.deepcopy(data)
        with mock.patch('activity_report.depiction_svg', return_value=None):
            rendered = build_html(data, 'Long queries')
        elements = ReportElements()
        elements.feed(rendered)
        self.assertIn(data['queries'][0]['source'].lower(), elements.searches[0])
        self.assertEqual(data, original)

    def test_missing_backend_retains_visible_smiles_and_metadata(self):
        depiction_svg.cache_clear()
        self.addCleanup(depiction_svg.cache_clear)
        with mock.patch('activity_report._depiction_modules', return_value=None):
            rendered = examples_html([dict(identifier='<example>', activity=2,
                                          labelled_smiles='[1CH3]N', nhits=1)])
        self.assertIn('Depiction unavailable', rendered)
        self.assertIn('<code>[1CH3]N</code>', rendered)
        self.assertIn('&lt;example&gt; activity 2 hits 1', rendered)
        self.assertNotIn('<img', rendered)

    def test_import_and_ranking_do_not_load_rdkit(self):
        code = ('import sys; sys.path.insert(0, sys.argv[1]); '
                'import activity_report; '
                'activity_report.ranked_rows({"queries": []}); '
                'assert "rdkit" not in sys.modules')
        result = subprocess.run([sys.executable, '-S', '-c', code, str(BIN)],
                                capture_output=True, text=True, timeout=10)
        self.assertEqual(result.returncode, 0, result.stderr)


@unittest.skipUnless(importlib.util.find_spec('rdkit'), 'requires RDKit')
class StructureDepictionTest(unittest.TestCase):
    def setUp(self):
        depiction_svg.cache_clear()
        self.addCleanup(depiction_svg.cache_clear)

    def test_labelled_atoms_are_highlighted_and_marker_isotopes_hidden(self):
        svg = depiction_svg('[1CH3][1CH2]O', True)
        self.assertIsNotNone(svg)
        root = ET.fromstring(svg)
        ellipses = [node.attrib.get('class') for node in root.iter()
                    if node.tag.endswith('ellipse')]
        self.assertEqual(ellipses, ['atom-0', 'atom-1'])
        # With the marker isotopes hidden, these carbons have no atom-label glyphs.
        labels = [node.attrib.get('class') for node in root.iter()
                  if node.tag.endswith('path')]
        self.assertNotIn('atom-0', labels)
        self.assertNotIn('atom-1', labels)
        self.assertIn('atom-2', labels)  # Oxygen remains visible.

    def test_plain_smiles_has_no_synthetic_highlights(self):
        svg = depiction_svg('CCO')
        self.assertIsNotNone(svg)
        root = ET.fromstring(svg)
        self.assertFalse(any(node.tag.endswith('ellipse') for node in root.iter()))

    def test_real_isotopes_are_preserved(self):
        svg = depiction_svg('[1CH3][13CH2]O', True)
        self.assertIsNotNone(svg)
        root = ET.fromstring(svg)
        ellipses = [node.attrib.get('class') for node in root.iter()
                    if node.tag.endswith('ellipse')]
        self.assertEqual(ellipses, ['atom-0'])
        self.assertTrue(any(node.tag.endswith('path') and node.attrib.get('class') == 'atom-1'
                            for node in root.iter()))

    def test_self_contained_image_and_expandable_smiles(self):
        example = dict(identifier='example <&"', activity=1.25, nhits=2,
                       smiles='CCO', labelled_smiles='[1CH3][1CH2]O')
        rendered = examples_html([example])
        elements = ReportElements()
        elements.feed(rendered)
        self.assertEqual(len(elements.images), 1)
        image = elements.images[0]
        self.assertEqual(image['alt'], 'Structure of example <&"')
        prefix = 'data:image/svg+xml;base64,'
        self.assertTrue(image['src'].startswith(prefix))
        ET.fromstring(base64.b64decode(image['src'][len(prefix):]))
        self.assertIn('<summary>SMILES</summary>', rendered)
        self.assertIn('[1CH3][1CH2]O', rendered)
        self.assertIn('activity 1.25 hits 2', rendered)

    def test_invalid_smiles_falls_back_without_breaking_report(self):
        rendered = examples_html([dict(identifier='bad', activity=1,
                                       smiles='not a molecule <&>')])
        self.assertIn('Depiction unavailable', rendered)
        self.assertIn('not a molecule &lt;&amp;&gt;', rendered)
        self.assertNotIn('<img', rendered)

    def test_repeated_examples_reuse_depiction(self):
        first = depiction_svg('CCO')
        second = depiction_svg('CCO')
        self.assertEqual(first, second)
        self.assertEqual(depiction_svg.cache_info().hits, 1)


class PaginationTest(unittest.TestCase):
    def run_browser(self, page_size, actions, count=5):
        if not shutil.which('node'):
            self.skipTest('requires Node.js to execute the generated pagination code')
        data = report()
        names = ['alpha', 'beta', 'alpha_other', 'gamma', 'alpha_last']
        data['queries'] = [row(names[i] if count == 5 else f'row{i}', i) for i in range(count)]
        data['metadata']['queries'] = count
        for item in data['queries']:
            item['examples'] = []
        rendered = build_html(data, 'Pages', page_size=page_size)
        elements = ReportElements()
        elements.feed(rendered)
        script = re.search(r'<script>(.*?)</script>', rendered, re.S).group(1)
        harness = r"""
const vm = require('node:vm');
const input = JSON.parse(require('node:fs').readFileSync(0, 'utf8'));
const rows = input.rows.map((dataset, id) => {
  const row = {dataset, id, hidden: false};
  row.classList = {toggle: (name, value) => {row.hidden = value;}};
  return row;
});
const controls = {};
for (const id of ['search', 'status', 'page-status', 'previous-page', 'next-page', 'pagination']) {
  controls[id] = {value: '', textContent: '', disabled: false, hidden: false};
}
const tbody = {
  querySelectorAll: () => rows,
  appendChild: row => { rows.splice(rows.indexOf(row), 1); rows.push(row); }
};
const document = {
  getElementById: id => controls[id],
  querySelectorAll: () => rows,
  querySelector: () => tbody
};
const context = vm.createContext({document});
vm.runInContext(input.script, context);
function state() {
  return {visible: rows.filter(row => !row.hidden).map(row => row.id),
          status: controls.status.textContent, page: controls['page-status'].textContent,
          previous: controls['previous-page'].disabled, next: controls['next-page'].disabled,
          pagerHidden: controls.pagination.hidden};
}
const result = [state()];
for (const action of input.actions) { vm.runInContext(action, context); result.push(state()); }
process.stdout.write(JSON.stringify(result));
"""
        result = subprocess.run(['node', '-e', harness],
            input=json.dumps(dict(script=script, rows=elements.rows, actions=actions)),
            capture_output=True, text=True, timeout=10)
        self.assertEqual(result.returncode, 0, result.stderr)
        return json.loads(result.stdout)

    def test_page_navigation_and_boundaries(self):
        states = self.run_browser(50, ['changePage(1)', 'changePage(1)', 'changePage(1)',
                                       'changePage(-1)'], count=101)
        self.assertEqual(states[0]['visible'], list(range(50)))
        self.assertEqual(states[0]['page'], 'Page 1 of 3')
        self.assertTrue(states[0]['previous'])
        self.assertFalse(states[0]['pagerHidden'])
        self.assertEqual(states[1]['visible'], list(range(50, 100)))
        self.assertEqual(states[2]['visible'], [100])
        self.assertTrue(states[2]['next'])
        self.assertEqual(states[3]['visible'], [100])
        self.assertEqual(states[4]['visible'], list(range(50, 100)))

    def test_filtering_and_sorting_cover_all_rows_and_reset_page(self):
        states = self.run_browser(2, [
            'changePage(2)',
            'sortTable("delta_mean")',
            'changePage(1)',
            'document.getElementById("search").value = "alpha"; updateView()',
            'changePage(1)',
            'document.getElementById("search").value = "absent"; updateView()',
            'document.getElementById("search").value = ""; updateView()'])
        self.assertEqual(states[1]['visible'], [4])
        self.assertEqual(states[2]['visible'], [4, 3])
        self.assertEqual(states[2]['page'], 'Page 1 of 3')
        self.assertEqual(states[4]['visible'], [4, 2])
        self.assertEqual(states[4]['page'], 'Page 1 of 2')
        self.assertEqual(states[5]['visible'], [0])
        self.assertEqual(states[6]['visible'], [])
        self.assertEqual(states[6]['status'], '0 matching fragments')
        self.assertTrue(states[6]['previous'])
        self.assertTrue(states[6]['next'])
        self.assertTrue(states[6]['pagerHidden'])
        self.assertEqual(states[7]['visible'], [4, 3])

    def test_zero_disables_pagination_and_small_reports_hide_controls(self):
        for size in (0, 50):
            with self.subTest(size=size):
                state = self.run_browser(size, [])[0]
                self.assertEqual(state['visible'], list(range(5)))
                self.assertTrue(state['pagerHidden'])
        with self.assertRaises(ValueError):
            build_html(report(), 'Pages', page_size=-1)

    def test_dicer_cli_page_size_is_independent_of_top_n(self):
        command = [sys.executable, '-S', str(BIN / 'dicer_fragment_activity_html.py'),
                   '-', '-o', '-', '--top_n', '2']
        for size in (1, 0):
            result = subprocess.run(command + ['--page_size', str(size)],
                input=json.dumps(report()), capture_output=True, text=True, timeout=10)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn(f'const pageSize = {size};', result.stdout)
            self.assertIn('2 of 4 fragments', result.stdout)
            self.assertEqual(result.stdout.count('<td class="name">'), 2)
        result = subprocess.run(command + ['--page_size', '-1'], input=json.dumps(report()),
                                capture_output=True, text=True, timeout=10)
        self.assertNotEqual(result.returncode, 0)


if __name__ == '__main__':
    unittest.main()
