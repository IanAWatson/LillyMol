"""Shared activity-report rendering; RDKit is used only for optional depictions."""
import base64
from functools import lru_cache
import heapq
import html
import re
import sys


SORT_FIELDS = ('abs_median', 'abs_mean', 'median', 'mean', 'support')


def ranked_rows(data, sort_by='abs_median', ascending=False, top_n=100,
                min_matches=1, min_unmatched=1):
    """Select supported features with deterministic SMILES/name tie-breaking."""
    if sort_by not in SORT_FIELDS:
        raise ValueError(f'Unknown sort field: {sort_by}')
    if min(top_n, min_matches, min_unmatched) < 0:
        raise ValueError('Limits must be non-negative')

    def eligible():
        for row in data['queries']:
            if (row['matched']['summary']['n'] >= min_matches and
                    row['unmatched']['summary']['n'] >= min_unmatched):
                yield row

    def key(row):
        if sort_by == 'support':
            value = row['matched']['summary']['n']
        else:
            field = 'delta_median' if 'median' in sort_by else 'delta_mean'
            value = row.get(field)
        if value is None:
            return (1, 0, row['source'], row['name'])
        if sort_by.startswith('abs_'):
            value = abs(value)
        return (0, value if ascending else -value, row['source'], row['name'])

    if top_n:
        return heapq.nsmallest(top_n, eligible(), key=key)
    return sorted(eligible(), key=key)


def sort_value(value) -> float:
    """Numeric sort key for HTML data attributes; missing values sort last when descending."""
    return -1e99 if value is None else value


def fmt(value, digits: int = 3) -> str:
    if value is None:
        return ''
    if isinstance(value, int):
        return str(value)
    return f"{value:.{digits}f}".rstrip('0').rstrip('.')


def simple_histogram_svg(counts: list[int], width: int = 220, height: int = 44) -> str:
    max_count = max(counts + [1])
    nbins = max(len(counts), 1)
    slot = width / nbins
    bar_width = max(1.0, slot * 0.72)
    parts = [f'<svg viewBox="0 0 {width} {height}" class="hist global-hist" aria-hidden="true">']
    parts.append(f'<line x1="0" y1="{height - 1}" x2="{width}" y2="{height - 1}"/>')
    for i, count in enumerate(counts):
        h = (count / max_count) * (height - 6)
        x = i * slot + slot * 0.14
        y = height - h - 1
        parts.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{bar_width:.2f}" height="{h:.2f}"/>')
    parts.append('</svg>')
    return ''.join(parts)


def histogram_svg(matched: list[int], unmatched: list[int], width: int = 220, height: int = 54) -> str:
    max_count = max(matched + unmatched + [1])
    nbins = max(len(matched), 1)
    slot = width / nbins
    bar_width = max(1.0, slot * 0.38)
    parts = [f'<svg viewBox="0 0 {width} {height}" class="hist" aria-hidden="true">']
    parts.append(f'<line x1="0" y1="{height - 1}" x2="{width}" y2="{height - 1}"/>')
    for i, count in enumerate(unmatched):
        h = (count / max_count) * (height - 6)
        x = i * slot + slot * 0.12
        y = height - h - 1
        parts.append(f'<rect class="unmatched" x="{x:.2f}" y="{y:.2f}" width="{bar_width:.2f}" height="{h:.2f}"/>')
    for i, count in enumerate(matched):
        h = (count / max_count) * (height - 6)
        x = i * slot + slot * 0.50
        y = height - h - 1
        parts.append(f'<rect class="matched" x="{x:.2f}" y="{y:.2f}" width="{bar_width:.2f}" height="{h:.2f}"/>')
    parts.append('</svg>')
    return ''.join(parts)


def multiplicity_profiles_html(row: dict[str, object], show_median: bool = True) -> str:
    profiles = row['multiplicity']['profiles']
    parts = ['<div class="count-profiles">']
    for bucket, profile in profiles.items():
        s = profile['summary']
        if s['n'] == 0:
            continue
        statistic = f'med {fmt(s["median"])}' if show_median else f'mean {fmt(s["mean"])}'
        parts.append(f'<span><b>{html.escape(bucket)}</b> n {s["n"]} {statistic}</span>')
    parts.append('</div>')
    return ''.join(parts)



QUERY_PREVIEW_LENGTH = 100


def query_source_html(source: str) -> str:
    """Keep long query definitions compact without losing the full text."""
    if len(source) <= QUERY_PREVIEW_LENGTH:
        return f'<code>{html.escape(source)}</code>'
    preview = html.escape(source[:QUERY_PREVIEW_LENGTH] + '\u2026')
    return (f'<code class="query-preview">{preview}</code>'
            '<details class="query-details"><summary>Full query</summary>'
            f'<pre>{html.escape(source)}</pre></details>')


@lru_cache(maxsize=1)
def _depiction_modules():
    # Analysis and ranking must not depend on RDKit. Import it only when an
    # example is rendered, and report an unavailable backend once per process.
    try:
        from rdkit import Chem, rdBase
        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError:
        print('RDKit is unavailable; displaying SMILES instead of structure images. '
              'Install it with: uv pip install rdkit (or python -m pip install rdkit)',
              file=sys.stderr)
        return None
    return Chem, rdBase, rdDepictor, rdMolDraw2D


@lru_cache(maxsize=512)
def depiction_svg(smiles: str, highlight_matches: bool = False) -> str | None:
    """Draw stored examples; LillyMol has already performed all matching."""
    modules = _depiction_modules()
    if modules is None:
        return None
    Chem, rdBase, rdDepictor, rdMolDraw2D = modules
    try:
        # LillyMol supports some structures RDKit cannot depict. The HTML
        # fallback reports these without flooding stderr with parser messages.
        with rdBase.BlockLogs():
            molecule = Chem.MolFromSmiles(smiles)
        if molecule is None or molecule.GetNumAtoms() == 0:
            return None
        highlighted = []
        if highlight_matches:
            for atom in molecule.GetAtoms():
                if atom.GetIsotope() == 1:
                    highlighted.append(atom.GetIdx())
                    atom.SetIsotope(0)
        atom_set = set(highlighted)
        highlighted_bonds = [bond.GetIdx() for bond in molecule.GetBonds()
                             if bond.GetBeginAtomIdx() in atom_set and
                             bond.GetEndAtomIdx() in atom_set]
        rdDepictor.Compute2DCoords(molecule)
        drawer = rdMolDraw2D.MolDraw2DSVG(220, 160)
        options = drawer.drawOptions()
        options.highlightRadius = 0.4
        options.highlightBondWidthMultiplier = 16
        colour = (0.2, 0.8, 0.4)
        drawer.DrawMolecule(
            molecule, highlightAtoms=highlighted, highlightBonds=highlighted_bonds,
            highlightAtomColors={atom: colour for atom in highlighted},
            highlightBondColors={bond: colour for bond in highlighted_bonds})
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except (ValueError, RuntimeError):
        return None



@lru_cache(maxsize=1)
def _fragment_matching_modules():
    try:
        from lillymol import MolFromSmiles, SubstructureQuery, SubstructureResults
    except ImportError:
        print('LillyMol bindings are unavailable; dicer examples will have no fragment '
              'highlighting. Use the LillyMol run_python.sh wrapper with staged bindings.',
              file=sys.stderr)
        return None
    return MolFromSmiles, SubstructureQuery, SubstructureResults


@lru_cache(maxsize=1)
def _warn_attachment_binding_unavailable():
    print('LillyMol bindings lack the isotope-attachment query keyword; rebuild and '
          'stage the nanobind module to enable dicer fragment highlighting.', file=sys.stderr)


@lru_cache(maxsize=128)
def _fragment_query(smiles: str):
    modules = _fragment_matching_modules()
    if modules is None:
        return None
    MolFromSmiles, SubstructureQuery, _ = modules
    fragment = MolFromSmiles(smiles)
    if fragment is None or not fragment.number_isotopic_atoms():
        return None
    query = SubstructureQuery()
    try:
        # Dicer never breaks a ring bond, so ring membership is conserved: a
        # benzene fragment must not match naphthalene, nor propane cyclohexane.
        ok = query.build_from_molecule(fragment, substituents_only_at_isotopic_atoms=True,
                                       atoms_conserve_ring_membership=True,
                                       non_ring_atoms_become_nrings_0=True)
    except TypeError:
        _warn_attachment_binding_unavailable()
        return None
    if not ok:
        return None
    query.set_find_unique_embeddings_only(True)
    return query


def labelled_fragment_examples(fragment_smiles: str, examples: list[dict]) -> list[dict]:
    """Label all LillyMol embeddings in retained dicer parents without changing JSON."""
    # Isotope SMILES start an atom bracket with digits. This cheap gate avoids
    # loading bindings for fragments that have no attachment labels to interpret.
    if not re.search(r'\[[0-9]', fragment_smiles) or not any(
            not example.get('labelled_smiles') for example in examples):
        return examples
    query = _fragment_query(fragment_smiles)
    if query is None:
        return examples
    MolFromSmiles, _, SubstructureResults = _fragment_matching_modules()
    labelled = []
    for example in examples:
        if example.get('labelled_smiles'):
            labelled.append(example)
            continue
        parent = MolFromSmiles(str(example.get('smiles', '')))
        if parent is None or parent.number_isotopic_atoms():
            # Existing parent isotopes may mean something other than a match.
            labelled.append(example)
            continue
        results = SubstructureResults()
        if not query.substructure_search(parent, results):
            labelled.append(dict(example, highlight_note='Fragment match unavailable'))
            continue
        for embedding in results:
            parent.set_isotopes(embedding, 1)
        labelled.append(dict(example, labelled_smiles=str(parent.aromatic_smiles()), nhits=len(results)))
    return labelled


def examples_html(examples: list[dict[str, object]], fragment_smiles: str | None = None) -> str:
    if not examples:
        return '<span class="muted">No matched examples</span>'
    if fragment_smiles is not None:
        examples = labelled_fragment_examples(fragment_smiles, examples)
    parts = ['<div class="examples">']
    for example in examples:
        parts.append('<div class="example">')
        smiles = str(example.get('labelled_smiles') or example.get('smiles', ''))
        svg = depiction_svg(smiles, bool(example.get('labelled_smiles')))
        if svg is None:
            parts.append('<span class="depiction-unavailable">Depiction unavailable</span>')
            parts.append(f'<code>{html.escape(smiles)}</code>')
        else:
            encoded = base64.b64encode(svg.encode('utf-8')).decode('ascii')
            alt = html.escape(f"Structure of {example['identifier']}", quote=True)
            parts.append(f'<img class="structure" loading="lazy" width="220" height="160" '
                         f'alt="{alt}" src="data:image/svg+xml;base64,{encoded}">')
            parts.append('<details class="smiles"><summary>SMILES</summary>'
                         f'<code>{html.escape(smiles)}</code></details>')
        parts.append(
            '<span>' + html.escape(str(example['identifier'])) +
            f' activity {fmt(example["activity"])}' + '</span>'
        )
        if example.get('highlight_note'):
            parts.append('<span class="highlight-note">' + html.escape(example['highlight_note']) + '</span>')
        parts.append('</div>')
    parts.append('</div>')
    return ''.join(parts)


def build_html(data: dict[str, object], title: str, rows=None, *,
               show_median: bool = True, page_size: int = 50) -> str:
    if not isinstance(page_size, int) or page_size < 0:
        raise ValueError('page_size must be a non-negative integer')
    metadata = data['metadata']
    global_data = data['global']
    if rows is None:
        rows = data['queries']
    fragments = metadata.get('feature_kind') == 'fragment'
    feature_label = 'fragments' if fragments else 'queries'
    show_multiplicity = not fragments
    parts = ["""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
"""]
    parts.append(f'<title>{html.escape(title)}</title>\n')
    parts.append("""<style>
:root { --bg: #f6f7f9; --panel: #fff; --ink: #20242a; --muted: #69707d; --border: #d7dce3; --matched: #246b73; --unmatched: #b7c2cf; --warn: #9b5b00; }
* { box-sizing: border-box; }
body { margin: 0; font-family: Arial, Helvetica, sans-serif; background: var(--bg); color: var(--ink); font-size: 14px; }
header { padding: 18px 24px 12px; background: var(--panel); border-bottom: 1px solid var(--border); }
h1 { margin: 0 0 6px; font-size: 22px; line-height: 1.2; }
.summary { margin: 0; color: var(--muted); line-height: 1.4; }
.toolbar { position: sticky; top: 0; z-index: 20; display: flex; flex-wrap: wrap; gap: 10px; align-items: center; padding: 10px 24px; background: #edf1f4; border-bottom: 1px solid var(--border); }
.toolbar input { width: min(420px, 100%); padding: 7px 9px; border: 1px solid var(--border); border-radius: 4px; font: inherit; }
.toolbar button { border: 1px solid var(--border); background: #fff; border-radius: 4px; padding: 7px 10px; font: inherit; cursor: pointer; }
.status { margin-left: auto; color: var(--muted); }
.pagination { display: flex; align-items: center; gap: 10px; margin-bottom: 10px; color: var(--muted); }
.pagination[hidden] { display: none; }
.pagination button { padding: 6px 10px; cursor: pointer; }
.pagination button:disabled { cursor: default; }
main { padding: 16px 24px 32px; }
.legend { display: flex; gap: 14px; margin-bottom: 12px; color: var(--muted); }
.legend span::before { content: ''; display: inline-block; width: 12px; height: 12px; margin-right: 5px; vertical-align: -1px; }
.legend .m::before { background: var(--matched); }
.legend .u::before { background: var(--unmatched); }
.table-wrap { overflow-x: auto; border: 1px solid var(--border); background: var(--panel); }
table { width: 100%; border-collapse: collapse; min-width: 1180px; }
th, td { border-bottom: 1px solid var(--border); padding: 8px 10px; text-align: left; vertical-align: top; }
th { background: #dfe6eb; color: #3d4651; font-size: 12px; text-transform: uppercase; letter-spacing: 0; white-space: nowrap; }
th.sortable { cursor: pointer; }
td.num { text-align: right; font-variant-numeric: tabular-nums; }
.name { min-width: 180px; width: 240px; max-width: 300px; }
.name strong { display: block; overflow-wrap: anywhere; }
.name pre { white-space: pre-wrap; overflow-wrap: anywhere; max-width: 280px; margin: 6px 0; font-size: 12px; }
details summary { cursor: pointer; color: var(--muted); font-size: 12px; margin-top: 4px; }
.name code { display: block; color: var(--muted); font-size: 12px; margin-top: 3px; overflow-wrap: anywhere; }
.hist { width: 220px; height: 54px; }
.hist line { stroke: #c7cdd4; stroke-width: 1; }
.hist .matched { fill: var(--matched); }
.hist .unmatched { fill: var(--unmatched); }
.global-hist rect { fill: #7a8794; }
.global-panel { display: inline-flex; align-items: center; gap: 10px; margin-bottom: 12px; color: var(--muted); }
.count-profiles { display: grid; gap: 2px; color: var(--muted); font-size: 12px; margin-top: 4px; }
.warn { color: var(--warn); font-weight: 600; }
.muted { color: var(--muted); }
.examples { display: grid; grid-template-columns: repeat(3, 220px); gap: 6px; width: max-content; }
.example { border: 1px solid #e2e5e9; background: #fff; padding: 6px; }
.example code { display: block; overflow-wrap: anywhere; }
.structure { display: block; max-width: 100%; height: auto; }
.example .highlight-note { color: var(--warn); }
.example span { display: block; color: var(--muted); font-size: 12px; margin-top: 4px; }
.hidden { display: none; }
@media (max-width: 720px) { header, .toolbar, main { padding-left: 12px; padding-right: 12px; } .status { width: 100%; margin-left: 0; } }
</style>
</head>
<body>
""")
    statistic = 'median' if show_median else 'mean'
    parts.append('<header>')
    parts.append(f'<h1>{html.escape(title)}</h1>')
    parts.append(
        '<p class="summary">' +
        f'{len(rows)} of {metadata["queries"]} {feature_label}, {metadata["molecules_read_with_activity"]} molecules with activity. ' +
        f'Activity file: <code>{html.escape(str(metadata["activity"]))}</code>. ' +
        f'Global activity {statistic} {fmt(global_data["summary"][statistic])}.</p>'
    )
    parts.append('</header>')
    parts.append("""<div class="toolbar">
<input id="search" type="search" placeholder="Filter by query name or source" oninput="updateView()">
<button type="button" onclick="sortTable('matched_n')">Sort support</button>
""")
    parts.append(f'<button type="button" onclick="sortTable(\'delta_{statistic}\')">Sort delta {statistic}</button>')
    if show_multiplicity:
        parts.append('<button type="button" onclick="sortTable(\'multi\')">Sort multiple hits</button>')
    parts.append("""<span id="status" class="status"></span>
</div>
<main>
""")
    parts.append('<div class="global-panel"><span>Global activity distribution</span>')
    parts.append(simple_histogram_svg(global_data['histogram']))
    parts.append('</div>')
    parts.append('<div class="legend"><span class="m">Matched</span><span class="u">Unmatched</span></div>')
    parts.append("""<nav id="pagination" class="pagination" aria-label="Report pages">
<button id="previous-page" type="button" onclick="changePage(-1)">Previous</button>
<span id="page-status" aria-live="polite"></span>
<button id="next-page" type="button" onclick="changePage(1)">Next</button>
</nav>
<div class="table-wrap">
<table id="results">
<thead><tr>
""")
    parts.append('<th>Fragment</th>' if fragments else '<th>Query</th>')
    parts.append("""<th class="sortable" onclick="sortTable('matched_n')">Matched N</th>""")
    if show_median:
        parts.append('<th>Matched</th><th>Unmatched</th>')
    else:
        parts.append("""<th class="sortable" onclick="sortTable('matched_mean')">Matched Mean</th>"""
                     """<th class="sortable" onclick="sortTable('unmatched_mean')">Unmatched Mean</th>""")
    if show_median:
        parts.append("""<th class="sortable" onclick="sortTable('delta_median')">Delta Median</th>""")
    parts.append("""<th class="sortable" onclick="sortTable('delta_mean')">Delta Mean</th>""")
    if show_multiplicity:
        parts.append('<th class="sortable" onclick="sortTable(\'multi\')">Multiple</th>')
    parts.append("""<th>Distribution</th><th>Examples</th>
</tr></thead>
<tbody>
""")

    for row in rows:
        matched = row['matched']['summary']
        unmatched = row['unmatched']['summary']
        multiplicity = row.get('multiplicity', {})
        multi = multiplicity.get('molecules_with_multiple_matches', 0)
        max_hits = multiplicity.get('max_matches', 0)
        searchable = f"{row['name']} {row['source']}".lower()
        median_attribute = ''
        if show_median:
            value = row['delta_median'] if row['delta_median'] is not None else -1e99
            median_attribute = f'data-delta_median="{value}" '
        parts.append(
            f'<tr data-search="{html.escape(searchable, quote=True)}" '
            f'data-matched_n="{matched["n"]}" {median_attribute}'
            f'data-matched_mean="{sort_value(matched["mean"])}" '
            f'data-unmatched_mean="{sort_value(unmatched["mean"])}" '
            f'data-delta_mean="{row["delta_mean"] if row["delta_mean"] is not None else -1e99}" data-multi="{multi}">'
        )
        parts.append('<td class="name">')
        parts.append(f'<strong>{html.escape(row["name"])}</strong>')
        parts.append(query_source_html(row["source"]))
        parts.append('</td>')
        parts.append(f'<td class="num">{matched["n"]}</td>')
        for summary in (matched, unmatched):
            if show_median:
                parts.append(f'<td class="num">median {fmt(summary["median"])}<br>mean {fmt(summary["mean"])}</td>')
            else:
                parts.append(f'<td class="num">{fmt(summary["mean"])}</td>')
        if show_median:
            parts.append(f'<td class="num">{fmt(row["delta_median"])}</td>')
        parts.append(f'<td class="num">{fmt(row["delta_mean"])}</td>')
        multi_text = f'{multi}' if multi else '<span class="muted">0</span>'
        if multi:
            multi_text += f'<br><span class="muted">max {max_hits}</span>'
        if metadata.get('multiplicity_mode') == 'count':
            multi_text += multiplicity_profiles_html(row, show_median)
        if show_multiplicity:
            parts.append(f'<td class="num">{multi_text}</td>')
        parts.append('<td>' + histogram_svg(row['matched']['histogram'], row['unmatched']['histogram']) + '</td>')
        parts.append('<td>' + examples_html(row['examples'], row['source'] if fragments else None) + '</td>')
        parts.append('</tr>\n')
    parts.append('</tbody></table></div></main>')
    parts.append(f'<script>const featureLabel = "{feature_label}";\nconst pageSize = {page_size};\n')
    parts.append("""
let sortDescending = true;
let currentPage = 0;
function updateView(resetPage = true) {
  if (resetPage) currentPage = 0;
  const needle = document.getElementById('search').value.toLowerCase();
  const rows = Array.from(document.querySelectorAll('#results tbody tr'));
  const matching = rows.filter(row => row.dataset.search.includes(needle));
  const pageCount = pageSize ? Math.ceil(matching.length / pageSize) : (matching.length ? 1 : 0);
  currentPage = Math.max(0, Math.min(currentPage, pageCount - 1));
  const start = pageSize ? currentPage * pageSize : 0;
  const end = pageSize ? Math.min(start + pageSize, matching.length) : matching.length;
  const visible = new Set(matching.slice(start, end));
  rows.forEach(row => row.classList.toggle('hidden', !visible.has(row)));
  document.getElementById('status').textContent = matching.length
      ? (start + 1) + '\u2013' + end + ' of ' + matching.length + ' matching ' + featureLabel
      : '0 matching ' + featureLabel;
  document.getElementById('page-status').textContent =
      'Page ' + (matching.length ? currentPage + 1 : 0) + ' of ' + pageCount;
  document.getElementById('previous-page').disabled = currentPage === 0;
  document.getElementById('next-page').disabled = currentPage + 1 >= pageCount;
  document.getElementById('pagination').hidden = !pageSize || pageCount <= 1;
}
function changePage(offset) {
  currentPage += offset;
  updateView(false);
}
function sortTable(key) {
  const tbody = document.querySelector('#results tbody');
  const rows = Array.from(tbody.querySelectorAll('tr'));
  rows.sort((a, b) => {
    const av = Number(a.dataset[key]);
    const bv = Number(b.dataset[key]);
    return sortDescending ? bv - av : av - bv;
  });
  sortDescending = !sortDescending;
  rows.forEach(row => tbody.appendChild(row));
  updateView();
}
updateView();
</script>
</body>
</html>
""")
    return ''.join(parts)
