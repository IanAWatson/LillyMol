"""Molecule-independent HTML rendering and feature ranking for activity reports."""
import heapq
import html


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


def multiplicity_profiles_html(row: dict[str, object]) -> str:
    profiles = row['multiplicity']['profiles']
    parts = ['<div class="count-profiles">']
    for bucket, profile in profiles.items():
        s = profile['summary']
        if s['n'] == 0:
            continue
        parts.append(
            f'<span><b>{html.escape(bucket)}</b> n {s["n"]} med {fmt(s["median"])}</span>'
        )
    parts.append('</div>')
    return ''.join(parts)


def examples_html(examples: list[dict[str, object]]) -> str:
    if not examples:
        return '<span class="muted">No matched examples</span>'
    parts = ['<div class="examples">']
    for example in examples:
        parts.append('<div class="example">')
        parts.append(f'<code>{html.escape(str(example.get("labelled_smiles", example.get("smiles", ""))))}</code>')
        parts.append(
            '<span>' + html.escape(str(example['identifier'])) +
            f' activity {fmt(example["activity"])}' +
            (f' hits {example["nhits"]}' if 'nhits' in example else '') + '</span>'
        )
        parts.append('</div>')
    parts.append('</div>')
    return ''.join(parts)


def build_html(data: dict[str, object], title: str, rows=None) -> str:
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
.name strong { display: block; }
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
.examples { display: grid; grid-template-columns: repeat(3, minmax(160px, 1fr)); gap: 6px; max-width: 640px; }
.example { border: 1px solid #e2e5e9; background: #fff; padding: 6px; }
.example code { display: block; overflow-wrap: anywhere; }
.example span { display: block; color: var(--muted); font-size: 12px; margin-top: 4px; }
.hidden { display: none; }
@media (max-width: 720px) { header, .toolbar, main { padding-left: 12px; padding-right: 12px; } .status { width: 100%; margin-left: 0; } }
</style>
</head>
<body>
""")
    parts.append('<header>')
    parts.append(f'<h1>{html.escape(title)}</h1>')
    parts.append(
        '<p class="summary">' +
        f'{len(rows)} of {metadata["queries"]} {feature_label}, {metadata["molecules_read_with_activity"]} molecules with activity. ' +
        f'Activity file: <code>{html.escape(str(metadata["activity"]))}</code>. ' +
        f'Global activity median {fmt(global_data["summary"]["median"])}.</p>'
    )
    parts.append('</header>')
    parts.append("""<div class="toolbar">
<input id="search" type="search" placeholder="Filter by query name or source" oninput="updateView()">
<button type="button" onclick="sortTable('delta_median')">Sort delta median</button>
<button type="button" onclick="sortTable('matched_n')">Sort support</button>
""")
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
    parts.append("""<div class="table-wrap">
<table id="results">
<thead><tr>
""")
    parts.append('<th>Fragment</th>' if fragments else '<th>Query</th>')
    parts.append("""<th class="sortable" onclick="sortTable('matched_n')">Matched N</th><th>Matched</th><th>Unmatched</th><th class="sortable" onclick="sortTable('delta_median')">Delta Median</th><th class="sortable" onclick="sortTable('delta_mean')">Delta Mean</th>""")
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
        parts.append(
            f'<tr data-search="{html.escape(searchable, quote=True)}" '
            f'data-matched_n="{matched["n"]}" data-delta_median="{row["delta_median"] if row["delta_median"] is not None else -1e99}" '
            f'data-delta_mean="{row["delta_mean"] if row["delta_mean"] is not None else -1e99}" data-multi="{multi}">'
        )
        parts.append('<td class="name">')
        parts.append(f'<strong>{html.escape(row["name"])}</strong>')
        parts.append(f'<code>{html.escape(row["source"])}</code>')
        parts.append('</td>')
        parts.append(f'<td class="num">{matched["n"]}</td>')
        parts.append(f'<td class="num">median {fmt(matched["median"])}<br>mean {fmt(matched["mean"])}</td>')
        parts.append(f'<td class="num">median {fmt(unmatched["median"])}<br>mean {fmt(unmatched["mean"])}</td>')
        parts.append(f'<td class="num">{fmt(row["delta_median"])}</td>')
        parts.append(f'<td class="num">{fmt(row["delta_mean"])}</td>')
        multi_text = f'{multi}' if multi else '<span class="muted">0</span>'
        if multi:
            multi_text += f'<br><span class="muted">max {max_hits}</span>'
        if metadata.get('multiplicity_mode') == 'count':
            multi_text += multiplicity_profiles_html(row)
        if show_multiplicity:
            parts.append(f'<td class="num">{multi_text}</td>')
        parts.append('<td>' + histogram_svg(row['matched']['histogram'], row['unmatched']['histogram']) + '</td>')
        parts.append('<td>' + examples_html(row['examples']) + '</td>')
        parts.append('</tr>\n')
    parts.append('</tbody></table></div></main>')
    parts.append(f'<script>const featureLabel = "{feature_label}";\n')
    parts.append("""
let sortDescending = true;
function updateView() {
  const needle = document.getElementById('search').value.toLowerCase();
  let shown = 0;
  document.querySelectorAll('#results tbody tr').forEach(row => {
    const match = row.dataset.search.includes(needle);
    row.classList.toggle('hidden', !match);
    if (match) shown += 1;
  });
  document.getElementById('status').textContent = shown + ' ' + featureLabel + ' shown';
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
