#!/usr/bin/env python
"""Quick qualitative ranking of iwdescr descriptors against an activity.

Each descriptor is split into quantile buckets and the activity variation across the
buckets is summarised as eta^2 (fraction of activity variance explained by the bucket
means). A random forest gives a second, multivariate, view (impurity importance, with
an out-of-bag R^2 to say whether the forest learned anything). The two are combined by
averaging their percentile ranks.

Correlated descriptors (Spearman |rho| >= --cluster_rho) are clustered. Only the best
descriptor of each cluster is plotted, with the others listed beneath it. A collapsed
member that still explains activity variation the representative cannot (residual eta^2)
is promoted to its own plot. Everything is written to one self-contained HTML file, and
optionally a CSV with every descriptor.

  descriptor_ranking.py -A activity.txt -o report.html descriptors.txt
"""
import argparse
import base64
import html
import io
import os
import sys
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    from scipy.stats import rankdata
    from sklearn.ensemble import RandomForestRegressor
except ImportError as exc:
    raise SystemExit(f'Missing dependency: {exc}. Install numpy pandas scipy scikit-learn matplotlib') from None


def read_table(path):
    """First column is the identifier. .csv/.tsv by suffix, otherwise whitespace separated."""
    path = Path(path)
    sep = {'.csv': ',', '.tsv': '\t'}.get(path.suffix, r'\s+')
    first = pd.read_csv(path, sep=sep, nrows=0, comment='#').columns[0]
    df = pd.read_csv(path, sep=sep, comment='#', dtype={first: str}, index_col=0)
    if df.index.has_duplicates:
        raise SystemExit(f'{path}: duplicate identifiers, e.g. {df.index[df.index.duplicated()][0]!r}')
    return df.apply(pd.to_numeric, errors='coerce')   # '.' and friends become NaN


def load(args):
    act = read_table(args.activity)
    column = args.activity_column or act.columns[0]
    if column not in act.columns:
        raise SystemExit(f'activity column {column!r} not in {args.activity}')
    y = act[column].dropna()
    desc = read_table(args.descriptors)
    ids = y.index.intersection(desc.index)
    if len(ids) < 20:
        raise SystemExit(f'only {len(ids)} molecules in common between the two files')
    y, desc = y.loc[ids], desc.loc[ids]
    if args.transform != 'none':
        if (y <= 0).any():
            raise SystemExit('log transforms need strictly positive activities')
        y = np.log10(y) * (-1 if args.transform == 'neglog10' else 1)
    if y.nunique() < 2:
        raise SystemExit('activity does not vary')
    return y, desc, column, len(desc.index.difference(ids)) + len(act.index.difference(ids))


def filter_descriptors(desc, args):
    keep, dropped = [], []
    for name in desc.columns:
        x = desc[name]
        if x.isna().mean() > args.max_missing:
            dropped.append((name, 'missing'))
        elif x.nunique() < 2 or x.value_counts(normalize=True).iloc[0] >= args.max_mode:
            dropped.append((name, 'constant / dominant value'))
        else:
            keep.append(name)
    return desc[keep].fillna(desc[keep].median()), dropped


def bucket_codes(x, k):
    """Quantile buckets; tied values are never split, so there can be fewer than k."""
    edges = np.unique(np.quantile(x, np.linspace(0, 1, k + 1)[1:-1]))
    return np.unique(np.searchsorted(edges, x, side='left'), return_inverse=True)[1]


def eta_squared(codes, y):
    n = np.bincount(codes)
    sums = np.bincount(codes, weights=y)
    centred = y - y.mean()
    total = centred @ centred
    return float((np.sum(sums ** 2 / n) - len(y) * y.mean() ** 2) / total) if total > 0 else 0.0


def bucket_means(codes, y):
    return np.bincount(codes, weights=y) / np.bincount(codes)


def cluster_descriptors(ranks, names, rho_cut):
    """Average-linkage clustering on 1-|Spearman|; returns {name: cluster id} and the rho matrix."""
    rho = np.nan_to_num(np.corrcoef(ranks, rowvar=False))
    if len(names) < 2:
        return {names[0]: 1}, rho
    dist = 1 - np.abs(rho)
    np.fill_diagonal(dist, 0)
    labels = fcluster(linkage(squareform(dist, checks=False), method='average'), 1 - rho_cut, criterion='distance')
    return dict(zip(names, labels)), rho


def iwdescr_descriptions():
    roots = [Path(__file__).resolve().parents[2]]
    if os.environ.get('LILLYMOL_HOME'):
        roots.insert(0, Path(os.environ['LILLYMOL_HOME']))
    for root in roots:
        path = root / 'docs/Molecule_Tools/iwdescr.md'
        if path.is_file():
            out = {}
            for line in path.read_text().splitlines():
                f = [s.strip() for s in line.split('|')]
                if line.startswith('|') and len(f) == 5 and f[1] not in ('Name', '----', ''):
                    out[f[1]] = f[2]
            return out
    return {}


def fmt(v):
    return f'{v:.3g}'


def plot_buckets(name, x, y, codes, title, rng, max_points=400):
    k = codes.max() + 1
    mean, count = bucket_means(codes, y), np.bincount(codes)
    sd = np.sqrt(np.maximum(np.bincount(codes, weights=y ** 2) / count - mean ** 2, 0) * count / np.maximum(count - 1, 1))
    ci = 1.96 * sd / np.sqrt(count)
    fig, ax = plt.subplots(figsize=(6.4, 3.3), dpi=100)
    pick = rng.choice(len(y), min(max_points, len(y)), replace=False)
    ax.scatter(codes[pick] + rng.uniform(-0.3, 0.3, len(pick)), y[pick], s=5, color='#9aa7b4', alpha=0.35, linewidths=0)
    ax.axhline(y.mean(), color='#888', ls='--', lw=1)
    ax.errorbar(range(k), mean, yerr=ci, fmt='o-', color='#b43e25', ecolor='#b43e25', capsize=3, lw=1.5, ms=5)
    lo = [x[codes == b].min() for b in range(k)]
    hi = [x[codes == b].max() for b in range(k)]
    ax.set_xticks(range(k))
    ax.set_xticklabels([f'{fmt(a)}–{fmt(b)}\n n={n}' if a != b else f'{fmt(a)}\n n={n}'
                        for a, b, n in zip(lo, hi, count)], fontsize=6.5)
    ax.set_xlabel(name, fontsize=9)
    ax.set_ylabel('activity (bucket mean, 95% CI)', fontsize=8)
    ax.set_title(title, fontsize=9, loc='left')
    ax.tick_params(axis='y', labelsize=8)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def histogram(y, name):
    fig, ax = plt.subplots(figsize=(4, 2.2), dpi=100)
    ax.hist(y, bins=30, color='#647e9e')
    ax.set_xlabel(name, fontsize=8)
    ax.tick_params(labelsize=7)
    for side in ('top', 'right', 'left'):
        ax.spines[side].set_visible(False)
    ax.set_yticks([])
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def analyse(args):
    y, desc, activity_name, unmatched = load(args)
    X, dropped = filter_descriptors(desc, args)
    names = list(X.columns)
    if not names:
        raise SystemExit('no descriptors survived filtering')
    yv, Xv = y.to_numpy(float), X.to_numpy(float)
    n = len(yv)

    # Bucket eta^2 per descriptor.
    codes = {name: bucket_codes(Xv[:, i], args.buckets) for i, name in enumerate(names)}
    eta = np.array([eta_squared(codes[name], yv) for name in names])

    # One random forest; impurity importance is free after the fit.
    rf = RandomForestRegressor(n_estimators=args.trees, min_samples_leaf=3, max_features=1 / 3,
                               oob_score=True, n_jobs=args.jobs, random_state=args.seed)
    rf.fit(Xv, yv)
    imp = rf.feature_importances_
    oob = float(rf.oob_score_)

    ranks = rankdata(Xv, axis=0)
    trend = np.array([np.corrcoef(ranks[:, i], rankdata(yv))[0, 1] for i in range(len(names))])
    table = pd.DataFrame(dict(eta2=eta, rf_importance=imp, spearman_with_activity=trend), index=names)
    table['score'] = (table.eta2.rank(pct=True) + table.rf_importance.rank(pct=True)) / 2

    cluster_of, rho = cluster_descriptors(ranks, names, args.cluster_rho)
    table['cluster'] = [cluster_of[nm] for nm in names]
    position = {nm: i for i, nm in enumerate(names)}
    prefer = [p for p in args.prefer.split(',') if p]

    chance = max(args.residual_eta, 3 * (args.buckets - 1) / (n - 1))
    table['representative'] = False
    table['promoted'] = False
    table['residual_eta2'] = np.nan
    items = []     # one plotted item per cluster or promoted member
    for cid, group in table.groupby('cluster'):
        best = group.score.max()
        near = group[group.score >= best - 0.05].index
        rep = next((p for p in prefer if p in near), group.score.idxmax() if len(near) == 0 else
                   group.loc[near].score.idxmax())
        table.loc[rep, 'representative'] = True
        members = [m for m in group.index if m != rep]
        # Remove the representative's effect with a finer fit than the display buckets, or
        # within-bucket trend of the representative leaks into its correlated neighbours.
        fine = bucket_codes(Xv[:, position[rep]], max(args.buckets, min(4 * args.buckets, n // 30)))
        residual = yv - bucket_means(fine, yv)[fine]
        promoted = []
        for m in members:
            # Partial association: remove the representative's effect from both the member
            # and the activity, otherwise the shared part swamps the member's own signal.
            xm = Xv[:, position[m]]
            xres = xm - bucket_means(fine, xm)[fine]
            r = eta_squared(bucket_codes(xres, args.buckets), residual) if np.ptp(xres) > 0 else 0.0
            table.loc[m, 'residual_eta2'] = r
            if r > chance:
                table.loc[m, 'promoted'] = True
                promoted.append(m)
        collapsed = [m for m in members if m not in promoted]
        items.append(dict(name=rep, members=group.index.tolist(), collapsed=collapsed,
                          eta=group.eta2.max(), rf=group.rf_importance.sum(), adds_to=None))
        for m in promoted:
            items.append(dict(name=m, members=[m], collapsed=[], eta=table.loc[m, 'eta2'],
                              rf=table.loc[m, 'rf_importance'], adds_to=rep))
    cards = pd.DataFrame(items)
    cards['adds_to'] = cards['adds_to'].fillna('')    # NaN is truthy; keep it a plain string
    cards['score'] = (cards.eta.rank(pct=True) + cards.rf.rank(pct=True)) / 2
    cards = cards.sort_values(['score', 'eta'], ascending=False).reset_index(drop=True)
    cards['eta_rank'] = cards.eta.rank(ascending=False, method='min').astype(int)
    cards['rf_rank'] = cards.rf.rank(ascending=False, method='min').astype(int)
    return dict(y=yv, X=Xv, names=names, table=table, cards=cards, codes=codes, rho=rho, position=position,
                oob=oob, n=n, activity=activity_name, dropped=dropped, unmatched=unmatched, chance=chance,
                total_descriptors=desc.shape[1])


CSS = '''body{font:14px system-ui,sans-serif;color:#233144;background:#f5f7fa;margin:24px;max-width:1100px}
h1{font-size:24px}h2{font-size:18px;margin-top:28px}table{border-collapse:collapse;background:#fff}
th,td{padding:5px 9px;border-bottom:1px solid #dde3ea;text-align:left;vertical-align:top}th{background:#e9eef5}
td.n,th.n{text-align:right}.card{background:#fff;border:1px solid #dde3ea;border-radius:6px;padding:12px;margin:14px 0}
.card h3{margin:0 0 4px 0;font-size:16px}.meta{color:#536170;margin:2px 0}.also{margin:6px 0 0 0;font-size:13px}
.warn{color:#95451a;font-weight:600}img{max-width:100%}.flex{display:flex;gap:24px;flex-wrap:wrap;align-items:center}
details{margin:12px 0}summary{cursor:pointer}code{font-size:12px}'''


TIPS = {
    'cluster size': 'Number of correlated descriptors (|Spearman| above the cluster threshold) represented by this plot.',
    'eta': 'eta-squared: fraction of the activity variance explained by the bucket means (0 to 1). '
           'Also detects non-monotonic effects. For a cluster, the best value among its members. '
           '(#k) is its rank among all plots.',
    'rf': 'Random forest importance: this descriptor\'s share of the forest\'s split credit. Importances of all descriptors '
          'sum to 1, so 0.10 means 10%. Summed over the cluster, because correlated descriptors share credit. '
          '(#k) is its rank among all plots.',
    'rho': 'Spearman rank correlation with activity. Positive: higher descriptor values go with higher activity. '
           'Negative: the reverse. Only measures monotonic trend, so a U-shaped effect can give a value near 0.',
    'score': 'Mean of the eta-squared and RF percentile ranks across all plots (1.0 = best on both). Used to order the report.',
    'resid': 'eta-squared of the activity left over after removing the representative descriptor\'s effect, '
             'against this descriptor with the representative\'s effect also removed. Large values mean it adds information.',
}


def th(label, key, numeric=True):
    return f'<th{" class=\"n\"" if numeric else ""} title="{html.escape(TIPS[key], quote=True)}">{label}</th>'


def render(res, args, descriptions):
    table, cards, rho, position = res['table'], res['cards'], res['rho'], res['position']
    rng = np.random.default_rng(args.seed)
    top = cards.head(args.top)
    out = [f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
           f'<title>{html.escape(args.title)}</title><style>{CSS}</style></head><body><h1>{html.escape(args.title)}</h1>']
    ncluster = table.cluster.nunique()
    out.append(f'<div class="flex"><div><p>{res["n"]} molecules; activity <b>{html.escape(res["activity"])}</b>'
               f'{"" if args.transform == "none" else " (" + args.transform + ")"}.<br>'
               f'{len(res["names"])} of {res["total_descriptors"]} descriptors analysed '
               f'({len(res["dropped"])} dropped as constant or mostly missing).<br>'
               f'{ncluster} clusters at |Spearman| &ge; {args.cluster_rho}; showing the top {len(top)} of {len(cards)} plots.<br>'
               f'Random forest out-of-bag R&sup2;: <b>{res["oob"]:.2f}</b></p>')
    if res['oob'] < 0.05:
        out.append('<p class="warn">The forest explains almost none of the activity out of bag, so treat the '
                   'importances (and this ranking) as weak evidence.</p>')
    out.append(f'</div><img src="data:image/png;base64,{histogram(res["y"], res["activity"])}" alt="activity distribution"></div>')
    out.append('<details><summary>How to read this report</summary>'
               f'<p>Each descriptor is cut into {args.buckets} quantile buckets (ties are never split, so there may be fewer). '
               '<b>&eta;&sup2;</b> is the fraction of activity variance explained by the bucket means; unlike a correlation it also '
               'detects non-monotonic effects. <b>RF</b> is the random forest impurity importance (summed over a cluster, since '
               'correlated descriptors share credit). The score is the mean of the two percentile ranks across clusters. '
               'It is for deciding what to look at, not a measure of causation. Impurity importance favours descriptors with many '
               'distinct values. &rho; is the Spearman correlation with activity (sign shows direction).</p>'
               f'<p>Plots show bucket means with 95% intervals (red) over a sample of the molecules (grey). '
               f'Correlated descriptors are listed under the plot of the best one. A collapsed descriptor that still explains '
               f'activity the representative cannot (residual &eta;&sup2; &gt; {res["chance"]:.3f}) gets its own plot marked '
               f'&ldquo;adds to&rdquo;.</p></details>')

    out.append('<h2>Ranking</h2><p class="meta">Hover over a column heading for its meaning.</p><table><tr><th class="n">#</th><th>Descriptor</th>'
               + th('cluster size', 'cluster size') + th('&eta;&sup2; (rank)', 'eta') + th('RF (rank)', 'rf')
               + th('&rho;', 'rho') + th('score', 'score') + '</tr>')
    for i, c in cards.iterrows():
        label = html.escape(c['name']) + (f' <small>(adds to {html.escape(c["adds_to"])})</small>' if c['adds_to'] else '')
        out.append(f'<tr><td class="n">{i + 1}</td><td><a href="#c{i}">{label}</a></td><td class="n">{len(c["members"])}</td>'
                   f'<td class="n">{c.eta:.3f} <small>(#{c.eta_rank})</small></td><td class="n">{c.rf:.3f} <small>(#{c.rf_rank})</small></td>'
                   f'<td class="n">{table.loc[c["name"], "spearman_with_activity"]:+.2f}</td><td class="n">{c.score:.2f}</td></tr>')
        if i + 1 == args.top and len(cards) > args.top:
            out.append('<tr><td colspan="7"><i>&mdash; not plotted below this line &mdash;</i></td></tr>')
    out.append('</table><h2>Plots</h2>')

    for i, c in top.iterrows():
        nm = c['name']
        r = table.loc[nm]
        j = position[nm]
        title = (f'#{i + 1}  {nm}    η²={c.eta:.3f} (#{c.eta_rank})   RF={c.rf:.3f} (#{c.rf_rank})'
                 f'   ρ={r.spearman_with_activity:+.2f}')
        png = plot_buckets(nm, res['X'][:, j], res['y'], res['codes'][nm], title, rng)
        desc = descriptions.get(nm) or descriptions.get(nm.removeprefix('w_'), '')
        out.append(f'<div class="card" id="c{i}"><h3>{i + 1}. {html.escape(nm)}'
                   + (f' <small>&mdash; adds to {html.escape(c["adds_to"])} (residual &eta;&sup2; {r.residual_eta2:.3f})</small>' if c['adds_to'] else '')
                   + '</h3>')
        if desc:
            out.append(f'<p class="meta">{html.escape(desc)}</p>')
        out.append(f'<img src="data:image/png;base64,{png}" alt="{html.escape(nm)} bucket plot">')
        if c['collapsed']:
            others = sorted(c['collapsed'], key=lambda m: -abs(rho[j, position[m]]))
            out.append('<p class="also"><b>also:</b> ' + ', '.join(
                f'{html.escape(m)} ({rho[j, position[m]]:+.2f})' for m in others) + '</p>')
        out.append('</div>')

    out.append('<details><summary>All descriptors</summary><table><tr><th>Descriptor</th><th>cluster</th><th></th>'
               + th('&eta;&sup2;', 'eta') + th('RF', 'rf') + th('&rho;', 'rho') + th('residual &eta;&sup2;', 'resid') + '</tr>')
    for nm, r in table.sort_values(['cluster', 'score'], ascending=[True, False]).iterrows():
        role = 'representative' if r.representative else 'promoted' if r.promoted else ''
        res_eta = '' if np.isnan(r.residual_eta2) else f'{r.residual_eta2:.3f}'
        out.append(f'<tr><td>{html.escape(nm)}</td><td>{r.cluster}</td><td>{role}</td><td class="n">{r.eta2:.3f}</td>'
                   f'<td class="n">{r.rf_importance:.4f}</td><td class="n">{r.spearman_with_activity:+.2f}</td>'
                   f'<td class="n">{res_eta}</td></tr>')
    out.append('</table></details>')
    if res['dropped']:
        out.append('<details><summary>Dropped descriptors (' + str(len(res['dropped'])) + ')</summary><p>'
                   + ', '.join(f'{html.escape(a)} ({b})' for a, b in res['dropped']) + '</p></details>')
    out.append('</body></html>')
    return '\n'.join(out)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('descriptors', type=Path, help='iwdescr output (header; first column is the identifier)')
    p.add_argument('-A', '--activity', type=Path, required=True, help='identifier and activity columns, with header')
    p.add_argument('--activity_column', help='activity column name (default: the second column)')
    p.add_argument('--transform', choices=('none', 'log10', 'neglog10'), default='none')
    p.add_argument('-o', '--output', type=Path, required=True, help='HTML report')
    p.add_argument('--csv', type=Path, help='also write every descriptor and its scores')
    p.add_argument('--buckets', type=int, default=8)
    p.add_argument('--trees', type=int, default=300)
    p.add_argument('--top', type=int, default=25, help='number of plots to draw')
    p.add_argument('--cluster_rho', type=float, default=0.9, help='|Spearman| at which descriptors are clustered')
    p.add_argument('--residual_eta', type=float, default=0.05,
                   help='minimum residual eta^2 for a clustered descriptor to get its own plot (raised automatically for small n)')
    p.add_argument('--prefer', default='', help='comma separated descriptors to favour as cluster representatives when scores are close')
    p.add_argument('--max_missing', type=float, default=0.2)
    p.add_argument('--max_mode', type=float, default=0.99, help='drop descriptors where one value has this fraction')
    p.add_argument('--jobs', type=int, default=-1)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--title', default='Descriptor ranking')
    args = p.parse_args(argv)
    if args.buckets < 2:
        p.error('--buckets must be at least 2')
    res = analyse(args)
    args.output.write_text(render(res, args, iwdescr_descriptions()))
    if args.csv:
        res['table'].sort_values('score', ascending=False).rename_axis('descriptor').to_csv(args.csv, float_format='%.5g')
    print(f'{len(res["names"])} descriptors, {res["table"].cluster.nunique()} clusters, '
          f'OOB R2 {res["oob"]:.2f}; wrote {args.output}', file=sys.stderr)


if __name__ == '__main__':
    main()
