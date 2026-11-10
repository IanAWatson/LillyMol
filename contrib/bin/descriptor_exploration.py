#!/usr/bin/env python3
"""Explore measured outcomes against an iwdescr table; write JSON and offline HTML."""
import argparse
import csv
import gzip
import json
import math
from pathlib import Path
import sys

try:
    import numpy as np
    from scipy import stats, __version__ as scipy_version
    from sklearn import __version__ as sklearn_version
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.feature_selection import mutual_info_regression
    from sklearn.impute import SimpleImputer
    from sklearn.inspection import permutation_importance
    from sklearn.linear_model import Ridge
    from sklearn.model_selection import KFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import SplineTransformer, StandardScaler
except ImportError as exc:
    raise SystemExit(f'Missing statistical dependency: {exc}. Install with: '
                     'uv pip install numpy scipy scikit-learn') from None


MISSING = {'', '.', 'NA', 'NaN', 'nan', 'null', 'None'}


def read_table(path):
    """Read a rectangular headered whitespace/CSV/TSV table, preserving string ids."""
    path = Path(path)
    opener = gzip.open if path.suffix == '.gz' else open
    suffix = path.with_suffix('').suffix if path.suffix == '.gz' else path.suffix
    with opener(path, 'rt', newline='') as inp:
        lines = (line for line in inp if line.strip() and not line.lstrip().startswith('#'))
        rows = csv.reader(lines, delimiter=',' if suffix == '.csv' else '\t') if suffix in ('.csv', '.tsv') else (line.split() for line in lines)
        header = next(rows, None)
        if header:
            header = [name.strip() for name in header]
        if not header or len(header) < 2 or not all(header) or len(set(header)) != len(header):
            raise ValueError(f'{path}: need a header with distinct column names and at least two columns')
        values = []
        seen = set()
        for number, row in enumerate(rows, 2):
            if len(row) != len(header):
                raise ValueError(f'{path}:record {number}: expected {len(header)} columns, got {len(row)}')
            row = [value.strip() for value in row]
            if not row[0] or row[0] in seen:
                raise ValueError(f'{path}:record {number}: empty or duplicate identifier {row[0]!r}')
            seen.add(row[0])
            values.append(row)
    if not values:
        raise ValueError(f'{path}: no data records')
    return header, values


def numeric(value, context, allow_missing=False):
    if allow_missing and value in MISSING:
        return np.nan
    try:
        result = float(value)
    except ValueError:
        raise ValueError(f'{context}: expected a numeric value, got {value!r}; censored values need explicit preprocessing') from None
    if not math.isfinite(result):
        raise ValueError(f'{context}: non-finite value {value!r}')
    return result


def load_data(args):
    activity_header, activity_rows = read_table(args.activity)
    column = args.activity_column or activity_header[1]
    if column not in activity_header[1:]:
        raise ValueError(f'Activity column {column!r} is not a measurement column')
    index = activity_header.index(column)
    activity = {row[0]: numeric(row[index], f'{args.activity}:{row[0]}') for row in activity_rows}
    header, rows = read_table(args.descriptors)
    names = header[1:]
    missing = [row[0] for row in rows if row[0] not in activity]
    if missing and not args.allow_missing_activity:
        raise ValueError(f'{len(missing)} descriptor identifiers lack activity, including {missing[:3]}; use --allow_missing_activity to skip')
    rows = sorted((row for row in rows if row[0] in activity), key=lambda row: row[0])
    ids = [row[0] for row in rows]
    y = np.array([activity[identifier] for identifier in ids], dtype=float)
    if len(y) < args.min_samples:
        raise ValueError(f'Only {len(y)} matched molecules; need at least --min_samples {args.min_samples}')
    if args.transform != 'none':
        if np.any(y <= 0):
            raise ValueError('Log transforms require strictly positive measurements')
        y = np.log10(y) * (-1 if args.transform == 'neglog10' else 1)
    if np.ptp(y) == 0:
        raise ValueError('Measurements do not vary after transformation')
    x = np.array([[numeric(value, f'{args.descriptors}:{row[0]}:{name}', True)
                   for name, value in zip(names, row[1:])] for row in rows], dtype=float)
    meta = dict(activity_column=column, matched=len(ids), skipped_missing_activity=len(missing),
                unused_activity_identifiers=len(set(activity) - set(ids)),
                transform=args.transform, direction=args.direction, units=args.units,
                descriptors=str(args.descriptors), activity=str(args.activity))
    return ids, names, x, y, meta


def bucket_edges(x, requested, minimum):
    """Adjacent groups of whole values; merge sparse groups without splitting ties."""
    values, counts = np.unique(x, return_counts=True)
    if len(values) < 2:
        return np.array([], dtype=float)
    if len(values) <= requested:
        groups = [[float(value), float(value), int(count)] for value, count in zip(values, counts)]
    else:
        # Quantiles are placed between distinct values, never inside a tie.
        cumulative = np.cumsum(counts)
        cuts = sorted(set(min(int(np.searchsorted(cumulative, target)), len(values) - 2)
                          for target in np.linspace(0, len(x), requested + 1)[1:-1]))
        groups = []
        start = 0
        for end in cuts + [len(values) - 1]:
            if end >= start:
                groups.append([float(values[start]), float(values[end]), int(sum(counts[start:end + 1]))])
                start = end + 1
    while len(groups) > 1 and min(group[2] for group in groups) < minimum:
        i = min(range(len(groups)), key=lambda j: (groups[j][2], j))
        if i == 0:
            neighbor = 1
        elif i == len(groups) - 1:
            neighbor = i - 1
        else:
            neighbor = min((i - 1, i + 1), key=lambda j: (groups[j][2], j))
        a, b = sorted((i, neighbor))
        groups[a:b + 1] = [[groups[a][0], groups[b][1], groups[a][2] + groups[b][2]]]
    return np.array([a[1] * 0.5 + b[0] * 0.5 for a, b in zip(groups, groups[1:])])


def bin_index(x, edges):
    return np.searchsorted(edges, x, side='left')


def eta_squared(x, y_matrix, edges):
    labels = bin_index(x, edges)
    mean = y_matrix.mean(axis=0)
    denominator = np.sum((y_matrix - mean) ** 2, axis=0)
    between = np.zeros(y_matrix.shape[1])
    for label in np.unique(labels):
        group = y_matrix[labels == label]
        between += len(group) * (group.mean(axis=0) - mean) ** 2
    return np.divide(between, denominator, out=np.zeros_like(between), where=denominator > 0)


def spline_model(x):
    unique = len(np.unique(x))
    if unique < 2:
        return make_pipeline(StandardScaler(), Ridge(alpha=1.0))
    knots = np.unique(np.quantile(np.asarray(x).ravel(), np.linspace(0, 1, min(5, unique))))
    # Tied integer values can collapse quantile knots. Distinct endpoints remain
    # even when most molecules have the same descriptor value.
    return make_pipeline(SplineTransformer(degree=min(3, unique - 1), knots=knots[:, None], include_bias=False),
                         Ridge(alpha=1.0))


def cv_predictions(x, targets, splits, args):
    """Fit all shuffled targets together; every preprocessing fit uses training x only."""
    per_repeat = []
    linear_repeat = []
    bucket_repeat = {count: [] for count in args.buckets}
    for folds in splits:
        error = np.zeros(targets.shape[1])
        linear_error = np.zeros_like(error)
        baseline = np.zeros_like(error)
        bucket_error = {count: np.zeros_like(error) for count in args.buckets}
        for train, test in folds:
            train_x, test_x = x[train, None], x[test, None]
            train_y, test_y = targets[train], targets[test]
            average = train_y.mean(axis=0)
            baseline += np.sum((test_y - average) ** 2, axis=0)
            model = spline_model(train_x)
            model.fit(train_x, train_y)
            error += np.sum((test_y - model.predict(test_x)) ** 2, axis=0)
            linear = make_pipeline(StandardScaler(), Ridge(alpha=1.0))
            linear.fit(train_x, train_y)
            linear_error += np.sum((test_y - linear.predict(test_x)) ** 2, axis=0)
            for count in args.buckets:
                edges = bucket_edges(x[train], count, args.min_bucket_size)
                train_bin, test_bin = bin_index(x[train], edges), bin_index(x[test], edges)
                prediction = np.tile(average, (len(test), 1))
                for group in np.unique(train_bin):
                    prediction[test_bin == group] = train_y[train_bin == group].mean(axis=0)
                bucket_error[count] += np.sum((test_y - prediction) ** 2, axis=0)
        def gain(errors):
            return np.divide(baseline - errors, baseline, out=np.zeros_like(errors), where=baseline > 0)
        per_repeat.append(gain(error))
        linear_repeat.append(gain(linear_error))
        for count in args.buckets:
            bucket_repeat[count].append(gain(bucket_error[count]))
    return (np.asarray(per_repeat), np.asarray(linear_repeat),
            {count: np.asarray(values) for count, values in bucket_repeat.items()})


def correlations(x, targets):
    centered = x - x.mean()
    y_centered = targets - targets.mean(axis=0)
    denominator = np.sqrt(np.sum(centered ** 2) * np.sum(y_centered ** 2, axis=0))
    return np.divide(centered @ y_centered, denominator,
                     out=np.zeros(targets.shape[1]), where=denominator > 0)


def composite(pearson, spearman, mi, eta, cv_gain):
    """Inspection priority, not a probability: 20% trend, 20% dependence, 60% prediction."""
    trend = np.maximum(pearson ** 2, spearman ** 2)
    dependence = (1.0 - np.exp(-2.0 * np.maximum(mi, 0)) + eta) * 0.5
    return 100 * (0.2 * trend + 0.2 * dependence + 0.6 * np.clip(cv_gain, 0, 1))


def bucket_summary(x, y, edges, indices):
    result = []
    labels = bin_index(x, edges)
    for group in np.unique(labels):
        selected = labels == group
        ys, xs = y[selected], x[selected]
        half = float(stats.t.ppf(0.975, len(ys) - 1) * stats.sem(ys)) if len(ys) > 1 else None
        result.append(dict(low=float(xs.min()), high=float(xs.max()), n=len(ys),
                           mean=float(ys.mean()), median=float(np.median(ys)),
                           sd=float(np.std(ys, ddof=1)) if len(ys) > 1 else None,
                           mean_ci=[float(ys.mean() - half), float(ys.mean() + half)] if half is not None else None,
                           members=indices[selected].tolist()))
    return result


def feature_profile(x, args):
    finite = np.isfinite(x)
    values = x[finite]
    unique, counts = np.unique(values, return_counts=True)
    profile = dict(n=len(values), missing=int(np.sum(~finite)), distinct=len(unique),
                   mode_fraction=float(max(counts) / len(values)) if len(values) else 1.0,
                   integer_valued=bool(len(values) and np.all(values == np.round(values))))
    reason = None
    if len(values) < args.min_samples:
        reason = 'insufficient nonmissing observations'
    elif len(values) / len(x) < 1 - args.max_missing_fraction:
        reason = 'too many missing values'
    elif len(unique) < 2:
        reason = 'constant'
    elif profile['mode_fraction'] >= args.max_mode_fraction:
        reason = 'near constant (dominant value)'
    elif np.std(values) <= args.min_relative_sd * max(np.max(np.abs(values)), np.finfo(float).tiny):
        reason = 'near constant (relative spread)'
    return finite, profile, reason


def correlated_descriptors(names, matrix, args):
    pairs = []
    # Most iwdescr tables are complete. Rank those columns once rather than
    # sorting them again for every pair. Pairs involving missing values still
    # need ranks recomputed on their common observed subset.
    complete = np.all(np.isfinite(matrix), axis=0)
    complete_indices = np.flatnonzero(complete)
    positions = {column: i for i, column in enumerate(complete_indices)}
    complete_correlation = None
    if len(complete_indices) > 1:
        complete_correlation = np.corrcoef(stats.rankdata(matrix[:, complete], axis=0), rowvar=False)
    # Pairwise-complete Spearman correlations are specific to this dataset.
    for i, name in enumerate(names):
        for j in range(i):
            if complete[i] and complete[j]:
                rho = float(complete_correlation[positions[i], positions[j]])
                count = len(matrix)
            else:
                good = np.isfinite(matrix[:, i]) & np.isfinite(matrix[:, j])
                count = int(np.sum(good))
                if count < args.min_samples:
                    continue
                a, b = matrix[good, i], matrix[good, j]
                if np.ptp(a) == 0 or np.ptp(b) == 0:
                    continue
                rho = float(stats.spearmanr(a, b).statistic)
            if abs(rho) >= args.correlation_threshold:
                pairs.append(dict(a=names[j], b=name, spearman=rho, n=count))
    return sorted(pairs, key=lambda pair: (-abs(pair['spearman']), pair['a'], pair['b']))


def forest_analysis(matrix, y, names, args):
    if args.no_forest or not names:
        return None
    values = {name: [] for name in names}
    skills = []
    for repeat in range(args.repeats):
        folds = KFold(n_splits=min(args.folds, len(y) // 2), shuffle=True, random_state=args.seed + repeat)
        for fold, (train, test) in enumerate(folds.split(y)):
            # Median imputation is learned inside the training fold.
            imputer = SimpleImputer(strategy='median', keep_empty_features=True)
            train_x = imputer.fit_transform(matrix[train])
            test_x = imputer.transform(matrix[test])
            model = RandomForestRegressor(n_estimators=args.trees, min_samples_leaf=3,
                                           random_state=args.seed + repeat, n_jobs=args.jobs)
            model.fit(train_x, y[train])
            baseline = np.mean((y[test] - np.mean(y[train])) ** 2)
            error = np.mean((y[test] - model.predict(test_x)) ** 2)
            skills.append(1 - error / baseline if baseline > 0 else 0.0)
            importance = permutation_importance(model, test_x, y[test], scoring='neg_mean_squared_error',
                                                n_repeats=3, random_state=args.seed + repeat + fold,
                                                n_jobs=args.jobs)
            variance = baseline if baseline > 0 else 1.0
            for i, name in enumerate(names):
                values[name].append(float(importance.importances_mean[i] / variance))
    return dict(cv_gain=float(np.mean(skills)), cv_gains=skills,
                importance={name: dict(mean=float(np.mean(scores)), sd=float(np.std(scores)))
                            for name, scores in values.items()},
                note='Held-out permutation MSE increase divided by fold training-mean baseline MSE. Correlated substitutes can each appear unimportant.')


def analyse(args):
    ids, names, matrix, y, metadata = load_data(args)
    rng = np.random.default_rng(args.seed)
    # One common shuffle per replicate preserves the across-descriptor dependence
    # when calibrating the best score found anywhere in the table.
    targets = np.column_stack([y] + [rng.permutation(y) for _ in range(args.permutations)])
    rows, omitted, retained = [], [], []
    for column, name in enumerate(names):
        x = matrix[:, column]
        finite, profile, reason = feature_profile(x, args)
        if reason:
            omitted.append(dict(name=name, reason=reason, **profile))
            continue
        local_x, local_targets = x[finite], targets[finite]
        # If a feature's observed subset has no activity variation, its association
        # cannot be assessed. Other descriptors can still be analysed.
        if np.ptp(local_targets[:, 0]) == 0:
            omitted.append(dict(name=name, reason='constant activity in observed subset', **profile))
            continue
        print(f'Analysing {name} ({len(rows) + 1})', file=sys.stderr)
        folds = min(args.folds, len(local_x) // 2)
        splits = [list(KFold(n_splits=folds, shuffle=True, random_state=args.seed + repeat).split(local_x))
                  for repeat in range(args.repeats)]
        pearson = correlations(local_x, local_targets)
        spearman = correlations(stats.rankdata(local_x), stats.rankdata(local_targets, axis=0))
        discrete = profile['integer_valued'] and profile['distinct'] <= max(20, int(math.sqrt(len(local_x))))
        mi = []
        for target in local_targets.T:
            # Discrete singleton levels cannot support the continuous-target MI
            # estimator. Treat high-cardinality counts as continuous instead.
            _, occurrences = np.unique(local_x, return_counts=True)
            use_discrete = discrete and np.all(occurrences > 1)
            mi.append(float(mutual_info_regression(local_x[:, None], target,
                          discrete_features=[bool(use_discrete)], n_neighbors=min(3, len(local_x) - 1),
                          random_state=args.seed)[0]))
        mi = np.asarray(mi)
        bucket_results = []
        eta = np.zeros(local_targets.shape[1])
        indices = np.flatnonzero(finite)
        for requested in args.buckets:
            edges = bucket_edges(local_x, requested, args.min_bucket_size)
            association = eta_squared(local_x, local_targets, edges)
            eta = np.maximum(eta, association)
            bucket_results.append(dict(requested=requested, actual=len(edges) + 1,
                                       eta_squared=float(association[0]),
                                       groups=bucket_summary(local_x, local_targets[:, 0], edges, indices)))
        gains, linear, bucket_cv = cv_predictions(local_x, local_targets, splits, args)
        average_gain = gains.mean(axis=0)
        scores = composite(pearson, spearman, mi, eta, average_gain)
        full_model = spline_model(local_x)
        full_model.fit(local_x[:, None], local_targets[:, 0])
        grid = np.linspace(local_x.min(), local_x.max(), 100)
        curve = full_model.predict(grid[:, None])
        # Limit plotted points without affecting any calculations or exports.
        order = np.argsort(local_x, kind='stable')
        plotted = order[np.unique(np.linspace(0, len(order) - 1, min(args.max_points, len(order))).astype(int))]
        qx = np.quantile(local_x, [0.25, 0.75])
        qpred = full_model.predict(qx[:, None])
        for result in bucket_results:
            result['cv_gain'] = float(bucket_cv[result['requested']][:, 0].mean())
        row = dict(name=name, cv_folds=folds, **profile, pearson=float(pearson[0]), spearman=float(spearman[0]),
                   mutual_information=float(mi[0]), mi_discrete=bool(use_discrete),
                   score=float(scores[0]), null_scores=scores[1:].tolist(),
                   spline_cv_gain=float(average_gain[0]), spline_cv_repeats=gains[:, 0].tolist(),
                   linear_cv_gain=float(linear[:, 0].mean()),
                   effect_iqr=float(qpred[1] - qpred[0]),
                   curve_range=float(np.ptp(curve)), buckets=bucket_results,
                   curve=list(map(list, zip(grid.tolist(), curve.tolist()))),
                   points=[[float(local_x[i]), float(local_targets[i, 0]), int(indices[i])] for i in plotted],
                   neighbors=[])
        rows.append(row)
        retained.append(column)
    if not rows:
        raise ValueError('No descriptors survived filtering')
    null_max = np.max([row['null_scores'] for row in rows], axis=0) if args.permutations else np.array([])
    for row in rows:
        null = np.asarray(row['null_scores'])
        row['permutation_p'] = float((1 + sum(null >= row['score'])) / (len(null) + 1)) if len(null) else None
        row['familywise_p'] = float((1 + sum(null_max >= row['score'])) / (len(null_max) + 1)) if len(null_max) else None
        stable = min(row['spline_cv_repeats']) > 0
        meaningful = row['spline_cv_gain'] >= args.min_cv_gain and row['curve_range'] >= args.min_effect
        supported = row['familywise_p'] is not None and row['familywise_p'] <= 0.05
        row['evidence'] = ('supported' if supported and meaningful and stable and row['n'] >= 30
                           else 'tentative' if meaningful and stable else 'weak')
        row['notes'] = []
        if row['n'] < 30:
            row['notes'].append('Small sample: associations and uncertainty are unstable')
        if not meaningful:
            row['notes'].append('Little held-out improvement or practical effect')
        if not stable:
            row['notes'].append('Spline does not improve on the mean in every repeat')
        if not supported:
            row['notes'].append('Not established beyond the best shuffled descriptor')
    retained_names = [names[i] for i in retained]
    pairs = correlated_descriptors(retained_names, matrix[:, retained], args)
    by_name = {row['name']: row for row in rows}
    for pair in pairs:
        by_name[pair['a']]['neighbors'].append(dict(name=pair['b'], spearman=pair['spearman'], n=pair['n']))
        by_name[pair['b']]['neighbors'].append(dict(name=pair['a'], spearman=pair['spearman'], n=pair['n']))
    forest = forest_analysis(matrix[:, retained], y, retained_names, args)
    for row in rows:
        row['forest_importance'] = forest['importance'][row['name']] if forest else None
    rows.sort(key=lambda row: (-row['score'], row['name']))
    warnings = []
    if len(y) < 30:
        warnings.append('Small dataset: treat rankings as hypotheses, not reliable property effects.')
    if args.permutations < 19:
        warnings.append('Fewer than 19 permutations: a 0.05 permutation threshold cannot be resolved.')
    if rows[0]['evidence'] != 'supported':
        warnings.append('The highest-ranked descriptor is not supported by all evidence checks; a winner is not proof of an association.')
    if forest and forest['cv_gain'] <= 0:
        warnings.append('The random forest does not outperform the training-mean baseline on held-out molecules; its importances are not reliable evidence.')
    metadata.update(seed=args.seed, folds=args.folds, repeats=args.repeats, permutations=args.permutations,
                    min_samples=args.min_samples, requested_buckets=args.buckets,
                    forest_enabled=not args.no_forest, forest_trees=args.trees, jobs=args.jobs,
                    max_plot_points=args.max_points,
                    software=dict(numpy=np.__version__, scipy=scipy_version, sklearn=sklearn_version),
                    min_bucket_size=args.min_bucket_size, max_mode_fraction=args.max_mode_fraction,
                    min_relative_sd=args.min_relative_sd, max_missing_fraction=args.max_missing_fraction,
                    min_cv_gain=args.min_cv_gain, min_effect=args.min_effect,
                    correlation_threshold=args.correlation_threshold, descriptor_count=len(names),
                    analysed_descriptors=len(rows), warnings=warnings,
                    score_definition='100 * (0.2*max(Pearson^2,Spearman^2) + 0.2*((1-exp(-2*MI)+max_bucket_eta2)/2) + 0.6*clip(spline_CV_gain,0,1))',
                    cv_definition='1 - held-out squared error / squared error of training-fold mean; repeated folds, fixed spline complexity',
                    calibration='Common activity shuffles; familywise p compares with maximum composite score across retained descriptors in each shuffle')
    data = dict(format_version=1, metadata=metadata, observations=dict(identifiers=ids, activity=y.tolist()),
                summary=dict(n=len(y), mean=float(np.mean(y)), sd=float(np.std(y, ddof=1)),
                             min=float(np.min(y)), max=float(np.max(y))),
                features=rows, omitted=omitted, correlations=pairs, forest=forest,
                null_max_scores=null_max.tolist())
    if args.smiles:
        attach_structures(data, args.smiles)
    return data


def attach_structures(data, path):
    """Optional structure enrichment; chemistry remains in LillyMol."""
    try:
        from lillymol import MolFromSmiles
    except ImportError as exc:
        raise ValueError('--smiles requires LillyMol bindings; use run_python.sh with staged bindings') from exc
    lookup = {}
    with open(path) as inp:
        for number, line in enumerate(inp, 1):
            if not line.strip() or line.startswith('#'):
                continue
            words = line.split()
            if len(words) < 2 or words[1] in lookup:
                raise ValueError(f'{path}:{number}: need SMILES and a distinct identifier')
            mol = MolFromSmiles(words[0])
            if mol is None:
                raise ValueError(f'{path}:{number}: invalid SMILES')
            lookup[words[1]] = dict(smiles=words[0], natoms=mol.natoms())
    data['observations']['structures'] = [lookup.get(identifier) for identifier in data['observations']['identifiers']]


def positive_integer(value):
    result = int(value)
    if result < 1:
        raise argparse.ArgumentTypeError('must be positive')
    return result


def nonnegative_integer(value):
    result = int(value)
    if result < 0:
        raise argparse.ArgumentTypeError('must be nonnegative')
    return result


def add_analysis_options(parser):
    parser.add_argument('descriptors', type=Path, help='Headered iwdescr table; first column is identifier')
    parser.add_argument('-A', '--activity', type=Path, required=True, help='Headered measurements; first column is identifier')
    parser.add_argument('--activity_column', help='Measurement header name (default: second column)')
    parser.add_argument('--allow_missing_activity', action='store_true')
    parser.add_argument('--transform', choices=('none', 'log10', 'neglog10'), default='none')
    parser.add_argument('--direction', choices=('higher', 'lower'), default='higher', help='Desirable direction on the transformed activity scale')
    parser.add_argument('--units', default='', help='Units/label of transformed activity')
    parser.add_argument('--smiles', type=Path, help='Optional SMILES/id file for examples; needs LillyMol bindings')
    parser.add_argument('--min_samples', type=positive_integer, default=8)
    parser.add_argument('--max_missing_fraction', type=float, default=0.2)
    parser.add_argument('--max_mode_fraction', type=float, default=0.99)
    parser.add_argument('--min_relative_sd', type=float, default=1e-12)
    parser.add_argument('--buckets', default='2,4,8', help='Comma-separated requested bucket counts')
    parser.add_argument('--min_bucket_size', type=positive_integer, default=5)
    parser.add_argument('--folds', type=positive_integer, default=5)
    parser.add_argument('--repeats', type=positive_integer, default=3)
    parser.add_argument('--permutations', type=nonnegative_integer, default=19, help='Activity shuffles (default 19; 99+ for finer p resolution; 0 disables)')
    parser.add_argument('--seed', type=int, default=2026)
    parser.add_argument('--min_cv_gain', type=float, default=0.02)
    parser.add_argument('--min_effect', type=float, default=0, help='Minimum descriptive curve range in transformed activity units')
    parser.add_argument('--correlation_threshold', type=float, default=0.85)
    parser.add_argument('--no_forest', action='store_true')
    parser.add_argument('--trees', type=positive_integer, default=100)
    parser.add_argument('--jobs', type=positive_integer, default=1)
    parser.add_argument('--max_points', type=positive_integer, default=400)
    parser.add_argument('-o', '--output', required=True, help='JSON report, or - for stdout')
    parser.add_argument('--html', help='Also write offline HTML')


def check_outputs(outputs, inputs):
    paths = [Path(path).resolve() for path in outputs if path and path != '-']
    if len(paths) != len(set(paths)):
        raise ValueError('Output paths must be distinct')
    for path in paths:
        if path in [Path(inp).resolve() for inp in inputs if inp]:
            raise ValueError(f'Output would overwrite input: {path}')


def write_text(path, text):
    if path == '-':
        sys.stdout.write(text)
    else:
        Path(path).write_text(text)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    analysis = commands.add_parser('analyse', help='Analyse descriptors and measurements')
    add_analysis_options(analysis)
    render = commands.add_parser('html', help='Render an existing JSON report')
    render.add_argument('input', type=Path)
    render.add_argument('-o', '--output', required=True)
    for command in (analysis, render):
        command.add_argument('--top_n', type=nonnegative_integer, default=50, help='HTML feature limit; 0 retains all')
        command.add_argument('--page_size', type=positive_integer, default=20)
        command.add_argument('--nomedian', action='store_true')
        command.add_argument('--max_structures', type=nonnegative_integer, default=100, help='Maximum embedded molecule images; other examples retain SMILES')
        command.add_argument('--title', default='Descriptor Activity Exploration')
    args = parser.parse_args(argv)
    try:
        if args.command == 'analyse':
            args.buckets = sorted(set(int(token) for token in args.buckets.split(',')))
            if not args.buckets or min(args.buckets) < 2:
                raise ValueError('--buckets must contain integers >=2')
            if args.min_samples < 4 or args.folds < 2:
                raise ValueError('--min_samples must be >=4 and --folds >=2')
            if (not 0 <= args.max_missing_fraction < 1 or not 0 < args.max_mode_fraction <= 1 or
                    not 0 <= args.correlation_threshold <= 1):
                raise ValueError('Invalid missingness, mode-frequency, or correlation threshold')
            for key in ('min_relative_sd', 'min_cv_gain', 'min_effect'):
                value = getattr(args, key)
                if not math.isfinite(value) or value < 0:
                    raise ValueError(f'--{key} must be finite and nonnegative')
            if args.output == '-' and args.html == '-':
                raise ValueError('JSON and HTML cannot both use stdout')
            check_outputs([args.output, args.html], [args.descriptors, args.activity, args.smiles])
            data = analyse(args)
        else:
            check_outputs([args.output], [args.input])
            data = json.loads(args.input.read_text())
        from descriptor_exploration_report import build_html
        if args.command == 'analyse':
            text = json.dumps(data, indent=2, allow_nan=False) + '\n'
            document = build_html(data, args) if args.html else None
            write_text(args.output, text)
            if document is not None:
                write_text(args.html, document)
        else:
            write_text(args.output, build_html(data, args))
    except (OSError, ValueError, KeyError) as exc:
        print(f'descriptor_exploration: {exc}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
