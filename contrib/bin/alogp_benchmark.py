#!/usr/bin/env python
"""Benchmark and compare LillyMol alogp with RDKit Crippen MolLogP.

Starting from a file of smiles (smiles, whitespace, identifier), this times
  - RDKit MolFromSmiles and Crippen.MolLogP
  - LillyMol alogp python bindings (parse, and ALogP.logp)
  - the LillyMol `alogp` command line tool (whole run)
each with LillyMol default settings and with the RDKit compatibility settings
(-Y alcacid -Y RDKIT.N+ -Y RDKIT.HP, or their python equivalents), and compares the
values with RDKit: correlation, mean difference, counts of molecules differing by
more than a threshold, the most common difference values, and the largest outliers.

  alogp_benchmark.py -o outdir molecules.smi

Use -s and -n to work on part of a file; the same range is given to the tool with its
-i skip= and -i do= options, so no temporary file is written.

Needs rdkit and numpy, and the LillyMol python bindings (run it with run_python.sh, or
with PYTHONPATH and LD_LIBRARY_PATH set to a directory staged by
copy_shared_libraries.sh). To test a freshly built tool or bindings without touching
the installed ones, stage the bindings in a scratch directory and set ALOGP to the
new `alogp` executable (or use -a).

Timings are single threaded wall clock for the whole file, and are from one run, so
treat ratios as approximate. Both the CLI and RDKit timings include parsing.
"""
import argparse
import collections
import json
import os
import subprocess
import sys
import time

try:
    import numpy as np
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Crippen
    from lillymol import ALogP, Molecule
except ImportError as e:
    sys.exit(f"alogp_benchmark.py: cannot import {e.name}; see the documentation at the top of this file")

RDKIT_FLAGS = ["-Y", "RDKIT.N+", "-Y", "RDKIT.HP", "-Y", "alcacid"]


def read_smiles(fname, skip, max_molecules):
    smiles, ids = [], []
    seen = 0
    with open(fname) as f:
        for line in f:
            tokens = line.split()
            if not tokens:
                continue
            seen += 1
            if seen <= skip:
                continue
            smiles.append(tokens[0])
            ids.append(tokens[1] if len(tokens) > 1 else str(len(ids)))
            if max_molecules and len(smiles) >= max_molecules:
                break
    return smiles, ids


def run_cli(alogp, flags, fname, ids):
    """Run the alogp tool on fname, return (seconds, values aligned with ids)."""
    t = time.perf_counter()
    p = subprocess.run([alogp, *flags, fname], capture_output=True, text=True)
    elapsed = time.perf_counter() - t
    if p.returncode != 0:
        sys.exit(f"alogp failed: {p.stderr[:500]}")
    values = {}
    for line in p.stdout.splitlines()[1:]:
        f = line.split()
        if len(f) >= 2:
            values[f[0]] = f[1]
    out = np.full(len(ids), np.nan)
    for i, name in enumerate(ids):
        v = values.get(name)
        if v is not None and v != ".":
            out[i] = float(v)
    return elapsed, out


def run_python_lillymol(smiles, rdkit_compatible):
    calc = ALogP()
    if rdkit_compatible:
        calc.set_rdkit_phoshoric_acid_hydrogen(1)
        calc.set_use_alcohol_for_acid(1)
    t = time.perf_counter()
    mols = []
    for s in smiles:
        m = Molecule()
        mols.append(m if m.build_from_smiles(s) else None)
    parse = time.perf_counter() - t
    t = time.perf_counter()
    values = [calc.logp(m) if m is not None else None for m in mols]
    logp = time.perf_counter() - t
    return parse, logp, np.array([np.nan if v is None else v for v in values], dtype=float)


def compare(name, x, ref, threshold):
    ok = ~np.isnan(x) & ~np.isnan(ref)
    a, b = x[ok], ref[ok]
    df = a - b
    return {
        "name": name, "n": int(ok.sum()),
        "pearson": float(np.corrcoef(a, b)[0, 1]),
        "mean_abs_diff": float(np.abs(df).mean()),
        "bias": float(df.mean()),
        "rmse": float(np.sqrt((df ** 2).mean())),
        "over_threshold": int((np.abs(df) > threshold).sum()),
        "over_1": int((np.abs(df) > 1.0).sum()),
        "max_abs_diff": float(np.abs(df).max()),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("smiles_file", help="smiles, whitespace, identifier")
    ap.add_argument("-o", "--outdir", required=True, help="directory for values.tsv and results.json")
    ap.add_argument("-a", "--alogp", default=os.environ.get("ALOGP"),
                    help="alogp executable (default $ALOGP, else $LILLYMOL_HOME/bin/Linux/alogp)")
    ap.add_argument("-s", "--skip", type=int, default=0, help="skip the first n molecules (passed to alogp as -i skip=n)")
    ap.add_argument("-n", "--max_molecules", type=int, default=0,
                    help="only use n molecules, after any skip (passed to alogp as -i do=n)")
    ap.add_argument("-t", "--threshold", type=float, default=0.1, help="difference counted as significant (default 0.1)")
    ap.add_argument("-k", "--top", type=int, default=15, help="number of largest outliers to list (default 15)")
    args = ap.parse_args()

    alogp = args.alogp
    if alogp is None:
        home = os.environ.get("LILLYMOL_HOME")
        if home is None:
            sys.exit("set -a, $ALOGP or $LILLYMOL_HOME to locate the alogp executable")
        alogp = os.path.join(home, "bin", "Linux", "alogp")
    os.makedirs(args.outdir, exist_ok=True)
    RDLogger.DisableLog("rdApp.*")

    smiles, ids = read_smiles(args.smiles_file, args.skip, args.max_molecules)
    fname = args.smiles_file
    # The tool reads the same molecules via its input options, no temporary file needed.
    subset = []
    if args.skip:
        subset += ["-i", f"skip={args.skip}"]
    if args.max_molecules:
        subset += ["-i", f"do={args.max_molecules}"]
    n = len(smiles)
    res = {"n": n, "alogp": alogp, "smiles_file": args.smiles_file}

    t_cli_def, cli_def = run_cli(alogp, subset, fname, ids)
    t_cli_rd, cli_rd = run_cli(alogp, subset + RDKIT_FLAGS, fname, ids)
    res["cli_default_s"], res["cli_rdkit_s"] = t_cli_def, t_cli_rd

    t = time.perf_counter()
    rmols = [Chem.MolFromSmiles(s) for s in smiles]
    res["rdkit_parse_s"] = time.perf_counter() - t
    t = time.perf_counter()
    rd = np.array([Crippen.MolLogP(m) if m is not None else np.nan for m in rmols])
    res["rdkit_logp_s"] = time.perf_counter() - t
    res["rdkit_failed"] = int(np.isnan(rd).sum())

    res["lly_parse_s"], res["lly_logp_default_s"], py_def = run_python_lillymol(smiles, False)
    _, res["lly_logp_rdkit_s"], py_rd = run_python_lillymol(smiles, True)

    rdkit_total = res["rdkit_parse_s"] + res["rdkit_logp_s"]
    res["speedup_logp_only"] = res["rdkit_logp_s"] / res["lly_logp_default_s"]
    res["speedup_python_total"] = rdkit_total / (res["lly_parse_s"] + res["lly_logp_default_s"])
    res["speedup_cli_total"] = rdkit_total / res["cli_default_s"]

    cols = {"rdkit": rd, "py_default": py_def, "py_rdkit_compat": py_rd, "cli_default": cli_def, "cli_rdkit_compat": cli_rd}
    with open(os.path.join(args.outdir, "values.tsv"), "w") as f:
        f.write("id\tsmiles\t" + "\t".join(cols) + "\n")
        for i in range(n):
            f.write(f"{ids[i]}\t{smiles[i]}\t" + "\t".join("" if np.isnan(c[i]) else f"{c[i]:.4f}" for c in cols.values()) + "\n")

    res["comparisons"] = [compare(k, v, rd, args.threshold) for k, v in cols.items() if k != "rdkit"]
    json.dump(res, open(os.path.join(args.outdir, "results.json"), "w"), indent=1)

    print(f"{n} molecules from {args.smiles_file}; alogp {alogp}")
    print(f"\nTimings (seconds, single thread)")
    print(f"  RDKit        parse {res['rdkit_parse_s']:8.2f}  logp {res['rdkit_logp_s']:8.2f}  ({res['rdkit_failed']} failed to parse)")
    print(f"  LillyMol py  parse {res['lly_parse_s']:8.2f}  logp {res['lly_logp_default_s']:8.2f}  (rdkit compat {res['lly_logp_rdkit_s']:.2f})")
    print(f"  LillyMol CLI whole run {res['cli_default_s']:.2f}  (rdkit compat {res['cli_rdkit_s']:.2f})")
    print(f"  speed up: logp only {res['speedup_logp_only']:.1f}x, parse+logp python {res['speedup_python_total']:.1f}x, CLI {res['speedup_cli_total']:.1f}x")

    print(f"\nAgreement with RDKit (LillyMol - RDKit); threshold {args.threshold}")
    print(f"  {'':18s} {'n':>7s} {'pearson':>8s} {'MAE':>7s} {'bias':>8s} {'RMSE':>7s} {'>thr':>6s} {'>1':>6s} {'max':>7s}")
    for c in res["comparisons"]:
        print(f"  {c['name']:18s} {c['n']:7d} {c['pearson']:8.5f} {c['mean_abs_diff']:7.4f} {c['bias']:+8.4f} {c['rmse']:7.4f} "
              f"{c['over_threshold']:6d} {c['over_1']:6d} {c['max_abs_diff']:7.3f}")

    df = py_rd - rd
    sig = np.abs(df) > args.threshold
    print(f"\nMost common differences among the {int(np.nansum(sig))} molecules over {args.threshold} (python, rdkit compat settings)")
    counts = collections.Counter(np.round(df[sig], 2))
    for value, count in counts.most_common(10):
        example = next(i for i in range(n) if sig[i] and round(df[i], 2) == value)
        print(f"  {value:+6.2f} {count:5d}  {smiles[example]}")

    print(f"\nLargest {args.top} differences")
    for i in np.argsort(-np.nan_to_num(np.abs(df)))[:args.top]:
        print(f"  {df[i]:+7.3f} lillymol {py_rd[i]:7.3f} rdkit {rd[i]:7.3f}  {ids[i]} {smiles[i]}")
    print(f"\nWrote {args.outdir}/values.tsv and results.json")


if __name__ == "__main__":
    main()
