"""Local, isolated-process comparison of recursive SMARTS parsing and searching."""
import argparse
import json
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time

PANEL = [
 ('CC(=O)Oc1ccccc1C(=O)O', True),
 ('O=C(O)c1ccccc1O', False),
 ('CC(=O)Oc1ccccc1C(=O)OC', False),
 ('CC(=O)Oc1ccccc1', False),
 ('CC(=O)Oc1cccc(C(=O)O)c1', False),
 ('CC(=O)Oc1ccc(C(=O)O)cc1', False),
 ('CCC(=O)Oc1ccccc1C(=O)O', False),
 ('CC(=O)Oc1ccccc1C(=O)[O-]', False),
 ('CC(=O)Oc1ccccc1C(=O)O.[Na+]', True),
 ('[13CH3]C(=O)Oc1ccccc1C(=O)O', False),
 ('CC(=O)Oc1ccc(C)cc1C(=O)O', False),
]


def worker(engine, query_path):
    import resource
    text = Path(query_path).read_text()
    if engine == 'rdkit':
        from rdkit import Chem, __version__
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.error')
        molecules = [Chem.MolFromSmiles(s) for s, _ in PANEL]
        parameters = Chem.SubstructMatchParameters()
        parameters.maxMatches = 1
        parameters.numThreads = 1
        start = time.perf_counter()
        query = Chem.MolFromSmarts(text)
        parse = time.perf_counter() - start
        if query is None:
            raise ValueError('RDKit rejected generated SMARTS')
        match = lambda mol: mol.HasSubstructMatch(query, parameters)
        version = __version__
    else:
        from lillymol import MolFromSmiles, SubstructureQuery
        molecules = [MolFromSmiles(s) for s, _ in PANEL]
        for mol in molecules:
            mol.compute_aromaticity_if_needed()
            mol.nrings()
        query = SubstructureQuery()
        start = time.perf_counter()
        ok = query.build_from_smarts(text)
        parse = time.perf_counter() - start
        if not ok:
            raise ValueError('LillyMol rejected generated SMARTS')
        query.set_max_matches_to_find(1)
        match = lambda mol: query.substructure_search(mol) > 0
        version = 'staged nanobind'
    start = time.perf_counter()
    assert match(molecules[0])
    first = time.perf_counter() - start
    warm = []
    for _ in range(5):
        start = time.perf_counter()
        assert match(molecules[0])
        warm.append(time.perf_counter() - start)
    mixed = []
    for _ in range(3):
        start = time.perf_counter()
        result = [match(mol) for mol in molecules]
        elapsed = time.perf_counter() - start
        assert result == [expected for _, expected in PANEL], result
        mixed.append(elapsed)
    print(json.dumps(dict(engine=engine, version=version, parse_s=parse, first_s=first,
                         warm_s=statistics.median(warm), panel_s=statistics.median(mixed),
                         rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss /
                         (1024 ** 2 if sys.platform == 'darwin' else 1024))), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output_dir', type=Path, required=True,
                        help='Dedicated directory for generated queries and results; existing files are replaced')
    args = parser.parse_args()
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    repo = Path(__file__).resolve().parents[3]
    sys.path.insert(0, str(repo / 'contrib/bin'))
    from mkaspirin import AspirinModel, Generator
    model = AspirinModel()
    files = []
    for extensions in (False, True):
        generator = Generator(model=model, seed=19, lillymol_extensions=extensions)
        for goal in (1000, 10000, 100000, 1000000, 3000000):
            while generator.size < goal:
                generator.expand()
            path = root / (('extensions' if extensions else 'portable') + f'_{goal}.smarts')
            text = generator.smarts()
            assert len(text) == generator.size
            path.write_text(text)
            files.append(dict(path=str(path), mode='extensions' if extensions else 'portable',
                              bytes=len(text), depth=generator.depth, expansions=generator.expansions))
    output = dict(python=sys.version, platform=platform.platform(), records=[])
    for item in files:
        for engine in (('lillymol', 'rdkit') if item['mode'] == 'portable' else ('lillymol',)):
            runs = []
            for repeat in range(3):
                cmd=[sys.executable, str(Path(__file__).resolve()), '--worker', engine, item['path']]
                try:
                    result=subprocess.run(cmd, capture_output=True, text=True, timeout=60)
                except subprocess.TimeoutExpired:
                    runs.append(dict(error='60-second timeout'))
                    break
                if result.returncode:
                    runs.append(dict(error=f'exit {result.returncode}: {result.stderr[-1000:]}'))
                    break
                runs.append(json.loads(result.stdout))
            record=dict(item,engine=engine,runs=runs)
            if all('error' not in row for row in runs):
                for key in ('parse_s','first_s','warm_s','panel_s','rss_mib'):
                    record[key]=statistics.median(row[key] for row in runs)
                print(item['mode'],engine,item['bytes'],
                      'parse',round(record['parse_s'],6),
                      'first',round(record['first_s'],6),
                      'warm',round(record['warm_s'],6),
                      'panel',round(record['panel_s'],6),
                      'MiB',round(record['rss_mib'],1),flush=True)
            else:
                print(item['mode'],engine,item['bytes'],runs[-1]['error'],flush=True)
            output['records'].append(record)
            (root/'results.json').write_text(json.dumps(output,indent=2)+'\n')


if __name__ == '__main__':
    if len(sys.argv)>1 and sys.argv[1]=='--worker':
        worker(sys.argv[2],sys.argv[3])
    else:
        main()
