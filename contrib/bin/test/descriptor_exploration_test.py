"""Run with Python containing numpy, scipy and scikit-learn."""
from contextlib import redirect_stderr
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

import numpy as np

BIN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BIN))
import descriptor_exploration as exploration
from descriptor_exploration_report import build_html, SCRIPT


def options(work, **overrides):
    result = dict(descriptors=work/'descriptors.w', activity=work/'activity.txt',
                  activity_column=None, allow_missing_activity=False, transform='none',
                  direction='higher', units='', min_samples=8, max_missing_fraction=0.2,
                  max_mode_fraction=0.99, min_relative_sd=1e-12, buckets=[2, 4, 8],
                  min_bucket_size=5, folds=3, repeats=2, permutations=3, seed=123,
                  min_cv_gain=0.02, min_effect=0, correlation_threshold=0.85,
                  no_forest=True, trees=15, jobs=1, max_points=40, smiles=None,
                  top_n=50, page_size=20, nomedian=False, title='Descriptor report', max_structures=0)
    result.update(overrides)
    return SimpleNamespace(**result)


def write_data(work, x, y, names=None):
    if names is None:
        names = [f'x{i}' for i in range(x.shape[1])]
    (work/'descriptors.w').write_text('Name '+' '.join(names)+'\n'+''.join(
        f'm{i} '+' '.join('.' if np.isnan(v) else str(v) for v in row)+'\n'
        for i, row in enumerate(x)))
    (work/'activity.txt').write_text('ID activity\n'+''.join(f'm{i} {value}\n' for i,value in enumerate(y)))


class BucketTest(unittest.TestCase):
    def test_integer_values_are_never_split_and_sparse_neighbors_merge(self):
        x = np.repeat([0, 1, 2, 3], [40, 2, 2, 20])
        edges = exploration.bucket_edges(x, 8, 5)
        labels = exploration.bin_index(x, edges)
        for value in np.unique(x):
            self.assertEqual(len(np.unique(labels[x == value])), 1)
        self.assertTrue(all(n >= 5 for n in np.bincount(labels)))
        self.assertLess(len(edges)+1, 4)
        self.assertTrue(np.all(np.diff(labels) >= 0))

    def test_quantile_buckets_and_constant_input(self):
        x = np.arange(80)
        edges = exploration.bucket_edges(x, 8, 5)
        self.assertEqual(len(edges), 7)
        self.assertTrue(np.all(np.bincount(exploration.bin_index(x, edges)) >= 5))
        self.assertEqual(len(exploration.bucket_edges(np.ones(10), 8, 5)), 0)

    def test_eta_detects_nonlinear_bucket_means(self):
        x = np.repeat([0, 1, 2, 3], 10)
        y = np.repeat([4., 0., 0., 4.], 10)
        eta = exploration.eta_squared(x, y[:,None], exploration.bucket_edges(x, 4, 5))
        self.assertAlmostEqual(eta[0], 1.)
        self.assertAlmostEqual(exploration.correlations(x, y[:,None])[0], 0.)


class AnalysisTest(unittest.TestCase):
    def analyse(self, x, y, names=None, **extra):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            write_data(work, x, y, names)
            with redirect_stderr(io.StringIO()):
                return exploration.analyse(options(work, **extra))

    def test_linear_signal_and_correlated_descriptors(self):
        rng = np.random.default_rng(4)
        x = np.linspace(-3, 3, 100)
        result = self.analyse(np.column_stack([x, 2*x, rng.normal(size=len(x))]),
                              2*x+rng.normal(0, 0.1, len(x)), ['linear', 'duplicate', 'noise'])
        by_name = {r['name']:r for r in result['features']}
        self.assertGreater(by_name['linear']['spline_cv_gain'], .95)
        self.assertGreater(by_name['linear']['score'], by_name['noise']['score'])
        self.assertAlmostEqual(result['correlations'][0]['spearman'], 1)
        self.assertEqual(by_name['linear']['neighbors'][0]['name'], 'duplicate')
        np.testing.assert_allclose(result['null_max_scores'],
                                   np.max([r['null_scores'] for r in result['features']], axis=0))
        self.assertIn('fewer than 19', ' '.join(result['metadata']['warnings']).lower())

    def test_u_shape_detected_when_signed_correlation_is_zero(self):
        x = np.linspace(-3,3,120)
        result = self.analyse(x[:,None], x*x)
        row = result['features'][0]
        self.assertLess(abs(row['pearson']), .01)
        self.assertLess(abs(row['spearman']), .01)
        self.assertGreater(row['spline_cv_gain'], .9)
        self.assertGreater(row['score'], 50)
        self.assertEqual(len(row['points']), 40)

    def test_constant_and_near_constant_are_reported_as_omitted(self):
        rng = np.random.default_rng(2)
        x = rng.normal(size=100)
        rare = np.zeros(100);rare[0]=1
        result = self.analyse(np.column_stack([x,np.ones(100),rare]), x,
                              ['signal', 'constant', 'rare'])
        omitted = {r['name']:r['reason'] for r in result['omitted']}
        self.assertEqual(omitted['constant'], 'constant')
        self.assertIn('near constant', omitted['rare'])
        self.assertEqual([r['name'] for r in result['features']], ['signal'])

    def test_missing_values_excluded_and_counts_export_correctly(self):
        x=np.arange(30,dtype=float)
        y=x.copy();x[:4]=np.nan
        result=self.analyse(x[:,None],y)
        row=result['features'][0]
        self.assertEqual(row['n'],26)
        self.assertEqual(row['missing'],4)
        members=[i for g in row['buckets'][0]['groups'] for i in g['members']]
        obs=result['observations']['identifiers']
        self.assertEqual({obs[i] for i in members},{f'm{i}' for i in range(4,30)})

    def test_small_noise_dataset_does_not_get_supported_winner(self):
        rng=np.random.default_rng(7)
        result=self.analyse(rng.normal(size=(12,3)),rng.normal(size=12),permutations=19)
        self.assertTrue(all(r['evidence']!='supported' for r in result['features']))
        self.assertTrue(any('Small dataset' in message for message in result['metadata']['warnings']))
        self.assertTrue(any('winner' in message for message in result['metadata']['warnings']))

    def test_full_size_noise_is_marked_weak(self):
        rng=np.random.default_rng(19)
        result=self.analyse(rng.normal(size=(150,3)),rng.normal(size=150),permutations=19)
        self.assertEqual(result['features'][0]['evidence'],'weak')
        self.assertTrue(any('winner' in message for message in result['metadata']['warnings']))

    def test_pairwise_correlations_use_only_common_observations(self):
        x=np.arange(30,dtype=float)
        matrix=np.column_stack([x,2*x,-x])
        matrix[:3,0]=np.nan
        matrix[3:6,1]=np.nan
        pairs=exploration.correlated_descriptors(['a','b','c'],matrix,options(Path('/tmp')))
        pair=next(p for p in pairs if p['a']=='a' and p['b']=='b')
        self.assertEqual(pair['n'],24)
        self.assertAlmostEqual(pair['spearman'],1.)
        pair=next(p for p in pairs if p['a']=='a' and p['b']=='c')
        self.assertEqual(pair['n'],27)
        self.assertAlmostEqual(pair['spearman'],-1.)

    def test_filtering_records_missingness_and_relative_spread(self):
        args=options(Path('/tmp'),min_relative_sd=1e-6)
        x=np.arange(100,dtype=float);x[:30]=np.nan
        _,profile,reason=exploration.feature_profile(x,args)
        self.assertEqual(profile['missing'],30)
        self.assertEqual(reason,'too many missing values')
        _,_,reason=exploration.feature_profile(1+np.linspace(0,1e-9,100),args)
        self.assertIn('relative spread',reason)

    def test_discrete_mi_and_no_permutations(self):
        x=np.repeat(np.arange(4),15)
        rng=np.random.default_rng(1)
        result=self.analyse(x[:,None],x+rng.normal(0,.2,len(x)),permutations=0)
        row=result['features'][0]
        self.assertTrue(row['mi_discrete'])
        self.assertIsNone(row['familywise_p'])
        self.assertNotEqual(row['evidence'],'supported')
        self.assertEqual(result['null_max_scores'],[])

    def test_permutation_calibration_and_determinism(self):
        x=np.linspace(-2,2,60);y=x*x
        a=self.analyse(x[:,None],y,permutations=19)
        b=self.analyse(x[:,None],y,permutations=19)
        self.assertEqual(a['features'],b['features'])
        row=a['features'][0]
        self.assertEqual(row['familywise_p'],.05)
        self.assertEqual(row['evidence'],'supported')
        # A null maximum is evaluated across all descriptor scores for each common shuffle.
        self.assertEqual(a['null_max_scores'],row['null_scores'])

    def test_forest_reports_held_out_permutation_importance(self):
        rng=np.random.default_rng(11)
        x=rng.normal(size=(60,2));y=3*x[:,0]+rng.normal(0,.1,60)
        x[0,1]=np.nan
        result=self.analyse(x,y,['signal','noise'],no_forest=False,repeats=1,permutations=0)
        forest=result['forest']
        self.assertGreater(forest['cv_gain'],.5)
        self.assertGreater(forest['importance']['signal']['mean'],forest['importance']['noise']['mean'])
        self.assertEqual(len(forest['cv_gains']),3)

    def test_three_integer_levels_can_have_a_middle_optimum(self):
        x=np.repeat([0.,1.,2.],30)
        result=self.analyse(x[:,None],(x-1)**2)
        row=result['features'][0]
        self.assertAlmostEqual(row['pearson'],0.)
        self.assertGreater(row['spline_cv_gain'],.8)

    def test_sparse_integer_levels_do_not_create_duplicate_spline_knots(self):
        x=np.array([0.]*80+[1.]*6+[2.]*6+[3.]*6+[4.]*2)
        result=self.analyse(x[:,None],x+np.random.default_rng(6).normal(0,.1,len(x)))
        self.assertEqual(len(result['features']),1)
        self.assertTrue(np.isfinite(result['features'][0]['spline_cv_gain']))

    def test_cv_transform_is_fitted_on_training_only(self):
        x=np.arange(12,dtype=float);targets=np.column_stack([x,x[::-1]])
        splits=[[(np.arange(8),np.arange(8,12)),(np.arange(4,12),np.arange(4))]]
        original=exploration.SplineTransformer.fit
        fitted=[]
        def capture(transformer,values,*args,**kwargs):
            fitted.append(values.copy())
            return original(transformer,values,*args,**kwargs)
        with mock.patch.object(exploration.SplineTransformer,'fit',capture):
            exploration.cv_predictions(x,targets,splits,options(Path('/tmp')))
        self.assertEqual([values.ravel().tolist() for values in fitted],
                         [list(range(8)),list(range(4,12))])


class InputAndReportTest(unittest.TestCase):
    def test_input_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            work=Path(tmp)
            write_data(work,np.arange(10)[:,None],np.arange(10))
            (work/'activity.txt').write_text('ID activity\nm0 1\nm0 2\n')
            with self.assertRaisesRegex(ValueError,'duplicate'):
                exploration.load_data(options(work))
            (work/'activity.txt').write_text('ID activity\n'+''.join(f'm{i} >10\n' for i in range(10)))
            with self.assertRaisesRegex(ValueError,'censored'):
                exploration.load_data(options(work))
            (work/'activity.txt').write_text('ID activity\n'+''.join(f'm{i} inf\n' for i in range(10)))
            with self.assertRaisesRegex(ValueError,'non-finite'):
                exploration.load_data(options(work))

    def test_transforms_header_selection_and_join_by_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            work=Path(tmp)
            write_data(work,np.arange(10)[:,None],np.arange(10))
            (work/'activity.txt').write_text('ID unused measured\n'+''.join(f'm{i} 0 {10.**i}\n' for i in range(9,-1,-1)))
            ids,_,_,y,_=exploration.load_data(options(work,activity_column='measured',transform='neglog10'))
            self.assertEqual(ids,[f'm{i}' for i in range(10)])
            np.testing.assert_allclose(y,-np.arange(10))

    def test_missing_activity_requires_explicit_skip(self):
        with tempfile.TemporaryDirectory() as tmp:
            work=Path(tmp)
            write_data(work,np.arange(10)[:,None],np.arange(10))
            (work/'activity.txt').write_text('ID activity\n'+''.join(f'm{i} {i}\n' for i in range(9)))
            with self.assertRaisesRegex(ValueError,'lack activity'):
                exploration.load_data(options(work))
            self.assertEqual(len(exploration.load_data(options(work,allow_missing_activity=True))[0]),9)

    def test_cli_json_html_and_no_input_clobber(self):
        with tempfile.TemporaryDirectory(prefix='descriptor exploration ') as tmp:
            work=Path(tmp)
            x=np.linspace(-3,3,40)
            write_data(work,x[:,None],x*x,names=['</script><script>bad()</script>'])
            output,report=work/'report.json',work/'report.html'
            command=['analyse',str(work/'descriptors.w'),'-A',str(work/'activity.txt'),
                     '-o',str(output),'--html',str(report),'--no_forest','--permutations','0','--repeats','1','--nomedian']
            with redirect_stderr(io.StringIO()):
                self.assertEqual(exploration.main(command),0)
            data=json.loads(output.read_text());html=report.read_text()
            self.assertEqual(data['summary']['n'],40)
            self.assertIn('Inspect / export',html)
            self.assertNotIn('</script><script>bad()',html)
            self.assertIn('\\u003c/script>',html)
            self.assertIn('"show_median": false',html)
            self.assertIn('mean_ci',html)
            import shutil
            if shutil.which('node'):
                subprocess.run(['node', str(BIN/'test/descriptor_exploration_ui_test.js'), str(report)],
                               check=True, capture_output=True, text=True)
            original=(work/'descriptors.w').read_bytes()
            with redirect_stderr(io.StringIO()):
                self.assertEqual(exploration.main(command+['-o',str(work/'descriptors.w')]),1)
                self.assertEqual(exploration.main(command+['--max_missing_fraction','nan']),1)
            self.assertEqual((work/'descriptors.w').read_bytes(),original)
            other=work/'other.html'
            self.assertEqual(exploration.main(['html',str(output),'-o',str(other),'--top_n','1']),0)
            self.assertTrue(other.exists())

    def test_javascript_syntax(self):
        import shutil
        if not shutil.which('node'):
            self.skipTest('node is not installed')
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'report.js';path.write_text(SCRIPT)
            subprocess.run(['node','--check',str(path)],check=True,capture_output=True)


if __name__=='__main__':
    unittest.main()
