import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(__file__))

from lillymol_tools import GFPList


def _write_file(contents):
    fd, fname = tempfile.mkstemp(suffix='.gfp')
    with os.fdopen(fd, 'w') as output:
        output.write(contents)
    return fname


def _testdata_file(fname):
    candidates = [
        os.path.join(os.path.dirname(__file__), 'testdata', fname),
    ]
    test_srcdir = os.environ.get('TEST_SRCDIR')
    if test_srcdir:
        candidates.extend([
            os.path.join(test_srcdir, '_main', 'pybind', 'testdata', fname),
            os.path.join(test_srcdir, 'pybind', 'testdata', fname),
        ])
    runfiles_dir = os.environ.get('RUNFILES_DIR')
    if runfiles_dir:
        candidates.extend([
            os.path.join(runfiles_dir, '_main', 'pybind', 'testdata', fname),
            os.path.join(runfiles_dir, 'pybind', 'testdata', fname),
        ])
    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate
    raise FileNotFoundError(fname)


CARBON_GFP = """
$SMI<C>
PCN<Methane>
FCTS<.E..........2;1;1;1;1>
|
$SMI<CC>
PCN<Ethane>
FCTS<.U..........2;1;1;1;1>
|
"""


class TestGFPList(unittest.TestCase):

    def test_read_and_distance(self):
        fname = _write_file(CARBON_GFP)
        try:
            gfp = GFPList.from_file(fname)
            self.assertEqual(len(gfp), 2)
            self.assertEqual(gfp.size(), 2)
            self.assertEqual(gfp.tags(), ['FCTS<'])
            self.assertEqual(gfp.smiles(0), 'C')
            self.assertEqual(gfp.id(1), 'Ethane')
            self.assertAlmostEqual(gfp.distance(0, 0), 0.0, places=6)
            self.assertAlmostEqual(gfp.distance(0, 1), 0.5, places=6)
            self.assertAlmostEqual(gfp.distance(1, 0), 0.5, places=6)
        finally:
            os.remove(fname)

    def test_nearest_neighbours(self):
        fname = _write_file(CARBON_GFP)
        try:
            gfp = GFPList.from_file(fname)
            hits = gfp.nearest_neighbours(0, 1)
            self.assertEqual(len(hits), 1)
            self.assertEqual(hits[0].index, 1)
            self.assertAlmostEqual(hits[0].distance, 0.5, places=6)

            close = gfp.nearest_neighbours_within_distance(0, 0.5)
            self.assertEqual(len(close), 1)
            self.assertEqual(close[0].index, 1)
            self.assertAlmostEqual(close[0].distance, 0.5, places=6)

            none = gfp.nearest_neighbours_within_distance(0, 0.49)
            self.assertEqual(none, [])
        finally:
            os.remove(fname)

    def test_complex_gfp_file(self):
        gfp = GFPList.from_file(_testdata_file('rand10.complex.gfp'))
        self.assertEqual(len(gfp), 10)
        tags = gfp.tags()
        self.assertGreaterEqual(len(tags), 9)
        self.assertIn('FPIW<', tags)
        self.assertIn('FPMK<', tags)
        self.assertIn('MPR<', tags)

        self.assertAlmostEqual(gfp.distance(0, 0), 0.0, places=6)
        d01 = gfp.distance(0, 1)
        self.assertGreaterEqual(d01, 0.0)
        self.assertLessEqual(d01, 1.0)
        self.assertAlmostEqual(d01, gfp.distance(1, 0), places=6)

        hits = gfp.nearest_neighbours(0, 3)
        self.assertEqual(len(hits), 3)
        self.assertNotIn(0, [hit.index for hit in hits])
        self.assertEqual([hit.distance for hit in hits], sorted(hit.distance for hit in hits))

        within = gfp.nearest_neighbours_within_distance(0, 1.0)
        self.assertEqual(len(within), 9)
        self.assertNotIn(0, [hit.index for hit in within])

    def test_use_only_and_use_all(self):
        fname = _write_file(CARBON_GFP)
        try:
            gfp = GFPList.from_file(fname)
            gfp.use_only(['FCTS<'])
            self.assertAlmostEqual(gfp.distance(0, 1), 0.5, places=6)
            gfp.use_all()
            self.assertAlmostEqual(gfp.distance(0, 1), 0.5, places=6)
        finally:
            os.remove(fname)


if __name__ == '__main__':
    unittest.main()
