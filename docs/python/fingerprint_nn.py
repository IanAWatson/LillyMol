# Demonstrate similarity searching with EC fingerprints and python

from typing import List

import numpy as np

from lillymol import *
from lillymol_fingerprint import *
from lillymol_io import *

from absl import app
from absl import flags
from absl import logging

FLAGS = flags.FLAGS

flags.DEFINE_integer("nbrs", 10, "Number Neighbours")
flags.DEFINE_integer("nbits", 512, "Number Neighbours")
flags.DEFINE_string("atype", "UST:APT", "Number Neighbours")
flags.DEFINE_boolean("verbose", False, 'Verbose output')
flags.DEFINE_enum("fp", "LINEAR", ["EC", "AP", "LINEAR"], "Kind of fingerprint to use")

def usage():
  print("Finds nearest neighbours within a file of smiles")
  print("-nbrs <nbrs>      number of neighbours to find")
  print("-nbits <bits>     number of bits in each fingerprint")
  print("-atype <atype>    tom typing to use, UST:ACHY for example")
  exit(1)

def nearneighbours(mols: List[Molecule],
                   bits:List,
                   ndx:int,
                   nbrs:int):
  """Find the `nbrs` nearest neighbours of bits[ndx]
    Args:
      mols: list of molecules - holds the names and smiles.
      bits: list of numpy arrays of fingerprints.
      ndx: the index in `mols` and `nds` of the needle.
      nbrs: number of neighbours of `ndx` to write.
  """
  # Compute distances wrt `ndx`. This is the expensive part.
  dists = [tanimoto(bits[ndx], b) for b in bits]

  # Get sort order - it is sorted lowest similarity to highest.
  # Should be able to use a partial sort.
  order = np.argsort(dists)

  # Print the needle.
  print(mols[ndx].smiles() + ' ' + mols[ndx].name())

  # Write neighbours.
  for i in range(0, nbrs):
    # Don't write self neighbours
    i = order[len(mols) - i - 1]
    if i == ndx:
      continue
    print(f"{mols[i].smiles()} {mols[i].name()} {dists[i]:.3f}")

def main(argv):
  """Demo application for using EC type fingerprints for a similarity search
  """
  if len(argv) == 1:
    logging.info("Must specify input smiles file")
    usage()

  mols = []
  with ReaderContext(argv[1]) as reader:
    for mol in reader:
      mols.append(mol)
  logging.info("Read %d moledules", len(mols))

  # Choose which kind of fingerprint to use.
  if FLAGS.fp == "LINEAR":
    fp_creator = LinearFingerprintCreator(FLAGS.nbits)
  elif FLAGS.fp == "EC":
    fp_creator = ECFingerprintCreator(FLAGS.nbits)
  elif FLAGS.fp == "AP":
    fp_creator = AtomPairFingerprintCreator(FLAGS.nbits)

  fp_creator.set_atom_type(FLAGS.atype)

  bits = [fp_creator.fingerprint(mol) for mol in mols]
  logging.info("Generated %d fingerprints", len(bits))

  for i in range(0, len(mols)):
    nearneighbours(mols, bits, i, FLAGS.nbrs)

if __name__ == "__main__":
  app.run(main)
