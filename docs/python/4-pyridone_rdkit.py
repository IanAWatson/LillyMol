# Implementation of the 4-pyridone algorithm in RDKit.
# Not sure if this is the most efficient implementation

from absl import app
from absl import logging

from rdkit import Chem

def isRingAromatic(mol, bondRing):
  for id in bondRing:
    if not mol.GetBondWithIdx(id).GetIsAromatic():
      return False
  return True


def four_pyridone(mol):
  ri = mol.GetRingInfo()
  atomrings = ri.AtomRings()
  bondrings = ri.BondRings()
  print(bondrings)
  if len(bondrings) == 0:
    return

  for ring_number, ring in enumerate(bondrings):
    print(len(ring))
    if len(ring) != 6:
      continue
    if not isRingAromatic(mol, bondrings[ring_number]):
      continue

  return

def main(argv):
  suppl = Chem.SmilesMolSupplier(argv[1])
  molecules_read = 0
  for mol in suppl:
    if mol is None:
      logging.warning("Skipping None molecule")
      continue
    four_pyridone(mol)
#   molecules_read += 1
#   if molecules_read > 10000:
#     break


if __name__ == "__main__":
  app.run(main)
