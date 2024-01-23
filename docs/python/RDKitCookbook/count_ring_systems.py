from lillymol import *

def get_ring_systems(mol:Molecule, include_spiro=False):
  # Return a list of the atoms in each ring system.
  if include_spiro:
    sysid = m.label_atoms_by_ring_system_including_spiro_fused()
  else:
    sysid = m.label_atoms_by_ring_system()

  result = []
  nsys = max(sysid)

  # non ring atoms get ring system id 0 so start at 1.
  for sys_num in range(1, nsys + 1):
    result.append([i for i,s in enumerate(sysid) if s == sys_num])

  return result

m = MolFromSmiles("CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3")
rsys = get_ring_systems(m)
for r in rsys:
  print(set(r))

