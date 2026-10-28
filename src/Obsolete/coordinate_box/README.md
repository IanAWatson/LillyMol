# Coordinate box experiment

Archived from Molecule_Lib because representing a position in space with a
single integer did not prove useful enough to retain in the library.

The implementation, schema, and tests are preserved here for reference.
`smiles_integration.cc` contains the removed SMILES helper functions and output
option branch; `smiles_coordinate_box_test.cc` preserves their integration tests.
These files are not part of the active Bazel or CMake builds and are not expected
to compile independently.

The active library no longer reads boxed `{{B...}}` coordinates or exposes the
`cbox` SMILES output option. Ordinary `{{x,y,z}}` coordinate SMILES and the
`coords` output option remain supported.
