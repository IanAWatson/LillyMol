# Molecule Filter
This tool is designed for the common task of filtering molecules according
to a set of limits on computed properties.

While the tool is designed for efficiency, it is configurable to the point
where relatively expensive properties need to be requested. But this is
fully under the control of the user.

The properties that are available as filters are described in the proto
[molecule_filter.proto(molecule_filter.proto). The list includes obvious
candidates, min and max heavy atom count, rotatable bonds, Lipinski parameters,
xlogp and others. The relative expense of these is noted in the proto. 

## HowTo
The tool works via a textproto input that specifies the property ranges
that define the filter. A typical invocation might be
```
molecule_filter -F molecule_filter.textproto -g all -v -l -B /tmp/bad /tmp/rand100k.smi 
```
Chemical standardisation is turned on, and this is relatively expensive. If
the molecules have already been standardised, definitely omit this step.

The `-l` option reduces to the largest fragment, but again, if the molecules
have already been reduced to the largest frgament, omit this requirement.

The config textproto file is specified via the `-F` option. The version used
for the timing tests below looks like
```
# Example input configuration file for molecule_filter

min_natoms: 10
max_natoms: 40

min_nrings: 1
max_nrings: 5

min_heteroatom_count: 2

min_heteroatom_fraction: 0.0
max_heteroatom_fraction: 0.5

min_aromatic_ring_count: 1
max_aromatic_ring_count: 4

max_aliphatic_ring_count: 3

max_ring_system_size: 4
largest_ring_size: 7

exclude_non_organic: true
exclude_isotopes: true

min_rotatable_bonds: 1
max_rotatable_bonds: 9

min_tpsa: 30
max_tpsa: 160

min_xlogp: 0.0
max_xlogp: 6.0

min_hba: 1
max_hba: 10

# min_hbd: 0
max_hbd: 6

max_halogen_count: 6

max_distance: 20
min_sp3_carbon: 1
max_aromatic_density: 0.85
```
In all cases, excluding conditions that are not needed is always desirable.

## Details
The tool is written for efficiency, and at each step, care is taken to
try and do the minimum amount of work possible.

LillyMol tool `iwdescr` can also perform filtering based on the value
of computed descriptors, but in that case, it performs calculations
before doing any filtering, and individual calculations can generally
not be turned on and off.

When I originally conceived this, helping people process Enamine Real,
I thought that including substructure search filters would be desirable,
but you will be better off piping the output of this tool to `tsubstructure`
thereby gaining use of multiple cores - smiles interpretation of the
100k molecules takes about 0.6 seconds.

## Timing
All timing data reported here is done on 
```
model name      : Intel(R) Core(TM) i7-8700K CPU @ 3.70GHz

```
which is a consumer level system from 2017.

Start with 100k random molecules from chembl 33,
```
100000 molecules had between 1 and 633 atoms. Average 32.9515
```

With every property activated, and set to a value where it removes
only a small number of molecules, the time is about 18.7 seconds, or
about 1 Million molecules in 3 minutes.

Omitting the xlogp filter drops the time to just 10.8 seconds, and so
we see that this is the most expensive of the constraints that can 
be imposed.

Dropping any of the other constraints really does not change the
timing that much.

On a more realistic set, the 2.1M Enamine HTS molecules
```
2144530 molecules had between 3 and 63 atoms. Average 24.3974
```
This is processed in 112 seconds, including an xlogp cutoff, 
about 1.15M per second. This is substantially faster than
what we saw from Chembl. This is probably largely attributable
to the differing atom counts between the two collections, 32.9 vs 24.4.

