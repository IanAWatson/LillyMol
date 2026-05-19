option(LILLYMOL_USE_EXISTING_PROTOBUF_FILES
  "Use existing .pb.cc/.pb.h files instead of running protoc"
  ON
)

lillymol_proto_library(atom_type_ext_cc_proto atom_type_ext.proto)
lillymol_proto_library(donor_acceptor_cc_proto donor_acceptor.proto)
lillymol_proto_library(etrans_cc_proto etrans.proto)
lillymol_proto_library(mol2graph_cc_proto mol2graph.proto)
lillymol_proto_library(geometric_constraints_cc_proto geometric_constraints.proto)

target_include_directories(geometric_constraints_cc_proto
  PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}
)

lillymol_proto_library(substructure_cc_proto
  substructure.proto
  geometric_constraints_cc_proto
)

lillymol_proto_library(molecule_to_query_cc_proto
  molecule_to_query.proto
  substructure_cc_proto
)

lillymol_proto_library(pharmacophore_cc_proto pharmacophore.proto)
lillymol_proto_library(standardise_cc_proto standardise.proto)

lillymol_proto_library(toggle_kekule_form_cc_proto
  toggle_kekule_form.proto
  substructure_cc_proto
)

lillymol_proto_library(reaction_cc_proto
  reaction.proto
  substructure_cc_proto
  toggle_kekule_form_cc_proto
)

add_library(space_vector INTERFACE)

target_sources(space_vector
  INTERFACE
    space_vector.h
)

target_include_directories(space_vector
  INTERFACE
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(space_vector
  INTERFACE
    iwmisc
)

add_library(coordinate_box
  coordinate_box.cc
)

target_sources(coordinate_box
  PUBLIC
    coordinate_box.h
)

target_include_directories(coordinate_box
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(coordinate_box
  PUBLIC
    space_vector
    iwstring
)

add_library(iwmolecule
  ISIS_Atom_List.cc
  _qry_wstats.cc
  _substructure.cc
  allowed_elements.cc
  aromatic.cc
  atom.cc
  atom_alias.cc
  atom_typing.cc
  bond.cc
  bond_list.cc
  careful_frag.cc
  charge_assigner.cc
  charge_calculation.cc
  chiral_centre.cc
  cis_trans_bond.cc
  coordinate_box.cc
  coordinates.cc
  count_atoms_in_smiles.cc
  csubstructure.cc
  dihedral.cc
  donor_acceptor.cc
  down_the_bond.cc
  element.cc
  element_hits_needed.cc
  ematch.cc
  etrans.cc
  frag.cc
  geometric_constraints.cc
  is_actually_chiral.cc
  iwrcb.cc
  iwrnm.cc
  iwsubstructure.cc
  mdl_atom_record.cc
  mdl_file_data.cc
  mdl_molecule.cc
  misc2.cc
  mol2graph.cc
  mol2graph_proto.cc
  molecule.cc
  molecule_smarts.cc
  molecule_to_query.cc
  moleculeb.cc
  moleculed.cc
  moleculeh.cc
  moleculeio.cc
  moleculer.cc
  numass.cc
  ostream_and_type.cc
  output.cc
  parse_smarts_tmp.cc
  path.cc
  path_around_ring.cc
  path_scoring.cc
  pearlman.cc
  qry_wstats.cc
  ring_bond_iterator.cc
  rmele.cc
  rotbond_common.cc
  rwsubstructure.cc
  rwsubstructure_proto.cc
  set_of_atoms.cc
  smiles.cc
  smiles_support.cc
  sort_atoms.cc
  ss_atom_env.cc
  ss_bonds.cc
  ss_ring.cc
  ss_ring_base.cc
  ss_ring_sys.cc
  standardise.cc
  substructure_a.cc
  substructure_chiral.cc
  substructure_env.cc
  substructure_nmab.cc
  substructure_query_match.cc
  substructure_results.cc
  substructure_spec.cc
  symm_class_can_rank.cc
  symmetry.cc
  target.cc
  temp_detach_atoms.cc
  toggle_kekule_form.cc
  tokenise_atomic_smarts.cc
  u3b.c
  unique.cc
)

target_include_directories(iwmolecule
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(iwmolecule
  PUBLIC
    atom_type_ext_cc_proto
    donor_acceptor_cc_proto
    etrans_cc_proto
    geometric_constraints_cc_proto
    mol2graph_cc_proto
    molecule_to_query_cc_proto
    pharmacophore_cc_proto
    space_vector
    standardise_cc_proto
    substructure_cc_proto
    toggle_kekule_form_cc_proto
    accumulator
    iwcmdline
    iwstring_data_source
    iwtfdata_record
    iw_tdt
    iwbits
    iwmisc
    iwconfig
    iwre2
    iwqsort
    iwstring
    xmlparser
    protobuf::libprotobuf
)

add_library(moleculeio
    charmm.cc
    cif.cc
    csv.cc
    marvin.cc
    mdl.cc
    mdl_v30.cc
    mmod.cc
    moe.cc
    mrk.cc
    msi.cc
    pdb.cc
    rwmolecule.cc
    smi.cc
    tripos.cc
    xyz.cc
)

target_link_libraries(moleculeio
  PUBLIC
    geometric_constraints_cc_proto
    substructure_cc_proto
    molecule_to_query_cc_proto
    reaction_cc_proto
    standardise_cc_proto
    toggle_kekule_form_cc_proto
    iwmolecule
)

target_include_directories(moleculeio
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

# Logical molecule IO target.
# Currently moleculeio.cc is compiled as part of iwmolecule,
# matching Bazel source ownership.
