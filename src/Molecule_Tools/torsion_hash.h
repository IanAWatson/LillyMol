#ifndef TORSION_HASH_H
#define TORSION_HASH_H

#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iwstring_and_file_descriptor.h"

class Torsion_Hash : public IW_STL_Hash_Map_int
{
  private:

  public:
    int torsions_stored () const { return size ();}

    int write (std::ostream &) const;
    int do_write (iwstring::IWString_and_File_Descriptor &) const;

    int add (const IWString &);

    int unique_identifier (const IWString &);

    int fetch_unique_identifier (const IWString &, int &) const;
};

#endif
