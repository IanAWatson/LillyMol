#ifndef TVERSKY_H
#define TVERSKY_H

class Command_Line;

/*
  Need typedef for similarity_type_t so include dyfp.h
*/

#include "dyfp.h"

class Tversky
{
  private:
    float _a;
    float _b;

    int _active;

//  One possibility is to have the Tversky object function in an "optimistic"
//  mode. It computes the distances either way (swapping A and B), as well as
//  Tanimoto, and then return the shortest distance.

    int _optimistic_mode;

//  When we get to counted fingerprints, there are two possibilities for Tversky.
//  Treat the fingerprint equivalent to 0/1 or do counts

    int _treat_non_colliding_as_01;

//  Jun 2002. Fixed bugs that occurred with zero bits set and tversky metric. o
//  But what to do when one FP has 0 bits set and the other has some bits set.
//  I'm experimenting with different similarity values depending on nset

    int _nset_sensitive_zero_bit_similarity;

  public:
    Tversky ();

    int optimistic_mode () const { return _optimistic_mode;}

    float a () const { return _a;}
    float b () const { return _b;}

    int active () const { return _active;}

//  Set the parameters directly, without a command line. The weights must be
//  non-negative, and not both zero. If `scale_to_2` the values are rescaled so
//  that a+b=2, as happens with the -V option, so a=b=1 is Tanimoto and only the
//  ratio matters. Returns 0 if the values are invalid.
    int Set (float a, float b, bool scale_to_2 = true);

    int parse_command_line (Command_Line &, char, char, int);
    int parse_command_line (Command_Line &, char, int);

    int treat_non_colliding_as_01 () const { return _treat_non_colliding_as_01;}

    int nset_sensitive_zero_bit_similarity () const { return _nset_sensitive_zero_bit_similarity;}

    similarity_type_t tanimoto (int nset1, int nset2, int bic) const;
};

#include <iostream>

extern int display_standard_tversky_options (std::ostream &);
extern int display_standard_tversky_options (std::ostream &, char);

#endif
