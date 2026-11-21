#include <iostream>

#define ACCUMULATOR_IMPLEMENTATION
#include "accumulator.h"

template class Accumulator<float>;
template class Accumulator_Base<float, KahanSum>;

#ifdef __GNUG__

template std::ostream & operator << (std::ostream &, const Accumulator<float> &);

#else

template std::ostream & operator << (std::ostream &, const Accumulator<float> &);

#endif
