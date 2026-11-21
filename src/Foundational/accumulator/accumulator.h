#ifndef FOUNDATIONAL_ACCUMULATOR_ACCUMULATOR_H_
#define FOUNDATIONAL_ACCUMULATOR_ACCUMULATOR_H_

#include <iostream>
#include <type_traits>

#include "Foundational/iwstring/iwstring.h"

#include "kahan_sum.h"
#include "running_variance.h"

#ifdef BUILD_BAZEL
#include "Foundational/accumulator/accumulator.pb.h"
#else
#include "accumulator.pb.h"
#endif

namespace iwaccumulator {
template <>
class VarianceState<KahanSum> {
 protected:
  RunningVariance _variance;
};
}  // namespace iwaccumulator

template <typename T, typename SUMMER>
class Accumulator_Base : private iwaccumulator::VarianceState<SUMMER>
{
  static_assert(std::is_same_v<SUMMER, KahanSum> ||
                    (std::is_integral_v<T> && std::is_integral_v<SUMMER>),
                "Non-Kahan accumulators require integral value and summer types");

  protected:
    uint64_t _n;
    T _minval;
    T _maxval;
    SUMMER _xsum;
    SUMMER _x2sum;

//  protected functions

    void _reset();
    void _default_values ();

  public:
    Accumulator_Base ();

    int ok () const;

    uint64_t n () const { return _n;}

    bool empty() const {
      return _n == 0;
    }

    T minval () const { return _minval;}
    T maxval () const { return _maxval;}
    T range  () const { return _maxval - _minval;}

    void reset ();

    T sum () const { return _xsum;}
    T sum_of_squares () const { return _x2sum;}

    unsigned int extra (T);
    unsigned int extra (T, int);       // add N copies; non-positive N is a no-op
    Accumulator_Base<T, SUMMER> & operator = (const Accumulator_Base<T, SUMMER> &);
    void operator += (const Accumulator_Base<T, SUMMER> & rhs) { extra (rhs);}

    unsigned int extra (const T *, int);       // add an array of values
    unsigned int extra (const Accumulator_Base<T, SUMMER> &);

    void operator () (T e) { (void) extra(e);}

    int subtract_data (const Accumulator_Base<T, SUMMER> &, Accumulator_Base<T, SUMMER> &) const;

    double average () const;
    int    average (double &);

    double average_if_available_minval_if_not () const;

    // Sample variance. Floating accumulators maintain central moments in double
    // precision; their public sums and proto fields retain their existing form.
    double variance () const;
    int    variance (double &);
};

template <typename T>
class Accumulator : public Accumulator_Base<T, KahanSum>
{
#ifdef IW_TWO_PHASE_TEMPLATES
  protected:
    using Accumulator_Base<T, KahanSum>::_n;
    using Accumulator_Base<T, KahanSum>::_minval;
    using Accumulator_Base<T, KahanSum>::_maxval;
    using Accumulator_Base<T, KahanSum>::_xsum;
    using Accumulator_Base<T, KahanSum>::_x2sum;
#endif

  private:
//  double _xsum;
//  double _x2sum;

//  private functions

    void _default_values ();

  public:

    int BuildProto(accumulator::AccumulatorData& proto) const;

};

template <typename T>
class Accumulator_Int : public Accumulator_Base<T, T>
{
  static_assert(std::is_integral_v<T>, "Accumulator_Int requires an integral type");

#ifdef IW_TWO_PHASE_TEMPLATES
  protected:
    using Accumulator_Base<T, T>::_n;
    using Accumulator_Base<T, T>::_minval;
    using Accumulator_Base<T, T>::_maxval;
    using Accumulator_Base<T, T>::_xsum;
    using Accumulator_Base<T, T>::_x2sum;
#endif

  private:

  public:
};

template <typename T>
std::ostream & operator << (std::ostream &, const Accumulator<T> &);
template <typename T>
std::ostream & operator << (std::ostream &, const Accumulator_Int<T> &);

template <typename T>
class Accumulator_with_Missing_Values : public Accumulator<T>
{
  private:
    int _nmissing;

  public:
    Accumulator_with_Missing_Values ();

    void extra_missing_value () { _nmissing++;}

    int number_missing_values () const { return _nmissing;}

    void reset ();
};

#if defined(ACCUMULATOR_IMPLEMENTATION) || defined(IW_IMPLEMENTATIONS_EXPOSED)

#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <iomanip>

template <typename T, typename SUMMER>
void
Accumulator_Base<T, SUMMER>::_default_values ()
{
  _n = 0;
  _minval = static_cast<T>(0);
  _maxval = static_cast<T>(0);

  _xsum = static_cast<T>(0);
  _x2sum = static_cast<T>(0);
  if constexpr (std::is_same_v<SUMMER, KahanSum>) {
    this->_variance = iwaccumulator::RunningVariance();
  }

  return;
}

template <typename T, typename SUMMER>
Accumulator_Base<T, SUMMER>::Accumulator_Base ()
{
  _default_values();

  return;
}

template <typename T>
void
Accumulator<T>::_default_values ()
{
  Accumulator_Base<T, KahanSum>::_default_values();

  return;
}

template <typename T, typename SUMMER>
int
Accumulator_Base<T, SUMMER>::ok() const
{
  return 1;
}

template <typename T, typename SUMMER>
void
Accumulator_Base<T, SUMMER>::reset ()
{
  _default_values();

  return;
}

template <typename T, typename SUMMER>
unsigned int
Accumulator_Base<T, SUMMER>::extra (T x)
{
  if constexpr (std::is_same_v<SUMMER, KahanSum>) {
    const double value = static_cast<double>(x);
    this->_variance.Add(value, _n);
    _x2sum += value * value;
  } else {
    _x2sum += x * x;
  }
  _n++;
  _xsum += x;

  if (1 == _n)
  {
    _minval = x;
    _maxval = x;
  }
  else if (x < _minval)
  {
    _minval = x;
  }
  else if (x > _maxval)
  {
    _maxval = x;
  }

  return _n;
}

/*
  Add N copies of a value
*/

template <typename T, typename SUMMER>
unsigned int
Accumulator_Base<T, SUMMER>::extra (T x, int n)
{
  if (n <= 0) return _n;

  if (0 == _n)
  {
    _minval = x;
    _maxval = x;
  }
  else if (x < _minval)
  {
    _minval = x;
  }
  else if (x > _maxval)
  {
    _maxval = x;
  }

  if constexpr (std::is_same_v<SUMMER, KahanSum>) {
    const double value = static_cast<double>(x);
    this->_variance.AddCopies(value, n, _n);
    _xsum += static_cast<double>(n) * value;
    _x2sum += static_cast<double>(n) * value * value;
  } else {
    _xsum += n * x;
    _x2sum += n * (x * x);
  }
  _n += n;

  return _n;
}

template <typename T, typename SUMMER>
unsigned int
Accumulator_Base<T, SUMMER>::extra (const T * x, int n)
{
  for (int i = 0; i < n; ++i) {
    extra(x[i]);
  }

  return _n;
}

template <typename T, typename SUMMER>
unsigned int
Accumulator_Base<T, SUMMER>::extra (const Accumulator_Base<T, SUMMER> & rhs)
{
  if (0 == rhs._n)
    return _n;

  if constexpr (std::is_same_v<SUMMER, KahanSum>) {
    this->_variance.Merge(rhs._variance, _n, rhs._n);
  }

  if (0 == _n)
  {
    _n = rhs._n;
    _xsum  = rhs._xsum;
    _x2sum = rhs._x2sum;
    _minval = rhs._minval;
    _maxval = rhs._maxval;

    return _n;
  }

  _n += rhs._n;
  _xsum += rhs._xsum;
  _x2sum += rhs._x2sum;

  if (rhs._minval < _minval)
    _minval = rhs._minval;
  if (rhs._maxval > _maxval)
    _maxval = rhs._maxval;

  return _n;
}

/*
  Strictly speaking, _minval and _maxval are unknown, but we lie...
*/

template <typename T, typename SUMMER>
int
Accumulator_Base<T, SUMMER>::subtract_data (const Accumulator_Base<T, SUMMER> & rhs,
                                            Accumulator_Base<T, SUMMER> & zresult) const
{
  assert (_n >= rhs._n);

  Accumulator_Base<T, SUMMER> result;
  if (_n == rhs._n) {
    zresult = result;
    return 1;
  }
  if constexpr (std::is_same_v<SUMMER, KahanSum>) {
    result._variance = this->_variance.Without(rhs._variance, _n, rhs._n);
  }
  result._n     = _n - rhs._n;
  result._xsum  = _xsum - rhs._xsum;
  result._x2sum = _x2sum - rhs._x2sum;

  if (result._x2sum >= 0.0)    // exactly how things should be
    ;
  else if (result._x2sum > -1.0e-06)     // probably just roundoff errors
  {
    std::cerr << "Accumulator::subtract_data: trimming possible roundoff error " << result._x2sum << '\n';
    result._x2sum = 0.0;
  }
  else
  {
    std::cerr << "Accumulator::subtract_data: invalid x2sum\n";
    abort ();
  }

  if (_minval < rhs._minval)
    result._minval = _minval;
  else
    result._minval = rhs._minval;

  if (_maxval > rhs._maxval)
    result._maxval = _maxval;
  else
    result._maxval = rhs._maxval;

  zresult = result;

  return 1;
}

template <typename T, typename SUMMER>
double
Accumulator_Base<T, SUMMER>::average () const
{
  return static_cast<double>(_xsum) / static_cast<double>(_n);
}

template <typename T, typename SUMMER>
double
Accumulator_Base<T, SUMMER>::variance () const
{       
  assert (_n > 1);

  if (_n < 2) {
    return 0;
  }

  if constexpr (std::is_same_v<SUMMER, KahanSum>) {
    return this->_variance.Variance(_n);
  }

  double tave = Accumulator_Base<T, SUMMER>::average();
  double rc = _x2sum - _n * tave * tave;

  if (rc < 0.0)   // presumably some roundoff
  {
//  std::cerr << "Accumulator::variance: Warning, negative variance intermediate " << tmp << '\n';

    return 0.0;
  }

  rc = rc / static_cast<double>(_n - 1);

  return rc;
}       

template <typename T, typename SUMMER>
int
Accumulator_Base<T, SUMMER>::variance (double & v)
{
  if (_n < 2)
    return 0;

  v = variance();

  return 1;
}

template <typename T>
std::ostream &
operator << (std::ostream & os, const Accumulator<T> & ac)
{
  assert (ac.n() > 0);

  os << "Accumulator " << ac.n() << " values, average " << std::setw(8) << ac.average();
  if (ac.n() > 1)
    os << ", variance " << std::setw(8) << ac.variance();
  os << "\n";

  return os << "min = " << ac.minval() << " max = " << ac.maxval();
}

template <typename T>
std::ostream &
operator << (std::ostream & os, const Accumulator_Int<T> & ac)
{
  assert (ac.n() > 0);

  os << "Accumulator_Int " << ac.n() << " values, average " << ac.average();
  if (ac.n() > 1)
    os << ", variance " << ac.variance();
  os << "\n";

  return os << "min = " << ac.minval() << " max = " << ac.maxval();
}

template <typename T, typename SUMMER>
Accumulator_Base<T, SUMMER> &
Accumulator_Base<T, SUMMER>::operator = (const Accumulator_Base<T, SUMMER> & rhs)
{
  _n = rhs._n;

  _xsum = rhs._xsum;
  _x2sum = rhs._x2sum;
  _minval = rhs._minval;
  _maxval = rhs._maxval;
  if constexpr (std::is_same_v<SUMMER, KahanSum>) {
    this->_variance = rhs._variance;
  }

  return *this;
}

template <typename T, typename SUMMER>
double
Accumulator_Base<T, SUMMER>::average_if_available_minval_if_not () const
{
  if (_n > 1)
    return average();

  if (_n > 0)
    return _minval;

  std::cerr << "Accumulator_Base::average_if_available_minval_if_not: no data!\n";
  return 0.0;
}

template <typename T>
int
Accumulator<T>::BuildProto(accumulator::AccumulatorData& proto) const {
  proto.set_n(_n);
  proto.set_minval(_minval);
  proto.set_maxval(_maxval);
  proto.set_xsum(_xsum);
  proto.set_x2sum(_x2sum);

  return 1;
}

#endif

#ifdef ACCUMULATOR_W_MISSING_IMPLEMENTATION

template <typename T>
Accumulator_with_Missing_Values<T>::Accumulator_with_Missing_Values ()
{
  _nmissing = 0;

  return;
}

template <typename T>
void
Accumulator_with_Missing_Values<T>::reset ()
{
  Accumulator<T>::reset();

  _nmissing = 0;

  return;
}

#endif

#ifdef ACCUMULATOR_INT_IMPLEMENTATION
#endif

#endif  // FOUNDATIONAL_ACCUMULATOR_ACCUMULATOR_H_
