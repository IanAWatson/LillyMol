#ifndef IW_ACCUMULATOR_V2H
#define IW_ACCUMULATOR_V2H

#include <iostream>

#include "Foundational/iwstring/iwstring.h"

#include "kahan_sum.h"
#include "running_variance.h"

template <typename COUNTER, typename T>
class Accumulator_V2
{
  private:
    COUNTER _n;
//  double _xsum;
//  double _x2sum;
    KahanSum _xsum;
    KahanSum _x2sum;
    iwaccumulator::RunningVariance _variance;
    T _minval;
    T _maxval;

//  private functions

    void _default_values ();

  public:
    Accumulator_V2 ();

    int ok () const;

    Accumulator_V2<COUNTER, T> & operator = (const Accumulator_V2<COUNTER, T> &);
    void operator += (const Accumulator_V2<COUNTER, T> & rhs) { extra (rhs);}

    COUNTER extra (T);
    COUNTER extra (T, int);       // add N copies of a value
    COUNTER extra (const T *, int);       // add an array of values
    COUNTER extra (const Accumulator_V2<COUNTER, T> &);

    COUNTER n () const { return _n;}

    double sum () const { return _xsum;}
    double sum_of_squares () const { return _x2sum;}

    double average () const;
    int    average (double &);
    double variance ();
    double variance () const;
    int    variance (double &);

    T minval () const { return _minval;}
    T maxval () const { return _maxval;}
    T range  () const { return _maxval - _minval;}

    void reset ();

    double average_if_available_minval_if_not () const;

    int subtract_data (const Accumulator_V2<COUNTER, T> &, Accumulator_V2<COUNTER, T> &) const;
};

template <typename COUNTER, typename T>
std::ostream & operator << (std::ostream &, const Accumulator_V2<COUNTER, T> &);

template <typename COUNTER, typename T>
class Accumulator_with_Missing_Values_V2 : public Accumulator_V2<COUNTER, T>
{
  private:
    COUNTER _nmissing;

  public:
    Accumulator_with_Missing_Values_V2 ();

    void extra_missing_value () { _nmissing++;}

    int number_missing_values () const { return _nmissing;}

    void reset ();
};

#ifdef ACCUMULATOR_V2_IMPLEMENTATION

#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <iomanip>

template <typename COUNTER, typename T>
void
Accumulator_V2<COUNTER, T>::_default_values ()
{
  _n = 0;
  _xsum = 0.0;
  _x2sum = 0.0;
  _variance = iwaccumulator::RunningVariance();
  _minval = 0;
  _maxval = 0;
}

template <typename COUNTER, typename T>
Accumulator_V2<COUNTER, T>::Accumulator_V2 ()
{
  _default_values ();

  return;
}

template <typename COUNTER, typename T>
int
Accumulator_V2<COUNTER, T>::ok () const
{
  return 1;
}

template <typename COUNTER, typename T>
void
Accumulator_V2<COUNTER, T>::reset ()
{
  _default_values ();

  return;
}

template <typename COUNTER, typename T>
COUNTER
Accumulator_V2<COUNTER, T>::extra (T x)
{
  const double value = static_cast<double>(x);
  _variance.Add(value, _n);
  _n++;
  _xsum += value;
  _x2sum += value * value;

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

template <typename COUNTER, typename T>
COUNTER
Accumulator_V2<COUNTER, T>::extra (T x, int n)
{
  if (n <= 0) return _n;
  const double value = static_cast<double>(x);
  _variance.AddCopies(value, n, _n);
  _xsum += static_cast<double>(n) * value;
  _x2sum += static_cast<double>(n) * value * value;

  if (_n == 0) {
    _minval = x;
    _maxval = x;
  } else if (x < _minval) {
    _minval = x;
  } else if (x > _maxval) {
    _maxval = x;
  }
  _n += n;

  return _n;
}

template <typename COUNTER, typename T>
COUNTER
Accumulator_V2<COUNTER, T>::extra (const T * x, int n)
{
  for (int i = 0; i < n; ++i) {
    extra(x[i]);
  }

  return _n;
}

template <typename COUNTER, typename T>
COUNTER
Accumulator_V2<COUNTER, T>::extra (const Accumulator_V2<COUNTER, T> & rhs)
{
  if (0 == rhs._n)
    return _n;

  _variance.Merge(rhs._variance, _n, rhs._n);

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

template <typename COUNTER, typename T>
int
Accumulator_V2<COUNTER, T>::subtract_data (const Accumulator_V2<COUNTER, T> & rhs,
                               Accumulator_V2<COUNTER, T> & zresult) const
{
  assert (_n >= rhs._n);

  Accumulator_V2<COUNTER, T> result;
  if (_n == rhs._n) {
    zresult = result;
    return 1;
  }
  result._variance = _variance.Without(rhs._variance, _n, rhs._n);
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

template <typename COUNTER, typename T>
double
Accumulator_V2<COUNTER, T>::average () const
{
  return static_cast<double> (_xsum) / static_cast<double> (_n);
}

template <typename COUNTER, typename T>
double
Accumulator_V2<COUNTER, T>::variance ()
{       
  assert (_n > 1);

  return _variance.Variance(_n);
}       

template <typename COUNTER, typename T>
double
Accumulator_V2<COUNTER, T>::variance () const
{       
  assert (_n > 1);

  return _variance.Variance(_n);
}       

template <typename COUNTER, typename T>
int
Accumulator_V2<COUNTER, T>::variance (double & v)
{
  if (_n < 2)
    return 0;

  v = variance ();

  return 1;
}

template <typename COUNTER, typename T>
std::ostream &
operator << (std::ostream & os, const Accumulator_V2<COUNTER, T> & ac)
{
  assert (ac.n () > 0);
  using std::setw;

  os << "Accumulator " << ac.n () << " values, average " << setw(8) << ac.average ();
  if (ac.n() > 1)
    os << ", variance " << setw (8) << ac.variance ();
  os << "\n";
  return os << "min = " << ac.minval () << " max = " << ac.maxval ();
}

template <typename COUNTER, typename T>
Accumulator_V2<COUNTER, T> &
Accumulator_V2<COUNTER, T>::operator = (const Accumulator_V2<COUNTER, T> & rhs)
{
  _n = rhs._n;

  _xsum = rhs._xsum;
  _x2sum = rhs._x2sum;
  _minval = rhs._minval;
  _maxval = rhs._maxval;
  _variance = rhs._variance;

  return *this;
}

template <typename COUNTER, typename T>
double
Accumulator_V2<COUNTER, T>::average_if_available_minval_if_not () const
{
  if (_n > 1)
    return average ();

  if (_n > 0)
    return _minval;

  std::cerr << "Accumulator::average_if_available_minval_if_not: no data!\n";
  return 0.0;
}

#endif

#ifdef ACCUMULATOR_W_MISSING_IMPLEMENTATION

template <typename COUNTER, typename T>
Accumulator_with_Missing_Values_V2<COUNTER, T>::Accumulator_with_Missing_Values_V2 ()
{
  _nmissing = 0;

  return;
}

template <typename COUNTER, typename T>
void
Accumulator_with_Missing_Values_V2<COUNTER, T>::reset ()
{
  Accumulator_V2<COUNTER, T>::reset ();

  _nmissing = 0;

  return;
}

#endif

#endif
