#include <algorithm>
#include <functional>
#include <stdlib.h>

#include "bsquared.h"

#include "Foundational/accumulator/accumulator.h"

template <typename T>
int
ComputeBSquared(const T * predicted,
                  T * ybest,
                  int n,
                  double & Bsquared)
{
  assert (n > 1);

  Bsquared = 0.0;

  Accumulator<double> acc;

  for (int i = 0; i < n; i++)
  {
    ybest[i] = predicted[i];

    acc.extra (predicted[i]);
  }

  double ymean = acc.average ();

  std::sort(ybest, ybest + n, std::greater<T>());

  double predicted_csum = 0.0;
  double ybest_csum = 0.0;
  double yworst_csum = 0.0;

  double last_non_zero_weight = 0.0;

  double sumwt = 0.0;

  for (int i = 0; i < n; i++)
  {
    double yworst = ybest[n - i - 1];

    double i_plus_one = static_cast<double> (i + 1);

    predicted_csum += predicted[i];
    double y2 = predicted_csum / i_plus_one;

    ybest_csum += ybest[i];
    double y2best = ybest_csum / i_plus_one;

    yworst_csum += yworst;
    double y2worst = yworst_csum / i_plus_one;

    double weight = y2best - y2worst;

    if (weight > 0.0)
      last_non_zero_weight = weight;
    else
      weight = last_non_zero_weight;

    if (y2 == ymean)
      continue;

    double myweight;
    if (y2 > ymean)
      myweight = y2best - ymean;
    else
      myweight = ymean - y2worst;

    if (myweight < 1.0e-05)
      continue;

    Bsquared += (y2 - ymean) / myweight * weight;

    sumwt += weight;
  }

  if (sumwt != 0.0)
    Bsquared = Bsquared / sumwt;

  return 1;
}

template <typename T>
int
ComputeBSquared (const T * predicted,
                   T * ybest,
                   int n,
                   double * Bsquared_array)
{
  assert (n > 1);

  Accumulator<double> acc;

  for (int i = 0; i < n; i++)
  {
    ybest[i] = predicted[i];

    acc.extra (predicted[i]);
  }

  double ymean = acc.average ();

  std::sort(ybest, ybest + n, std::greater<T>());

  double range = ybest[0] - ybest[n - 1];     // the dynamic range of the data

  assert (range >= 0.0);
  if (0.0 == range)
    range = 1.0;

  double too_small = range * 1.0e-05;

  double predicted_csum = 0.0;
  double ybest_csum = 0.0;
  double yworst_csum = 0.0;

  double last_non_zero_weight = 0.0;

  double sumwt = 0.0;
  double Bsquared = 0.0;

  for (int i = 0; i < n; i++)
  {
    double yworst = ybest[n - i - 1];

    double i_plus_one = static_cast<double> (i + 1);

    predicted_csum += predicted[i];
    double y2 = predicted_csum / i_plus_one;

    ybest_csum += ybest[i];
    double y2best = ybest_csum / i_plus_one;

    yworst_csum += yworst;
    double y2worst = yworst_csum / i_plus_one;

    double weight = y2best - y2worst;

    if (weight > 0.0)
      last_non_zero_weight = weight;
    else
      weight = last_non_zero_weight;

    sumwt += weight;

    double myweight;
    if (y2 > ymean)
      myweight = y2best - ymean;
    else if (y2 < ymean)
      myweight = ymean - y2worst;
    else
      myweight = 0.0;

    if (myweight > too_small && sumwt > 0.0)
    {
      Bsquared += (y2 - ymean) / myweight * weight;

      Bsquared_array[i] = Bsquared / sumwt;
    }
    else if (i > 0)
    {
      Bsquared_array[i] = Bsquared_array[i - 1];
    }
    else     // i == 0 and no contribution! let's hope this never happens
    {
      Bsquared_array[i] = 0.0;
    }
  }

  return 1;
}

int compute_b_squared(const float* values, float* work, int n, double& result) {
  return ComputeBSquared(values, work, n, result);
}
int compute_b_squared(const float* values, float* work, int n, double* result) {
  return ComputeBSquared(values, work, n, result);
}
int compute_b_squared(const double* values, double* work, int n, double& result) {
  return ComputeBSquared(values, work, n, result);
}
int compute_b_squared(const double* values, double* work, int n, double* result) {
  return ComputeBSquared(values, work, n, result);
}
