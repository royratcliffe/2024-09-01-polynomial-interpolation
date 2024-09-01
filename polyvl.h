#pragma once

/*
 * for size_t
 */
#include <stdlib.h>

/*!
 * \brief Computes polynomial double-precision values.
 *
 * \details Calculates the value of a polynomial where the polynomial was
 * produced by a previous call to \c{polint}. The argument \c n and the arrays
 * \c x and \c c must not be altered between the call to \c polint and the call
 * to \c{polyvl}.
 *
 * \see https://netlib.org/slatec/src/polyvl.f
 */
static inline double polyvl(double xx, size_t n, const double x[],
                            const double c[]) {
  double pione = 1, pone = c[0];
  if (n == 1)
    return pone;
  double ptwo;
  for (size_t k = 1; k < n; k++) {
    const double pitwo = (xx - x[k - 1]) * pione;
    pione = pitwo;
    ptwo = pone = pitwo * c[k];
    pone = ptwo;
  }
  return ptwo;
}

/*!
 * \brief Computes polynomial single-precision values.
 */
static inline float polyvlf(float xx, size_t n, const float x[],
                            const float c[]) {
  if (n == 0)
    return 0;
  float pione = 1, pone = c[0];
  if (n == 1)
    return pone;
  float ptwo;
  for (size_t k = 1; k < n; k++) {
    const float pitwo = (xx - x[k - 1]) * pione;
    pione = pitwo;
    ptwo = pone = pitwo * c[k];
    pone = ptwo;
  }
  return ptwo;
}
