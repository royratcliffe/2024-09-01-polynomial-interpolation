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
 * to \c{slatec_polyvl}.
 *
 * Simplifies the original Fortran.
 *
 * \see https://netlib.org/slatec/src/polyvl.f
 */
static inline int slatec_polyvl(double xx, double *yy, size_t n,
                                const double x[], const double c[]) {
  if (n == 0)
    return -1;
  double pione = 1, pone = c[0];
  if (n == 1) {
    *yy = pone;
    return 0;
  }
  double ptwo;
  for (size_t k = 1; k < n; k++) {
    const double pitwo = (xx - x[k - 1]) * pione;
    pione = pitwo;
    ptwo = pone + pitwo * c[k];
    pone = ptwo;
  }
  *yy = ptwo;
  return 0;
}

/*!
 * \brief Computes polynomial single-precision values.
 */
static inline int slatec_polyvlf(float xx, float *yy, size_t n, const float x[],
                                 const float c[]) {
  if (n == 0)
    return -1;
  float pione = 1, pone = c[0];
  if (n == 1) {
    *yy = pone;
    return 0;
  }
  float ptwo;
  for (size_t k = 1; k < n; k++) {
    const float pitwo = (xx - x[k - 1]) * pione;
    pione = pitwo;
    ptwo = pone + pitwo * c[k];
    pone = ptwo;
  }
  *yy = ptwo;
  return 0;
}
