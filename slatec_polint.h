#pragma once

/*
 * for size_t
 */
#include <stddef.h>

enum slatec_polint_status {
  slatec_polint_success,
  slatec_polint_failure = -1,
  slatec_polint_abscissae_not_distinct = -2
};

/*!
 * \brief Double-precision polynomial interpolation.
 *
 * \details Function \c slatec_polint generates a polynomial that interpolates
 * the vectors \c x[i] by \c y[i] where \c i ranges from 1 to \c{n}. It prepares
 * information in the array \c c that can be utilized by the subroutine \c
 * polyvl to calculate the polynomial and its derivatives.
 *
 * \see https://netlib.org/slatec/src/polint.f
 */
static inline enum slatec_polint_status
slatec_polint(size_t n, const double x[], const double y[], double c[]) {
  if (n == 0)
    return slatec_polint_failure;
  c[0] = y[0];
  if (n == 1)
    return slatec_polint_success;
  for (size_t k = 1; k < n; k++) {
    c[k] = y[k];
    for (size_t i = 0; i < k; i++) {
      const double dif = x[i] - x[k];
      if (dif == 0)
        return slatec_polint_abscissae_not_distinct;
      c[k] = (c[i] - c[k]) / dif;
    }
  }
  return slatec_polint_success;
}

/*!
 * \brief Single-precision polynomial interpolation.
 */
static inline enum slatec_polint_status
slatec_polintf(size_t n, const float x[], const float y[], float c[]) {
  if (n == 0)
    return slatec_polint_failure;
  c[0] = y[0];
  if (n == 1)
    return slatec_polint_success;
  for (size_t k = 1; k < n; k++) {
    c[k] = y[k];
    for (size_t i = 0; i < k; i++) {
      const float dif = x[i] - x[k];
      if (dif == 0)
        return slatec_polint_abscissae_not_distinct;
      c[k] = (c[i] - c[k]) / dif;
    }
  }
  return slatec_polint_success;
}
