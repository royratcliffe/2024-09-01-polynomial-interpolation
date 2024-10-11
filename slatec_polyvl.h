/*!
 * \file slatec_polyvl.h
 * \copyright Roy Ratcliffe, Northumberland, United Kingdom
 *
 * SPDX-License-Identifier: MIT
 *
 * Permission is hereby granted, free of charge,  to any person obtaining a
 * copy  of  this  software  and    associated   documentation  files  (the
 * "Software"), to deal in  the   Software  without  restriction, including
 * without limitation the rights to  use,   copy,  modify,  merge, publish,
 * distribute, sublicense, and/or sell  copies  of   the  Software,  and to
 * permit persons to whom the Software is   furnished  to do so, subject to
 * the following conditions:
 *
 *     The above copyright notice and this permission notice shall be
 *     included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT  WARRANTY OF ANY KIND, EXPRESS
 * OR  IMPLIED,  INCLUDING  BUT  NOT   LIMITED    TO   THE   WARRANTIES  OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR   PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS  OR   COPYRIGHT  HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY,  WHETHER   IN  AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM,  OUT  OF   OR  IN  CONNECTION  WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#pragma once

/*
 * for size_t
 */
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

enum slatec_polyvl_status { slatec_polyvl_success, slatec_polyvl_failure = -1 };

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
static inline enum slatec_polyvl_status slatec_polyvl(double xx, double *yy,
                                                      size_t n,
                                                      const double x[],
                                                      const double c[]) {
  if (n == 0)
    return slatec_polyvl_failure;
  double pione = 1, pone = c[0];
  if (n == 1) {
    *yy = pone;
    return slatec_polyvl_success;
  }
  double ptwo;
  for (size_t k = 1; k < n; k++) {
    const double pitwo = (xx - x[k - 1]) * pione;
    pione = pitwo;
    ptwo = pone + pitwo * c[k];
    pone = ptwo;
  }
  *yy = ptwo;
  return slatec_polyvl_success;
}

/*!
 * \brief Computes polynomial single-precision values.
 */
static inline enum slatec_polyvl_status slatec_polyvlf(float xx, float *yy,
                                                       size_t n,
                                                       const float x[],
                                                       const float c[]) {
  if (n == 0)
    return slatec_polyvl_failure;
  float pione = 1, pone = c[0];
  if (n == 1) {
    *yy = pone;
    return slatec_polyvl_success;
  }
  float ptwo;
  for (size_t k = 1; k < n; k++) {
    const float pitwo = (xx - x[k - 1]) * pione;
    pione = pitwo;
    ptwo = pone + pitwo * c[k];
    pone = ptwo;
  }
  *yy = ptwo;
  return slatec_polyvl_success;
}

#ifdef __cplusplus
}
#endif
