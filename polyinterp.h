// -*- c++ -*- $Id: polyinterp.h,v 1.2 2000/11/01 11:19:05 royr Exp $
// SPDX-License-Identifier: MIT
//
// A polynomial interpolator.
//
// Operation add(x,y) adds point (x,y) to the interpolating polynomial.
// All the abscissae must be unique.  No two abscissae can be alike.
// How close can they be?  A threshold specifies the minimum distance.
// Points with abscissae closer than the minimum threshold merge
// at the arithmetic mean.

extern "C" {
#include "slatec_polint.h"
#include "slatec_polyvl.h"
}

#ifdef __cplusplus

#include <functional>
#include <limits>
#include <vector>

template <typename Scalar>
enum slatec_polint_status polint(size_t n, const Scalar x[], const Scalar y[],
                                 Scalar c[]);

template <typename Scalar>
enum slatec_polyvl_status polyvl(Scalar xx, Scalar *yy, size_t n,
                                 const Scalar x[], const Scalar c[]);

template <typename Scalar>
struct poly_interpolator // a unary functor
{
private:
  Scalar abscissaDeltaThres;
  std::vector<Scalar> X, Y, C;
  std::vector<int> N;

public:
  poly_interpolator() : abscissaDeltaThres(0) {}

  void set_abscissa_thres(Scalar const &x) {
    if (0 <= x)
      abscissaDeltaThres = x;
  }

  void add(Scalar const &x, Scalar const &y) {
    // sort x by insertion -- iterate while X[i] < x
    auto Xi = X.begin();
    for (; Xi != X.end() && *Xi < x; ++Xi)
      ;
    // Xi == X.end() || *Xi >= x
    auto i = std::distance(X.begin(), Xi);
    if (Xi != X.begin() && x - Xi[-1] <= abscissaDeltaThres) {
      --i;
      X[i] = (x + X[i] * N[i]) / (N[i] + 1);
      Y[i] = (y + Y[i] * N[i]) / (N[i] + 1);
      ++N[i];
    } else if (Xi != X.end() && Xi[0] - x <= abscissaDeltaThres) {
      X[i] = (x + X[i] * N[i]) / (N[i] + 1);
      Y[i] = (y + Y[i] * N[i]) / (N[i] + 1);
      ++N[i];
    } else {
      // Get the dangerous bit over with!  Throwing is the worry.
      // It's fine---except half way through an add op.  Since there
      // are four separate vectors to update.  Four memory
      // re-allocations aren't atomic!
      X.reserve(N.size() + 1);
      Y.reserve(N.size() + 1);
      C.reserve(N.size() + 1);
      N.reserve(N.size() + 1);

      Xi = X.begin();
      auto Yi = Y.begin();
      auto Ci = C.begin();
      auto Ni = N.begin();
      std::advance(Xi, i);
      std::advance(Yi, i);
      std::advance(Ci, i);
      std::advance(Ni, i);
      X.insert(Xi, x);
      Y.insert(Yi, y);
      C.insert(Ci, 0);
      N.insert(Ni, 1);
    }
  }

  void interpolate() {
    enum slatec_polint_status status =
        polint(N.size(), X.data(), Y.data(), C.data());
    if (status != slatec_polint_success)
      throw status;
  }

  Scalar operator()(const Scalar &x) const {
    if (N.size() == 0)
      return x;
    Scalar y;
    enum slatec_polyvl_status status =
        polyvl(x, &y, N.size(), X.data(), C.data());
    if (status != slatec_polyvl_success)
      throw status;
    return y;
  }

  size_t n() const   // How many interpolating
  {                  // points evaluate the
    return N.size(); // polynomial?
  }

  void clear() {
    X.clear();
    Y.clear();
    C.clear();
    N.clear();
  }
};

////////////////////////////////////////////////////////////////////////

#ifdef TEST

#include <stdlib.h>
#include <unistd.h>

#include <cstdio>

int main(int argc, char *argv[]) {
  poly_interpolator<float> poly;
  double a = -1, b = 1, c = 0.1;
  int opt;
  while ((opt = getopt(argc, argv, "a:b:c:d:")) != -1)
    switch (opt) {
    case 'a':
      a = atof(optarg);
      break;
    case 'b':
      b = atof(optarg);
      break;
    case 'c':
      c = atof(optarg);
      break;
    case 'd':
      poly.set_abscissa_thres(atof(optarg));
    }
  double x, y;
  while (optind < argc && sscanf(argv[optind++], " %lf,%lf", &x, &y) == 2)
    poly.add(x, y);
  poly.interpolate();
  for (x = a; x < b; x += c)
    printf("%lf,%lf\n", x, poly(x));
  return EXIT_SUCCESS;
}

// g++ -o polyinterp -O2 -DTEST -x c++ polyinterp.h

#endif

template <>
enum slatec_polint_status polint<double>(size_t n, const double x[],
                                         const double y[], double c[]) {
  return slatec_polint(n, x, y, c);
}

template <>
enum slatec_polyvl_status polyvl<double>(double xx, double *yy, size_t n,
                                         const double x[], const double c[]) {
  return slatec_polyvl(xx, yy, n, x, c);
}

template <>
enum slatec_polint_status polint<float>(size_t n, const float x[],
                                        const float y[], float c[]) {
  return slatec_polintf(n, x, y, c);
}

template <>
enum slatec_polyvl_status polyvl<float>(float xx, float *yy, size_t n,
                                        const float x[], const float c[]) {
  return slatec_polyvlf(xx, yy, n, x, c);
}

#endif // __cplusplus
