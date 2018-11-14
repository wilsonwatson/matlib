#ifndef FRESNEL_H_
#define FRESNEL_H_

#include "gauss_quad.h"

namespace matlib {

template <unsigned V = 16>
double fresnelC(double x, double sharpness = 1.0) {
  return gauss_legendre_integrate<V>([&](double x) { return cos(M_PI * sharpness * x * x / 2.0); }, 0,
                            x);
}

template <unsigned V = 16>
double fresnelS(double x, double sharpness = 1.0) {
  return gauss_legendre_integrate<V>([&](double x) { return sin(M_PI * sharpness * x * x / 2.0); }, 0,
                            x);
}

} // namespace fresnel

#endif
