#ifndef GAUSS_QUAD_H_
#define GAUSS_QUAD_H_

#include <math.h>
#include <tuple>

namespace matlib {

namespace internal {

template <class T, class Tuple, size_t... Is>
constexpr auto tuple_plus_1(Tuple &&tuple, std::index_sequence<Is...>, T t) {
  return std::make_tuple(std::get<Is>(std::forward<Tuple>(tuple))..., t);
}

template <class T, class Tuple>
constexpr auto tuple_plus_1(Tuple &&tuple, T t) {
  return tuple_plus_1<T>(
      std::forward<Tuple>(tuple),
      std::make_index_sequence<std::tuple_size<std::decay_t<Tuple>>::value>{},
      t);
}

template <unsigned Order>
inline static std::pair<double, double>
compute_legendre_polynomial_root(int index) {
  static auto eval = [](double x) {
    double vsub1 = x;
    double vsub2 = 1;
    double f = 1 / (x * x - 1);
    double _v, _d;

    for (int i = 2; i <= Order; ++i) {
      _v = ((2 * i - 1) * x * vsub1 - (i - 1) * vsub2) / i;
      _d = i * f * (x * _v - vsub1);

      vsub2 = vsub1;
      vsub1 = _v;
    }
    return std::make_pair(_v, _d);
  };
  double dr = 1;
  double x = cos(M_PI * ((index + 1) - 0.25) / (Order + 0.5));
  auto ev = eval(x);
  do {
    dr = ev.first / ev.second;
    x -= dr;
    ev = eval(x);
  } while (fabs(dr) > 2e-16);
  return std::make_pair(x, 2.0 / ((1.0 - x * x) * ev.second * ev.second));
}

template <unsigned Order, unsigned Index> struct legendre_polynomial_root {
  using type = decltype(tuple_plus_1(
      std::declval<typename legendre_polynomial_root<Order, Index - 1>::type>(),
      std::pair<double, double>{}));
  inline static type val() {
    return tuple_plus_1(legendre_polynomial_root<Order, Index - 1>::val(),
                        compute_legendre_polynomial_root<Order>(Index));
  }
};

template <unsigned Order> struct legendre_polynomial_root<Order, 0> {
  using type = std::tuple<std::pair<double, double>>;
  inline static type val() {
    return std::make_tuple(compute_legendre_polynomial_root<Order>(0));
  }
};

template <unsigned Order> struct legendre_polynomial {
  inline static typename legendre_polynomial_root<Order, Order - 1>::type
  roots() {
    return legendre_polynomial_root<Order, Order - 1>::val();
  }
};

template <unsigned I> struct accum {
  template <typename Tup, typename F>
  inline static double eval(const Tup &t, F f) {
    return f(I, std::get<I>(t)) + accum<I - 1>::eval(t, f);
  }
};

template <> struct accum<0> {
  template <typename Tup, typename F>
  inline static double eval(const Tup &t, F f) {
    return f(0, std::get<0>(t));
  }
};

template <unsigned I, typename Tup, typename F>
inline double accumulate(const Tup &t, F f) {
  return accum<I - 1>::eval(t, f);
}

} // namespace internal

template <unsigned numPoints = 16, typename F = void>
inline double gauss_legendre_integrate(F f, double a, double b) {
  static auto roots = internal::legendre_polynomial<numPoints>::roots();
  return ((b - a) / 2.0) *
         internal::accumulate<numPoints>(
             roots, [&](int, std::pair<double, double> k) {
               return f((b - a) / 2.0 * k.first + (a + b) / 2.0) * k.second;
             });
}

} // namespace matlib

#endif
