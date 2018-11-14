#ifndef FINITE_DIFF_H_
#define FINITE_DIFF_H_

namespace matlib {
  
  template <typename F, typename D>
  inline auto finite_difference_forward(F f, D d, D dt = 0.01){
    return (f(d + dt) - f(d)) / dt;
  }
  
  template <typename F, typename D>
  inline auto finite_difference_backward(F f, D d, D dt = 0.01){
    return (f(d) - f(d - dt)) / dt;
  }
  
  template <typename F, typename D>
  inline auto finite_difference(F f, D d, D dt = 0.01){
    return (f(d + dt) - f(d - dt)) / (2.0 * dt);
  }
  
} // namespace matlib

#endif
