#ifndef FASTMATH_H_
#define FASTMATH_H_

#include "tpool.h"

#include <functional>
#include <optional>
#include <vector>

namespace matlib {
  
int next_pow_2(int v){
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return v;
}

template <unsigned I = 0>
inline void
matmult_strassen(std::vector<double> left, int left_rows, int left_cols,
                 int left_offset, std::vector<double> right, int right_rows,
                 int right_cols, int right_offset, std::vector<double>& out, int out_rows, int out_cols, int out_offset,
                 thread_pool<I>& tp) {
  int left_major = left_rows > left_cols ? left_rows : left_cols;
  int right_major = right_rows > right_cols ? right_rows : right_cols;
  int out_major = out_rows > out_cols ? out_rows : out_cols;
  int major = next_pow_2(left_major > right_major ? left_major : right_major);
  printf("%i, %i, %i\n", left_major, right_major, major);
  if(major == 1){ // TODO implement strassens algorithm with threading
    out[out_offset] = left[left_offset] * right[right_offset];
    return;
  }
  
  std::vector<std::future<double>> parts;
}

}; // namespace matlib

#endif
