#include "matlib/matlib.h"

int main() {
  std::vector<double> x(10);
  for (int i = 0; i < 10; i++)
    x[i] = 0;
  auto rosen = [](std::vector<double> x) {
    double fx = 0.0;
    for (int i = 0; i < 10; i += 2) {
      double t1 = 1.0 - x[i];
      double t2 = 10 * (x[i + 1] - x[i] * x[i]);
      fx += t1 * t1 + t2 * t2;
    }
    return fx;
  };
  auto rosengrad = [](std::vector<double> x) {
    std::vector<double> grad(10);
    for (int i = 0; i < 10; i += 2) {
      double t1 = 1.0 - x[i];
      double t2 = 10 * (x[i + 1] - x[i] * x[i]);
      grad[i + 1] = 20 * t2;
      grad[i] = -2.0 * (x[i] * grad[i + 1] + t1);
    }
    return grad;
  };
  std::ofstream out;
  out.open("test.dat");
  auto [success, itr, calls] = bfgs(x, rosen, rosengrad, &out);
  printf("%s! x ended as <", success ? "Success": "Failure");
  for(int i = 0; i < 10; i++){
    printf("%s%f", (i == 0) ? "" : ", ", x[i]);
  }
  printf("> with %i iterations and %i function calls\n", itr, calls);
  std::ofstream script;
  script.open("script.plot");
  script << "set xlabel 'iteration'\nset ylabel 'value'\n";
  //script << "plot ";
  for(int i = 0; i < 10; i++){
    script << "plot 'test.dat' using 1:" << (i * 2) + 2 << " with lines title 'x_" << i << "'\npause -1 'press a key to continue'\n";
    //script << (i == 0 ? "" : ", ") << "'test.dat' using 1:" << (i * 2) + 2 << " with lines title 'x_" << i << "'";
  }
  script << "\n";
  return 0;
}
