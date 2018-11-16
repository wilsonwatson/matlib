#ifndef BFGS_H_
#define BFGS_H_

#include "constants.h"

#include <math.h>
#include <tuple>
#include <vector>
#include <stdio.h>
#include <optional>
#include <fstream>
#include <utility>
#include <functional>


#define CLAMP(x, min, max) (x > max ? max : (x < min ? min : x))

template <typename F1, typename F2, typename OptType = std::ofstream*>
inline std::tuple<bool, int, int> bfgs(std::vector<double>& x, F1 func, F2 deriv,
                                       std::optional<double> tol = std::nullopt, std::optional<int> maxItrOpt = std::nullopt, OptType file = nullptr) {
  double tolerance = tol.value_or(1.0e-6);
  int maxItr = maxItrOpt.value_or(1000);
  auto norm = [](std::vector<double> n) {
    double sum = 0;
    for (int i = 0; i < n.size(); i++) {
      sum += n[i] * n[i];
    }
    return sqrt(sum);
  };

  int funcNum = 0;
  auto step = [&](std::vector<double> xk, std::vector<double> dk,
                  int maxItr = 10) {
    double mu = 0.001, kUp = 0.5, kLow = 0.1, alpha = 1.0, alphaMin, alphaMax;
    double fNew, fk = func(xk);
    std::vector<double> xNew(xk.size()), gk = deriv(xk);
    double gd = 0;
    for (int i = 0; i < gk.size(); i++) {
      gd += gk[i] * dk[i];
    }
    for (int i = 0; i < maxItr; i++) {
      for (int j = 0; j < xk.size(); j++) {
        xNew[j] = xk[j] + alpha * dk[j];
      }
      fNew = func(xNew);
      funcNum++;
      if (fNew < fk + mu * alpha * gd) {
        return std::make_tuple(true, alpha);
      } else {
        alphaMin = kLow * alpha;
        alphaMax = kUp * alpha;

        alpha = -0.5 * alpha * alpha * gd / (fNew - fk - alpha * gd);
        alpha = CLAMP(alpha, alphaMin, alphaMax);
      }
    }
    if(fNew >= fk){
      return std::make_tuple(false, alpha);
    }
    return std::make_tuple(true, alpha);
  };
  int k = 0, cnt = 0, N = x.size();

  double ys, yHy, alpha;
  std::vector<double> d(N), s(N), y(N), v(N), Hy(N), gPrev(N);

  std::vector<double> H(N * N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      H[i * N + j] = (i == j) ? 1 : 0;
    }
  }

  double fx = func(x);
  funcNum++;

  std::vector<double> gnorm(maxItr), g = deriv(x);
  gnorm[k++] = norm(g);

  while ((gnorm[k - 1] > tolerance) && (k < maxItr)) {
    std::vector<double> d(g.size());
    for (int i = 0; i < N; i++) {
      d[i] = 0;
      for (int j = 0; j < N; j++) {
        d[i] -= H[i * N + j] * g[j];
      }
    }
    auto [success, dalpha] = step(x, d);
    alpha = dalpha;
    if(file){
      (*file) << (k) << " " << fx << " " << alpha << " ";
      for(int i = 0; i < N; i++){
        (*file) << (x[i]) << " " << g[i] << " ";
      }
      (*file) << "\n";
    }
    //printf("%f\n", alpha);
    if(!success){
      std::vector<double> HmE(H.size());
      for(int i = 0; i < H.size(); i++){
        int r = (i / N);
        double offset = ((i - r * N) == r) ? -1 : 0;
        HmE[i] = H[i] - offset;
      }
      if( norm(HmE) < EPS )
        break;
      else
      {
        for(int i = 0; i < N; i++){
          for(int j = 0; j < N; j++){
            H[i * N + j] = (j == i) ? 1 : 0;
          }
        }
        cnt++;
        if( cnt == maxItr )
          break;
      }
    } else {
      /*
       s = alpha * d;*
       x += s;
       fx = func(x);
       this->funcNum++;
       gPrev = g;
       g = func.grad(x);
       y = g - gPrev;
       
       Hy = H * y;
       ys = dotProd( y, s );
       yHy = dotProd( y, Hy );
       if( (ys < EPS) || (yHy < EPS) )
         H = eye( N, Dtype(1.0) );
       else
       {
       v = sqrt(yHy) * ( s/ys - Hy/yHy );
       H = H + multTr(s,s)/ys - multTr(Hy,Hy)/yHy + multTr(v,v);
       }
       gnorm[k++] = norm(g);
      */
      for(int i = 0; i < N; i++){
        s[i] = alpha * d[i];
        x[i] += s[i];
      }
      fx = func(x);
      funcNum++;
      gPrev = g;
      g = deriv(x);
      for(int i = 0; i < N; i++){
        y[i] = g[i] - gPrev[i];
      }
      ys = 0;
      yHy = 0;
      for (int i = 0; i < N; i++) {
        Hy[i] = 0;
        for (int j = 0; j < N; j++) {
          Hy[i] += H[i * N + j] * y[j];
        }
        ys += y[i] * s[i];
        yHy += Hy[i] * y[i];
      }
      if(ys < EPS || yHy < EPS){
        for(int i = 0; i < N; i++){
          for(int j = 0; j < N; j++){
            H[i * N + j] = (j == i) ? 1 : 0;
          }
        }
      } else {
        double sqrtyHy = sqrt(yHy);
        for(int i = 0; i < N; i++){
          v[i] = sqrtyHy * (s[i] / ys - Hy[i] / yHy);
        }
        std::vector<double> Hnew(N * N);
        for(int i = 0; i < N; i++){
          for(int j = 0; j < N; j++){
            Hnew[i * N + j] = H[i * N + j] + (s[i] * s[j]) / ys - (Hy[i] * Hy[j]) / yHy + (v[i] * v[j]);
          }
        }
        for(int i = 0; i < N*N; i++){
          H[i] = Hnew[i];
        }
      }
      gnorm[k++] = norm(g);
    }
  }

  bool success = gnorm[k - 1] <= tolerance;
  return std::make_tuple(success, k, funcNum);
}

template <typename F1, typename F2, typename OptType = std::ofstream*>
inline std::tuple<bool, int, int> bfgs(std::vector<double>& x, F1 func, F2 deriv,
                                       OptType file, std::optional<double> tol = std::nullopt, std::optional<int> maxItrOpt = std::nullopt){
  return bfgs(x, func, deriv, tol, maxItrOpt, file);
                                       }

#undef CLAMP

#endif
