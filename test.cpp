#include "matlib/matlib.h"

#include <iostream>
#include <sstream>
#include <unistd.h>

class koopman_nn {
public:
  inline koopman_nn(std::vector<unsigned> topology) : m_topology(topology) {
    unsigned num = 0;
    for (int i = 1; i < m_topology.size(); i++)
      num += m_topology[i] * m_topology[i - 1];
    m_weights.reserve(num);
  }
  inline ~koopman_nn() {}
  struct derivative {
    koopman_nn *koop;
    inline std::vector<double> operator()(std::vector<double> inputs){
      
      return inputs;
    }
  };
  struct cost {
    koopman_nn *koop;
    inline double operator()(std::vector<double> inputs){
      
      return 0;
    }
  };
  inline std::vector<double> operator()(std::vector<double> inputs = {}){
    
    return inputs;
  }
  
  
private:
  std::vector<unsigned> m_topology;
  std::vector<double> m_weights;
};

int test(int y){
  return y * y;
}

int main(int, char **) {
  thread_pool<8> tp;
  std::vector<double> A = {2, 3, 4}, B = {1, 2, 1, -1, 1, -1, 2};
  std::vector<double> C = {0, 0};
  matlib::matmult_strassen(A, 1, 3, 0, B, 3, 2, 0, C, 2, 1, 0, tp);
  std::cout << C[0] << std::endl;
  return 0;
}
