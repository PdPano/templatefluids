#include "input_output/options.hpp"
#include "time_integrators/solver.hpp"
#include <fstream>
#include <iostream>

int main() {
  std::ifstream config_file("./solver.cfg");
  Options opt(config_file);

  Solver solver(opt);
  solver.run();

  return 0;
}
