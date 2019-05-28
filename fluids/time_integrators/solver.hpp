#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "../utils/point_functions.hpp"
#include <memory>

class Options;
class Derivatives;
class Convection;
class Dissipation;
class Boundary;

class Solver {
public:
    Solver(Options& opt_in);
    void run();

private:
    Options& opt;
    PointFunctions pf;
    template <typename Grid>
    void setup_and_run();
};

#endif /* SOLVER_HPP */
