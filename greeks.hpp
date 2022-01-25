#ifndef GREEKS_HPP
#define GREEKS_HPP

#include "matrix.hpp"
#include "volatility.hpp"
#include "rate.hpp"
#include "solver.hpp"
#include <vector>
#include <math.h>


namespace dauphine
{
    double delta(const matrix& result, const std::vector<double>& x_axis, const int& index);
    double gamma(const matrix& result, const std::vector<double>& x_axis);
    double vega(const matrix& result, solver& solv, const volatility& vol, const rate& rate, const int& ndx);
    double rho(const matrix& result, solver& solv, const volatility& vol, const rate& rate, const int& ndx);
    double theta(const matrix& result, const double& dt, const double& index);

}

#endif