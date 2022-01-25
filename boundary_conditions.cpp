#include "boundary_conditions.hpp"


namespace dauphine

{

boundary::boundary(){}

dirichlet::dirichlet(dauphine::payoff& payoff, double s_min, double s_max)
{
    upper_b = payoff(s_max);
    lower_b = payoff(s_min);
}

double boundary::get_upper_b() const
{
    return upper_b;
}

double boundary::get_lower_b() const
{
    return lower_b;
}

}