#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include "payoff.hpp"

namespace dauphine

{
class boundary
{
    public:
        boundary();
        double get_upper_b() const;
        double get_lower_b() const;


    protected:
        double upper_b;
        double lower_b;        
};

class dirichlet : public boundary
{
    public:
        dirichlet(dauphine::payoff& payoff, double s_min, double s_max);
};

}


#endif