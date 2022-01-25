#include "payoff.hpp"
#include <algorithm> //for std::max function
#include <math.h>

namespace dauphine

{
payoff::payoff(){} 

vanilla_poff::vanilla_poff(const double& strike, const bool& is_call)
: p_strike(strike), p_is_call(is_call)
{
}

double vanilla_poff::operator() (const double& spot) const
{
    
    return std::max(p_is_call ? exp(spot) - p_strike : p_strike - exp(spot), 0.);
}

std::vector<double> vanilla_poff::operator() (const std::vector<double>& spot) const
{
    std::vector<double> payout(spot.size());
    for( int i=0; i<spot.size();i++){
        payout[i]=std::max(p_is_call ? exp(spot[i]) - p_strike : p_strike - exp(spot[i]), 0.);
    }
    return payout;
}

double vanilla_poff::get_strike() const{
    return p_strike;
}

bool vanilla_poff::get_type() const{
    return p_is_call;
}

}