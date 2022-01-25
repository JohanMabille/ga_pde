#ifndef PAYOFF_HPP
#define PAYOFF_HPP

#include <vector>

namespace dauphine

{

class payoff
{
    public :
        payoff();
        virtual double operator() (const double& spot) const = 0;
        virtual std::vector<double> operator() (const std::vector<double>& spot) const = 0;
};

class vanilla_poff : public payoff
{
    private : 
        double p_strike;
        bool p_is_call;

    public:
        vanilla_poff(const double& strike, const bool& is_call);
        double get_strike() const;
        bool get_type() const;
        virtual double operator() (const double& spot) const;
        virtual std::vector<double> operator() (const std::vector<double>& spot) const;

};


};



#endif