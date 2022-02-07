#ifndef VOLATILITY_HPP
#define VOLATILITY_HPP

#include <iostream>
#include <vector>
#include "mesh.hpp"


namespace dauphine {

class volatility{

    public : 
        volatility();
        volatility(std::vector<double> v);
        // you should return a const reference for performances
        std::vector<double> get_vol() const;

    protected:
        std::vector<double> p_vol;
};

// You don't really use inheritance for redefining the behavior
// of the base class here. These classes simply encapsulate
// different method for filling a volatility vector, free
// functions (make_vol_BS, make_vol_gen, make_vol_heston),
// would have been enough
class vol_BS : public volatility
{
    public:
        vol_BS(double sigma);
};

class vol_gen : public volatility
{
    public:
        vol_gen(double sigma, dauphine::mesh msh,  double alpha);
    
};

class vol_heston : public volatility
{
    public:
        vol_heston(const double& sigma,const double& kappa, const double& theta , const dauphine::mesh& msh);
};

};

#endif
