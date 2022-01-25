#include "volatility.hpp"

namespace dauphine
{

volatility::volatility(){};

volatility::volatility(std::vector<double> v){
    p_vol.resize(v.size());
    p_vol = v;
};

vol_BS::vol_BS(double sigma)
{
    p_vol.resize(1);
    p_vol[0] = sigma;
}

vol_gen::vol_gen(double sigma, dauphine::mesh msh,  double alpha)
{
    const int n_dt = msh.get_ndt();
    p_vol.resize(n_dt);
    p_vol[0] = sigma;
    for(int i=1; i<n_dt ; i++){
            p_vol[i] = p_vol[i-1] + alpha;
        };
}

vol_heston::vol_heston(const double& sigma,const double& kappa, const double& theta , const dauphine::mesh& msh)
{
    const int n_dt = msh.get_ndt();
    const double dt = msh.get_dt();
    p_vol.resize(n_dt);
    p_vol[0] = sigma;
    for(int i=1; i<n_dt ; i++){
        p_vol[i] = p_vol[i-1] + kappa * (theta - p_vol[i-1]) * dt;
    };


}

std::vector<double> volatility::get_vol() const
{
    return p_vol;
}
}