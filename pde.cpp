#include "pde.hpp"

namespace dauphine 
{
    pde::pde(){};

    pde_european::pde_european(const volatility& vol, const rate& rate){
        const int size_vol = vol.get_vol().size();
        const int size_rate = rate.get_rates().size();
        const int size = std::max(size_vol,size_rate);

        p_a.resize(size);
        p_b.resize(size);
        p_c.resize(size);
        p_d.resize(size);

        for (int i=0; i<size; i++)
        {
            int k = std::min(i,size_vol-1);
            int j = std::min(i,size_rate-1);
            p_a[i] = - 0.5 * pow(vol.get_vol()[k],2);
            p_b[i] = 0.5 * pow(vol.get_vol()[k],2) - rate.get_rates()[j];
            p_c[i] = rate.get_rates()[j];
            p_d[i]=0; 
        }
    };

    std::vector<double> pde_european::get_coeff_a() const{
        return p_a;
    }

    std::vector<double> pde_european::get_coeff_b() const{
        return p_b;
    }

    std::vector<double> pde_european::get_coeff_c() const{
        return p_c;
    }

    std::vector<double> pde_european::get_coeff_d() const{
        return p_d;
    }


}