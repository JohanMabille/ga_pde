#include "greeks.hpp"
#include "matrix.hpp"
#include <vector>


namespace dauphine {
    double delta(const matrix& result, const std::vector<double>& x_axis, const int& index){
        return (result(index+1,0)- result(index-1,0))/(exp(x_axis[index-1])-exp(x_axis[index+1]));
    }

    double gamma(const matrix& result, const std::vector<double>& x_axis){
        int index = x_axis.size()/2;
        return (delta(result,x_axis,index)-delta(result,x_axis,index-1))/(exp(x_axis[index-1])-exp(x_axis[index]));
    }

    double vega(const matrix& result, solver& solv, const volatility& vol, const rate& r, const int& ndx){
        std::vector<double> new_vol = vol.get_vol() + 0.01;
        volatility new_v(new_vol);
        pde_european pd(new_v,r);
        matrix new_price = solv.price(pd);
        return (new_price(ndx/2,0) - result(ndx/2,0))*100 ; 
    }

    double rho(const matrix& result, solver& solv, const volatility& vol, const rate& r, const int& ndx){
        std::vector<double> new_rate = r.get_rates() + 0.01;
        rate new_r(new_rate);
        pde_european pd(vol,new_r);
        matrix new_price = solv.price(pd);
        return (new_price(ndx/2,0) - result(ndx/2,0))*100 ; 
    }

    double theta(const matrix& result, const double& dt, const double& index){
        return (result(index,1)-result(index,0))/dt;
    }

}