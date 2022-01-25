#include "closed_form.hpp"
#include "payoff.hpp"
#include "boundary_conditions.hpp"
#include "volatility.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "rate.hpp"
#include "matrix.hpp"
#include "pde.hpp"
#include "greeks.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h> 


// Guidelines:
//
// 1] Equation
// Even if the volatility and the rate are constant in the BS model,
// we can prove that, under certain assumptions, we can use the same
// PDE equation with volatility surface and rate curves. We could also
// want to take into account the repo (even if it could theoretically
// be part of the r factor). Therefore, it is more generic to solve
// the following generic equation;
// df/dt = a(x,t)d2f/dx2 + b(x,t)df/dx + c(x, t)f + d(x, t).
// The pricer should ask the coefficients a, b, c and d to an
// abstract class that will be inherited by classes implementing
// different models.
// 
// 2] Payoff
// The pricer should be able to price exotic options, such as
// barriers or asian options. To do so, the pricer must be able
// to call an external function between each step. Define an API
// that allows to register a function or an abstract class modeling
// a payoff.

 /* print vector */
    std::ostream &operator<<(std::ostream &out, const std::vector<double> &input)
    {
        for(std::size_t i = 0; i < input.size(); ++i)
            {
                    out << input[i] << " ";
            }
            out << std::endl;
            return out;
        return out;
    }


int main(int argc, const char * argv[])
{
    // Parameters of the option
    std::cout << std::endl;
    std::cout << "While waiting for the next release, the pricer is only able to price vanilla options." << std::endl;
    std::cout << "Do you want to price a call or a put ? Enter p for put or c for call : ";
    char type[1];
    std::cin >> type;
    bool is_call;
    if (strcmp(type,"c")==0){
        is_call = true;
    }
    else {
        is_call = false;
    }
    std::cout << std::endl;
    std::cout << "Now please give the different parameters of the option" << std::endl;
    std::cout << "Spot : ";
    double spot;
    std::cin >> spot;
    std::cout << "Strike : ";
    double strike;
    std::cin >> strike;
    std::cout << "Maturity (in years) : ";
    double maturity;
    std::cin >> maturity;
    std::cout << "Initial implied volatility (ex 0.10): ";
    double v0 ;
    std::cin >> v0;
    std::cout << "Initial interest rate (ex 0.05): ";
    double r0 ;
    std::cin >> r0;
    std::cout << std::endl;
    
    // Parameters of the mesh space 
    std::cout << "Now please specify the parameters of the mesh space." << std::endl;
    std::cout << "Number of points for the price discretisation (please enter odd number, ex 101 ): ";
    int ndx ;// enter odd number 
    std::cin >> ndx;
    std::cout << "Number of points for the time discretisation (ex 100) : ";
    int ndt ;// enter odd number 
    std::cin >> ndt;
    std::cout << std::endl;

    // Theta scheme parameter
    std::cout << "The price will be computed using a finite difference method. " << std::endl;
    std::cout << "Please give a theta scheme parameter between 0 and 1 : ";
    double theta;
    std::cin >> theta;
    std::cout << std::endl;

    // Parameters to implement an heston model 
    std::cout << "Do you want to implement an heston model ? y for yes, n for no : " ;
    char vol_type[1];
    std::cin >> vol_type;
    bool is_heston;
    double kappa = 0 ;
    double heston_theta = 0;
    if (strcmp(vol_type,"y")==0){
        is_heston = true;
        std::cout << "Kappa of the volatility diffusion (ex 2) : ";
        std::cin >> kappa;
        std::cout << "Theta of the volatility diffusion (ex 0.12) : ";
        std::cin >> heston_theta;
    }
    else {
        is_heston = false;
    }

    //Parameters for increasing interest rate :
    std::cout << "Do you want to use a constant interest rate ? y for yes, n for no : " ;
    char rate_type[1];
    std::cin >> rate_type;
    double alpha = 0 ;
    if (strcmp(rate_type,"n")==0){
        std::cout << "Enter a rate of increase at each time step (ex 0.0001): ";
        std::cin >> alpha;
    } 
    std::cout << std::endl;
    // Mesh space definition
    dauphine::mesh msh(spot, maturity,ndx,ndt,v0);
    dauphine::vanilla_poff poff(strike, is_call); 
    dauphine::dirichlet bound(poff,msh.get_xmin(),msh.get_xmax());

    // Constant parameters definition
    dauphine::vol_BS vol(v0);

    // Non-constant parameters definition
    dauphine::vol_heston vol_h(v0,kappa,heston_theta,msh);
    dauphine::rate rate(r0, msh, alpha);
    
    // Building finite difference method parameters
    dauphine::solver solv(poff, msh, bound, theta);
    std::cout << solv << std::endl;

    // Results
    double cf  =  dauphine::bs_price(spot,strike,v0,maturity,is_call) ;
    std::cout << "Closed form price : " << cf << " (with assumption r=0)" << std::endl;

    dauphine::matrix result(ndx,ndt);
    if(is_heston){
        dauphine::pde_european pde(vol_h,rate);
        result+=solv.price(pde);
    }
    else{
        dauphine::pde_european pde(vol,rate);
        result+=solv.price(pde);
    }
    std::cout << "Finite difference method price : " << result((ndx-1)/2,0) << std::endl;
    std::cout << "Difference : " << cf - result((ndx-1)/2,0) << std::endl << std::endl;
    
    //Greeks
    double delta = dauphine::delta(result,msh.get_xaxis(),ndx/2);
    std::cout << "Delta : " << delta << std::endl;

    double gamma = dauphine::gamma(result,msh.get_xaxis());
    std::cout << "Gamma : " << gamma << std::endl;

    double rho = dauphine::rho(result,solv,vol,rate, ndx);
    std::cout << "Rho : " <<  rho << std::endl;

    double vega = dauphine::vega(result,solv,vol,rate, ndx);
    std::cout << "Vega : " <<  vega << std::endl;

    double t = dauphine::theta(result,msh.get_dt(),ndx/2);
    std::cout << "Theta : " << t << std::endl;
    std::cout << std::endl;

    std::cout << "Do you want to see the full mesh matrix ? y for yes, n for no : ";
    char mat[1];
    std::cin >> mat;
    
    if (strcmp(mat,"y")==0){
        std::cout << std::endl;
        std::cout << result;
    } 

    return 0;
}
