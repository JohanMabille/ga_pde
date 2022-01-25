#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "payoff.hpp"
#include "mesh.hpp"
#include "volatility.hpp"
#include "boundary_conditions.hpp"
#include "rate.hpp"
#include "matrix.hpp"
#include "pde.hpp"
#include <vector>
#include <ostream>
#include <iostream>
#include <math.h>


namespace dauphine
{

    class solver {
        public : 

        friend std::ostream& operator<<(std::ostream& os, solver const & s) {
            
            os << "Spot is : " << s.m_msh.get_spot() << std::endl;

            os << "Smax  : " << exp(s.m_msh.get_xmax()) << std::endl;
            os << "Smin  : " << exp(s.m_msh.get_xmin()) << std::endl;
            os << std::endl;

            os << "Upper Bound : " << s.m_bd.get_upper_b() << std::endl;
            os << "Lower Bound : " << s.m_bd.get_lower_b() << std::endl << std::endl;
            os << std::endl;

            os << "ndt : " << s.m_msh.get_ndt() << std::endl;
            os << "ndx : " << s.m_msh.get_ndx() << std::endl;
            os << std::endl;

    
        return os ;
         
        }

        solver(payoff& poff, mesh& msh, boundary& bd,double theta);

        matrix price(const pde_european& pde);
        matrix compute_price(double& a, double& b, double& c, double& d);
        matrix compute_price(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d);

        void solve_mesh(matrix& m_trans_1,matrix& m_trans_2,
                            double& a, double& b, double& c, double& d,
                            const int& ndx,
                            std::vector<double>& final_vect,std::vector<double>& final_poff,std::vector<double>& vect,
                            matrix& mesh_matrix);
        void solve_mesh(matrix& m_trans_1,matrix& m_trans_2,
                            std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d,
                            const int& ndx,
                            std::vector<double>& final_vect,std::vector<double>& final_poff,std::vector<double>& vect,
                            matrix& mesh_matrix);

        double get_coeff_a(const pde_european& pde);
        double get_coeff_b(const pde_european& pde);
        double get_coeff_c(const pde_european& pde);
        double get_coeff_d(const pde_european& pde);

        void transform_coeff(double& a, double& b, double& c, double& d);
        void transform_coeff(std::vector<double>& a,
                        std::vector<double>& b,
                        std::vector<double>& c,
                        std::vector<double>& d);

        void init_matrice_1(matrix& m_trans, const int& dim, double a, double b, double c, double d);
        void init_matrice_1(matrix& m_trans, const int& dim, const int& i,
                                std::vector<double>& a,
                                std::vector<double>& b,
                                std::vector<double>& c,
                                std::vector<double>& d);

        void init_matrice_2(matrix& m_trans, const int& dim, double a, double b, double c, double d);
        void init_matrice_2(matrix& m_trans,const int& dim, const int& i,
                                std::vector<double>& a,
                                std::vector<double>& b,
                                std::vector<double>& c,
                                std::vector<double>& d);

        void fill_matrix(matrix& mesh_matrix, int t, std::vector<double> vect);

        std::vector<double> solve_tridiag(const matrix& mat, std::vector<double>& d);

        private :
            mesh& m_msh;
            payoff& m_poff;
            boundary& m_bd;
            double m_theta;

    };

}

#endif