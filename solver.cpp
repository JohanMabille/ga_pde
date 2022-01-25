#include "solver.hpp"
#include <algorithm>

namespace dauphine
{

    solver::solver(payoff& poff, mesh& msh, boundary& bd, double theta)
    : m_poff(poff), m_msh(msh), m_bd(bd), m_theta(theta)
    {}; 

    void solver::transform_coeff(std::vector<double>& a,
                        std::vector<double>& b,
                        std::vector<double>& c,
                        std::vector<double>& d)
    {
        const double dt = m_msh.get_dt();
        const double dx = m_msh.get_dx();
        const double nu1 = dt / pow(dx,2);
        const double nu2 = dt / dx;
    
    
        std::vector<double> A = a;
        std::vector<double> B = b;
        std::vector<double> C = c;
        std::vector<double> D = d;

        std::transform(A.begin(),A.end(),B.begin(),a.begin(),[nu1,nu2](double xA,double xB)-> double 
        {return xA*nu1 - 0.5 * xB * nu2;});

        std::transform(C.begin(),C.end(),A.begin(),b.begin(),[nu1,nu2,dt](double xC,double xA)-> double 
        {return xC * dt - 2 * xA * nu1;});

        std::transform(A.begin(),A.end(),B.begin(),c.begin(),[nu1,nu2](double xA,double xB)-> double 
        {return xA * nu1  + 0.5 * xB * nu2;});

        std::transform(D.begin(),D.end(),d.begin(),[dt](double xD)-> double 
        {return xD * dt;});

    };

    void solver::transform_coeff(double& a, double& b, double& c, double& d)
    {
        const double nu1 = m_msh.get_dt() / pow(m_msh.get_dx(),2);
        const double nu2 = m_msh.get_dt() / m_msh.get_dx();

        double A = a;
        double B = b;
        double C = c;
        double D = d;

        a = A*nu1 - 0.5 * B * nu2;
        b = C * m_msh.get_dt() - 2 * A * nu1;
        c = A * nu1  + 0.5 * B * nu2;
        d = D * m_msh.get_dt();
    };

    void solver::init_matrice_1(matrix& m_trans,const int& dim, double a, double b, double c, double d)
    {
        m_trans(0,0) = m_theta*b+1;
        m_trans(0,1) = m_theta*c;

        for(int i=1; i<dim-1; i++){
            m_trans(i,i-1) = m_theta*a;
            m_trans(i,i) = m_theta*b+1;
            m_trans(i,i+1) = m_theta*c;
        };

        m_trans(dim-1,dim-2) = m_theta*a;
        m_trans(dim-1,dim-1) = m_theta*b+1;

        
    }

    void solver::init_matrice_1(matrix& m_trans, const int& dim, const int& j,
                        std::vector<double>& a,
                        std::vector<double>& b,
                        std::vector<double>& c,
                        std::vector<double>& d)
    {
        m_trans(0,0) = m_theta*b[j]+1;
        m_trans(0,1) = m_theta*c[j];

        for(int i=1; i<dim-1; i++){
            m_trans(i,i-1) = m_theta*a[j];
            m_trans(i,i) = m_theta*b[j]+1;
            m_trans(i,i+1) = m_theta*c[j];
        };

        m_trans(dim-1,dim-2) = m_theta*a[j];
        m_trans(dim-1,dim-1) = m_theta*b[j]+1;
    }

    void solver::init_matrice_2(matrix& m_trans,const int& dim, double a, double b, double c, double d)
    {
        m_trans(0,0) = 1-(1-m_theta)*b;
        m_trans(0,1) = -(1-m_theta)*c;

        for(int i=1; i<dim-1; i++){
            m_trans(i,i-1) = -(1-m_theta)*a;
            m_trans(i,i) = 1-(1-m_theta)*b;
            m_trans(i,i+1) = -(1-m_theta)*c;
        };

        m_trans(dim-1,dim-2) = -(1-m_theta)*a;
        m_trans(dim-1,dim-1) = 1-(1-m_theta)*b;        
    }

    void solver::init_matrice_2(matrix& m_trans,const int& dim, const int& j,
                                std::vector<double>& a,
                                std::vector<double>& b,
                                std::vector<double>& c,
                                std::vector<double>& d)
    {
        m_trans(0,0) = 1-(1-m_theta)*b[j+1];
        m_trans(0,1) = -(1-m_theta)*c[j+1];

        for(int i=1; i<dim-1; i++){
            m_trans(i,i-1) = -(1-m_theta)*a[j+1];
            m_trans(i,i) = 1-(1-m_theta)*b[j+1];
            m_trans(i,i+1) = -(1-m_theta)*c[j+1];
        };

        m_trans(dim-1,dim-2) = -(1-m_theta)*a[j+1];
        m_trans(dim-1,dim-1) = 1-(1-m_theta)*b[j+1];        
    }

    void solver::fill_matrix(matrix& mesh_matrix, int t, std::vector<double> vect)
    {
        for (int x=0; x<mesh_matrix.nb_rows();x++){
            mesh_matrix(x,t) = vect[mesh_matrix.nb_rows()-x-1];
        }
    }


    std::vector<double> solver::solve_tridiag(const matrix& mat, std::vector<double>& d){

        if(mat.nb_cols()!=mat.nb_rows()){
            throw std::domain_error("Not square matrix" );
        }

        if(mat.nb_cols()!=d.size()){
            throw std::domain_error("b of wrong size when solving Ax=b" );
        }
        
        std::size_t n = mat.nb_cols();
        std::vector<double> a(n-1);
        std::vector<double> b(n);
        std::vector<double> c(n-1);

        for(std::size_t i=0; i<n; i++){
            for( std::size_t j=0; j<n;j++){
                if(j==i+1){
                    c[i]=mat(i,j);
                }
                else if(j==i-1){
                    a[j]=mat(i,j);
                }
                else if (j==i){
                    b[i]=mat(i,j);
                }
            }
        }

        std::vector<double> x(n);

        for(int i=1; i<=n-1; i++){
            
            double w = double(a[i-1])/double(b[i-1]);
            b[i] -= w*c[i-1];
            d[i] -= w*d[i-1];
        }
        x[n-1] = d[n-1]/b[n-1];
        for(int i=n-1; i >= 1; i--)
            x[i-1] = (d[i-1]- c[i-1]*x[i])/b[i-1];

		return x;
        
    }

    void solver::solve_mesh(matrix& m_trans_1,matrix& m_trans_2,
                            double& a, double& b, double& c, double& d,
                            const int& ndx,
                            std::vector<double>& final_vect,std::vector<double>& final_poff,std::vector<double>& vect,
                            matrix& mesh_matrix)
    {
        this->init_matrice_1(m_trans_1,ndx-2,a,b,c,d); 
        this-> init_matrice_2(m_trans_2,ndx-2,a,b,c,d);

             
        vect[0] = -m_bd.get_lower_b() * a + d;
        vect[ndx-3] = - m_bd.get_upper_b()* a + d;

        for (int i=(m_msh.get_ndt()-2); i>=0;i--){

            final_vect = m_trans_2 * final_vect;
            final_vect = final_vect + vect;
            final_vect = this->solve_tridiag(m_trans_1,final_vect);

            std::copy(final_vect.begin(),final_vect.end(),final_poff.begin()+1);

            this->fill_matrix(mesh_matrix,i,final_poff);
        };
    }
    
    void solver::solve_mesh(matrix& m_trans_1,matrix& m_trans_2,
                            std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d,
                            const int& ndx,
                            std::vector<double>& final_vect,std::vector<double>& final_poff,std::vector<double>& vect,
                            matrix& mesh_matrix)
    {

        for (int i=(m_msh.get_ndt()-2); i>=0;i--){

            this->init_matrice_1(m_trans_1,ndx-2,i,a,b,c,d); 
            this->init_matrice_2(m_trans_2,ndx-2,i,a,b,c,d);

            
            vect[0] = -m_theta*(m_bd.get_lower_b() * a[i] + d[i]) 
                        - (1-m_theta)*(m_bd.get_lower_b() * a[i+1] + d[i+1]);

            vect[ndx-3] = -m_theta*(m_bd.get_upper_b()* a[i] + d[i])
                            - (1-m_theta)*(m_bd.get_upper_b()* a[i+1] + d[i+1]);

            final_vect = m_trans_2 * final_vect;
            final_vect = final_vect + vect;
            final_vect = this->solve_tridiag(m_trans_1,final_vect);

            std::copy(final_vect.begin(),final_vect.end(),final_poff.begin()+1);

            this->fill_matrix(mesh_matrix,i,final_poff);
        };
    }



    matrix solver::price(const pde_european& pde)
    {
        const int size = pde.get_coeff_a().size();

        if (size==1){
            double a = pde.get_coeff_a()[0];
            double b = pde.get_coeff_b()[0];
            double c = pde.get_coeff_c()[0];
            double d = pde.get_coeff_d()[0];
            return this->compute_price(a,b,c,d);
            }
        else{
            std::vector<double> a = pde.get_coeff_a();
            std::vector<double> b = pde.get_coeff_b();
            std::vector<double> c = pde.get_coeff_c();
            std::vector<double> d = pde.get_coeff_d();
            return this->compute_price(a,b,c,d);
            }
    }

    matrix solver::compute_price(double& a, double& b, double& c, double& d)
    {       
            this->transform_coeff(a,b,c,d);

            const int ndx = m_msh.get_ndx();
            const int ndt = m_msh.get_ndt();

            matrix m_trans_1(ndx-2,ndx-2);
            matrix m_trans_2(ndx-2,ndx-2);

            std::vector<double> vect(ndx-2); 

            matrix mesh_matrix(ndx,ndt);
        
            std::vector<double> final_poff = m_poff(m_msh.get_xaxis());
            std::vector<double> final_vect(ndx-2);

            this->fill_matrix(mesh_matrix,mesh_matrix.nb_cols()-1,final_poff); //Compute complete final cond

            std::copy(final_poff.begin()+1,final_poff.end()-1,final_vect.begin()); // extract only the core part in final_vect


            this->solve_mesh(m_trans_1,m_trans_2,a,b,c,d,ndx,final_vect,final_poff,vect,mesh_matrix);

            return mesh_matrix;
    
    }

    matrix solver::compute_price(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d)
    {       
            this->transform_coeff(a,b,c,d);

            const int ndx = m_msh.get_ndx();
            const int ndt = m_msh.get_ndt();

            matrix m_trans_1(ndx-2,ndx-2);
            matrix m_trans_2(ndx-2,ndx-2);

            std::vector<double> vect(ndx-2);   

            matrix mesh_matrix(ndx,ndt);
        
            std::vector<double> final_poff = m_poff(m_msh.get_xaxis());
            std::vector<double> final_vect(ndx-2);

            this->fill_matrix(mesh_matrix,mesh_matrix.nb_cols()-1,final_poff); //Compute complete final cond

            std::copy(final_poff.begin()+1,final_poff.end()-1,final_vect.begin()); // extract only the core part in final_vect


            this->solve_mesh(m_trans_1,m_trans_2,a,b,c,d,ndx,final_vect,final_poff,vect,mesh_matrix);

            return mesh_matrix;

            
    
    }


}
