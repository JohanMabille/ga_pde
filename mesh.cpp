#include "mesh.hpp"
#include <math.h>

namespace dauphine
{
    mesh::mesh(double spot, double T, int n_dx, int n_dt, double vol)
    :m_spot(spot),m_T(T), m_ndx(n_dx),m_ndt(n_dt)
    {
        double std= vol*sqrt(T);
        m_xmin=std::max(log(spot)-5*std,0.);
        m_xmax=2*log(spot)-m_xmin;
        m_dx=(m_xmax-m_xmin)/double(n_dx-1);
        m_dt=T/double(n_dt-1);

        m_xaxis.resize(n_dx);
        for(int i=0; i<n_dx; i++){
            m_xaxis[i]=m_xmin+i*m_dx;
        }
    }

    std::vector<double> mesh::get_xaxis() const{
        return m_xaxis;
    }

    double mesh::get_ndx() const{
        return m_ndx;
    }

    double mesh::get_ndt() const{
        return m_ndt;
    }

    double mesh::get_xmin() const{
        return m_xmin;
    }

    double mesh::get_xmax() const{
        return m_xmax;
    }

    double mesh::get_T() const{
        return m_T;
    }

    double mesh::get_dx() const{
        return m_dx;
    }

    double mesh::get_dt() const{
        return m_dt;
    }

    double mesh::get_spot() const{
        return m_spot;
    }

}