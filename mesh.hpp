#ifndef MESH_HPP
#define MESH_HPP

#include <iostream>
#include <vector>


namespace dauphine {

class mesh{
    public: 
    
    mesh(double spot, double T, int n_dx, int n_dt, double vol);

    std::vector<double> get_xaxis() const;
    double get_ndx() const;
    double get_ndt() const;
    double get_xmin() const;
    double get_xmax() const;
    double get_dx() const;
    double get_dt() const;
    double get_T() const;
    double get_spot() const;

    private:
    std::vector<double> m_xaxis;
    double m_spot;
    double m_xmin, m_xmax;
    int m_ndx, m_ndt;
    double m_dx, m_dt;
    double m_T;
};

}

#endif