#ifndef rate_hpp
#define rate_hpp

#include <vector>
#include "mesh.hpp"

namespace dauphine

{

    class rate 
    {
        public:
            rate();
            rate(double r);
            rate(std::vector<double> r);
            rate(double r, dauphine::mesh msh,  double alpha);
            std::vector<double> get_rates() const;
            
        protected:
            std::vector<double> p_rate;
    };

}

#endif