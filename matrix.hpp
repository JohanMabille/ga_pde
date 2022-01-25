#ifndef MATRIX_H
#define MATRIX_H
#include <string>
#include <vector>

#include <ostream>
#include <iostream>


namespace dauphine {
class matrix
    {
    public:
        
        matrix(std::size_t nb_rows, std::size_t nb_cols);
        matrix(const matrix&);
        std::size_t nb_rows() const;
        std::size_t nb_cols() const;
        void resize(std::size_t nb_rows, std::size_t nb_cols);
        double& operator()(std::size_t i, std::size_t j);
        const double& operator()(std::size_t i, std::size_t j) const;
        

        matrix& operator+=(const matrix& rhs);
        matrix& operator-=(const matrix& rhs);
        matrix& operator*=(const matrix& rhs);
        matrix& operator/=(const matrix& rhs);
        matrix& operator*=(const std::vector<double>& rhs);
        

        matrix& operator+=(double rhs);
        matrix& operator-=(double rhs);
        matrix& operator*=(double rhs);
        matrix& operator/=(double rhs);

    private:

        std::size_t m_nb_rows;
        std::size_t m_nb_cols;
        std::vector<double> m_data;
    };

    std::ostream& operator<<(std::ostream& out, const matrix& m);

    matrix operator+(const matrix& lhs, const matrix& rhs);
    matrix operator+(const matrix& lhs, double rhs);
    matrix operator+(double lhs, const matrix& rhs);

    std::vector<double> operator*(const matrix& lhs, const std::vector<double>& rhs);

    std::vector<double> operator*(const double& scal, const std::vector<double>& rhs);
    std::vector<double> operator*(const std::vector<double>& lhs, const double& scal);

    std::vector<double> operator+(const double& scal, const std::vector<double>& rhs);
    std::vector<double> operator+(const std::vector<double>& lhs, const double& scal);

    std::vector<double> operator-(const double& scal, const std::vector<double>& rhs);
    std::vector<double> operator-(const std::vector<double>& lhs, const double& scal);

    std::vector<double> operator+(const std::vector<double>& lhs, const std::vector<double>& rhs);
    std::vector<double> operator-(const std::vector<double>& lhs, const std::vector<double>& rhs);

}



#endif // MATRIX_H