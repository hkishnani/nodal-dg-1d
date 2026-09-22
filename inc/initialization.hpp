// initialization.hpp
#ifndef INITIALIZATION_HPP
#define INITIALIZATION_HPP

#include <cmath>
#include <string>
#include <vector>
#include <limits>
#include <eigen3/Eigen/Dense>
#include "quadrature.hpp"
#include "local_matrix_assembly.hpp"
#include "basis.hpp"

/**
 * @brief returns coefficient vector for particular initialization method
 * for piecewise polynomial approximation of Nth degree polynomial
 * @return exact vector of coefficients for given basis
 */

// u_h = sum_{l=0}^{N} u^{hat}_l psi^{hat}_l
// l2_projection in a single element for polynomial of degree N
inline void l2_projection(std::vector<double>& u_hat,
                          const size_t N,
                          const size_t Nq,
                          const std::string& basis,
                          const std::string& quadrature)
{
    // for barycentric lagrange polynomial definition
#ifndef __OPTIMIZE__
    if (N < 0)
        throw std::invalid_argument("degree of polynomial representing "
                                    "solution u_h cannot be less than 0");
#endif //__OPTIMIZE__
    // M_ml -> index :: m->test | l->trial
    
    std::vector<double> M;

    // Vandermonde Matrix [ Nq x (N + 1) ]
    std::vector<double> V_qj;

    // Diagonal Matrix with quadrature weights as Diagonal entries
    std::vector<double> W;



    LOCAL_MASS_MATRIX(M, N, Nq, basis, quadrature);

}

#endif // !DEBUG
