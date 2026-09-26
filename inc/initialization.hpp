// initialization.hpp
#ifndef INITIALIZATION_HPP
#define INITIALIZATION_HPP

#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <string>
#include <type_traits>
#include <vector>
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
template <typename Func>
inline void l2_projection(std::vector<double>& u_hat,
                          Func u0,
                          const double xl,
                          const double xr,
                          const size_t N,
                          const size_t Nq,
                          const std::string& basis,
                          const std::string& quadrature)
{

    static_assert(std::is_invocable_r_v<double, Func, double>,
                  "Passed function must accept a double and return a double");

#ifndef __OPTIMIZE__
    if (N < 0)
        throw std::invalid_argument("degree of polynomial representing "
                                    "solution u_h cannot be less than 0");
#endif //__OPTIMIZE__

    // ===std::vector -> Eigen::MatrixXd [for Matrix vector calculations]===

    // ===MASS [ (N + 1) x (N + 1) ] ==> M_ml{m->test, l->trial}===
    std::vector<double> M; // symmetric and square
    LOCAL_MASS_MATRIX(M, N, Nq, basis, quadrature);
    Eigen::Map<
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
        M_ml(M.data(), N + 1, N + 1);
    // =========================================================

    // ===VANDERMONDE [ Nq x (N + 1) ] ==> V_qm = { Ψ̂m (ζ𝑞) }===
    std::vector<double> V; // rectangular
    LOCAL_Vqj_MATRIX(V, Nq, N, basis, quadrature);
    Eigen::Map<
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
        V_qm(V.data(), Nq, N + 1);
    // =========================================================

    // ===quadrature weights W [ Nq x Nq ] and ZETA_q [Nq]===
    std::vector<double> wq, ZETA_q; // quadrature weights and abscissae
    if (quadrature == "GL")
        compute_GL_quadrature_weights_and_roots(Nq, wq, ZETA_q);
    if (quadrature == "GLL")
        compute_GLL_quadrature_weights_and_roots(Nq, wq, ZETA_q);
    Eigen::Map<Eigen::VectorXd> w(wq.data(), Nq);
    auto W_qq = w.asDiagonal();
    // =========================================================

    // u0(ZETA_q)
    Eigen::VectorXd U0_q = Eigen::VectorXd::Zero(ZETA_q.size());
    double xq = xl;
    for (size_t q = 0; q < Nq; q++)
    {
        xq = zeta_to_x(ZETA_q[q], xl, xr);
        U0_q[q] = u0(xq);
    }

    // coefficient vector
    u_hat.resize(N + 1);
    Eigen::VectorXd::Map(u_hat.data(), N + 1) =
        M_ml.inverse() * V_qm.transpose() * W_qq * U0_q;
}
#endif // !DEBUG
