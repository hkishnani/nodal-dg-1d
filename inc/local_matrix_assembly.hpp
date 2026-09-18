// local_matrix_assembly.hpp
#ifndef LOCAL_MATRIX_ASSEMBLY_HPP
#define LOCAL_MATRIX_ASSEMBLY_HPP

#include <cmath>
// #include <compare>
#include <string>
#include <vector>
// #include <limits>
#include "basis.hpp"
#include "quadrature.hpp"

/**
 * @brief return values for Mass matrix from standard element [-1.0, 1.0] for
 * 0<=j,l<=N Mass matrix only uses basis functions
 * @param N : degree of piecewise polynomial u_h
 *
 * Pre-compute once and save for recursive time step calls
 *
 * for diagonal Local Mass matrix, Nq = N+1
 *
 * but since the basis function are defined on GLL points are GL or GLL
 */
inline void LOCAL_MASS_MATRIX(std::vector<double>& M,
                              const size_t N,
                              const size_t Nq,
                              const std::string& basis,
                              const std::string& quadrature)
{
    // for barycentric lagrange polynomial definition
    std::vector<double>
        w_barycentric; // weights for evaluation of N^th degree polynomial
    std::vector<double>
        zeta_barycentric; // N + 1 GLL/GL nodes in normalized coordinates
    std::vector<double> w_quadrature;    // quadrature weights
    std::vector<double> zeta_quadrature; // quadrature abscissae

    // first get abscissae for lagrange polynomial basis function
    if (basis == "GL")
        compute_GL_quadrature_weights_and_roots(
            N + 1, w_barycentric, zeta_barycentric);

    if (basis == "GLL")
        compute_GLL_quadrature_weights_and_roots(
            N + 1, w_barycentric, zeta_barycentric);

    if (quadrature == "GL")
        compute_GL_quadrature_weights_and_roots(
            Nq, w_quadrature, zeta_quadrature);

    if (quadrature == "GLL")
        compute_GLL_quadrature_weights_and_roots(
            Nq, w_quadrature, zeta_quadrature);

    // barycentric weights need to be evaluated anyways
    w_barycentric.clear();
    compute_weights_for_barycentric_lagrange_polynomial(zeta_barycentric,
                                                        w_barycentric);

    // 0 <= l,m <= N ==> N+1 terms for l and m respectively
    M.resize((N + 1) * (N + 1));

    double val = 0.0;
    for (size_t row = 0; row < (N + 1); row++)
        for (size_t col = 0; col < (N + 1); col++)
            M[row * (N + 1) + col] = 0.0;

    // m_th row i.e. m is fixed along row
    for (size_t m = 0; m <= N; m++)
        for (size_t l = 0; l <= N; l++)
        {
            // l_th column i.e. l is fixed along column
            val = 0.0;
            for (size_t q = 0; q < Nq; q++)
            {
                double zeta_q = zeta_quadrature[q];
                val += w_quadrature[q] *
                       l_j(l, zeta_barycentric, w_barycentric, zeta_q) *
                       l_j(m, zeta_barycentric, w_barycentric, zeta_q);
            }
            M[m * (N + 1) + l] = val;
        }
}

// Vandermonde matrix   --> [ Nq x (N + 1) ]
inline void LOCAL_Vqj_MATRIX(std::vector<double>& V,
                             const size_t Nq,
                             const size_t N,
                             const std::string& basis,
                             const std::string& quadrature)
{
    std::vector<double> w_barycentric, zeta_barycentric;
    std::vector<double> w_quadrature, zeta_quadrature;

    if (basis == "GL")
        compute_GL_quadrature_weights_and_roots(
            N + 1, w_barycentric, zeta_barycentric);

    if (basis == "GLL")
        compute_GLL_quadrature_weights_and_roots(
            N + 1, w_barycentric, zeta_barycentric);

    if (quadrature == "GL")
        compute_GL_quadrature_weights_and_roots(
            Nq, w_quadrature, zeta_quadrature);

    if (quadrature == "GLL")
        compute_GLL_quadrature_weights_and_roots(
            Nq, w_quadrature, zeta_quadrature);

    w_barycentric.clear();
    compute_weights_for_barycentric_lagrange_polynomial(zeta_barycentric,
                                                        w_barycentric);

    // 0 <= q < Nq ==> Nq terms for q (quadrature index)
    // 0 <= j <=N  ==> N + 1 terms for j (index of basis function)
    V.resize(Nq * (N + 1)); // Nq x (N + 1)

    for (size_t q_row = 0; q_row < Nq; q_row++)
        for (size_t j_col = 0; j_col <= N; j_col++)
            V[q_row * (N + 1) + j_col] = 0.0;

    for (size_t q_row = 0; q_row < Nq; q_row++)
        for (size_t j_col = 0; j_col <= N; j_col++)
            V[q_row * (N + 1) + j_col] = l_j(
                j_col, zeta_barycentric, w_barycentric, zeta_quadrature[q_row]);
}

// Differentiation Matrix Dqj = Ψ̂'j (ζ𝑞)    [ Nq x (N + 1) ]
inline void LOCAL_Dqj_MATRIX(std::vector<double>& D,
                             const size_t Nq,
                             const size_t N,
                             const std::string& basis,
                             const std::string& quadrature)
{
    std::vector<double> w_barycentric, zeta_barycentric;
    std::vector<double> w_quadrature, zeta_quadrature;

    if (basis == "GL")
        compute_GL_quadrature_weights_and_roots(
            N + 1, w_barycentric, zeta_barycentric);
    if (basis == "GLL")
        compute_GLL_quadrature_weights_and_roots(
            N + 1, w_barycentric, zeta_barycentric);
    if (quadrature == "GL")
        compute_GL_quadrature_weights_and_roots(
            Nq, w_quadrature, zeta_quadrature);
    if (quadrature == "GLL")
        compute_GLL_quadrature_weights_and_roots(
            Nq, w_quadrature, zeta_quadrature);

    w_barycentric.clear();
    compute_weights_for_barycentric_lagrange_polynomial(zeta_barycentric,
                                                        w_barycentric);

    // 0 <= q < Nq ==> Nq terms for quadrature index
    // 0 <= j <= N ==> N + 1 terms for test function index

    D.resize(Nq * (N + 1)); // Nq x (N + 1)

    // Can zero be treated as a nonsense value?
    for (size_t q_row = 0; q_row < Nq; q_row++)
        for (size_t j_col = 0; j_col <= N; j_col++)
            D[q_row * (N + 1) + j_col] = 0.0;

    for (size_t q_row = 0; q_row < Nq; q_row++)
        for (size_t j_col = 0; j_col <= N; j_col++)
            D[q_row * (N + 1) + j_col] =
                l_j_prime(j_col,
                          zeta_barycentric,
                          w_barycentric,
                          zeta_quadrature[q_row]);
}

// R_{mq} = w_q * Ψ̂_m (ζ_q);    0<=m<=N  0<=q<=Nq
// R_{mq} = matmul(V_{qj}^T , w_q)    [w_j = column vector of weights]
inline void LOCAL_Rmq_MATRIX(std::vector<double>& R,
                             const size_t N,
                             const size_t Nq,
                             const std::string& basis,
                             const std::string& quadrature)
{
    std::vector<double> w_barycentric, zeta_barycentric;
    std::vector<double> w_quadrature, zeta_quadrature;

    if (basis == "GL")
        compute_GL_quadrature_weights_and_roots(
            N + 1, w_barycentric, zeta_barycentric);

    if (basis == "GLL")
        compute_GLL_quadrature_weights_and_roots(
            N + 1, w_barycentric, zeta_barycentric);

    if (quadrature == "GL")
        compute_GL_quadrature_weights_and_roots(
            Nq, w_quadrature, zeta_quadrature);

    if (quadrature == "GLL")
        compute_GLL_quadrature_weights_and_roots(
            Nq, w_quadrature, zeta_quadrature);

    w_barycentric.clear();
    compute_weights_for_barycentric_lagrange_polynomial(zeta_barycentric,
                                                        w_barycentric);

    // 0 <= m <= N ==> N + 1 terms for m (index of test function)
    // 0 <= q < Nq ==> Nq terms for q (quadrature index)
    R.resize((N + 1) * (Nq));

    for (size_t m_row = 0; m_row <= N; m_row++)
        for (size_t q_col = 0; q_col < Nq; q_col++)
            R[m_row * Nq + q_col] = 0.0;

    // m_th row
    for (size_t m_row = 0; m_row <= N; m_row++)
        for (size_t q_col = 0; q_col < Nq; q_col++)
            R[m_row * Nq + q_col] =
                w_quadrature[q_col] * l_j(m_row,
                                          zeta_barycentric,
                                          w_barycentric,
                                          zeta_quadrature[q_col]);
}

#endif
