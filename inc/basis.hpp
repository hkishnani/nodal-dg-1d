// basis.hpp
#ifndef BASIS_HPP
#define BASIS_HPP

#include <cmath>
#include <vector>
#include <limits>
#include "lin_alg.hpp"

// #define PI 3.141592653589793238462643383279502884

//===========reference cell calculations==============
// [xl, xr] -> [-1.0, 1.0]
inline constexpr double
x_to_zeta(const double x, const double xl, const double xr)
{
    return (2 * x - (xr + xl)) / (xr - xl);
}

inline constexpr double
zeta_to_x(const double zeta, const double xl, const double xr)
{
    return 0.5 * ((1.0 - zeta) * xl + (1.0 + zeta) * xr);
}

// Linearly spaced vector on [-1, 1]
inline void linspace_eta_vector(const double x0,
                                const double xL,
                                const size_t N,
                                std::vector<double>& eta)
{
    eta.resize(N);
    eta.clear();
    if (N == 1)
    {
        eta[0] = (x0 + xL) * 0.5;
        return;
    }

    for (size_t i = 0; i < N; ++i)
        eta[i] = x0 + (xL - x0) * static_cast<double>(i) /
                          static_cast<double>(N - 1);
}

//===========reference cell calculations============
// define nodal basis functions evaluation zeta in [-1.0, 1.0]
// Barycentric Lagrange form [C.Praveen notes Appendix-B]
// Evaluation of weights on zeta[n] --> test once
inline void compute_weights_for_barycentric_lagrange_polynomial(
    const std::vector<double>& ZETA,
    std::vector<double>& w)
{
    // 0 <= j <= (n-1)
    double t = 0.0;
    size_t n = ZETA.size();
    w.resize(n);
    for (size_t i = 0; i < n; ++i)
        w[i] = 1.0;

    for (size_t i = 1; i < n; ++i)
        for (size_t j = 0; j < i; ++j)
        {
            t = ZETA[i] - ZETA[j];
            w[i] *= t;
            w[j] *= -t;
        }

    for (size_t i = 0; i < n; ++i)
        w[i] = 1.0 / w[i];
} // verification done Aug 19 2026

// FOR LAGRANGE POLYNOMIAL OF DEGREE n-1 =>
// p(eta) = sum[j=0 to n] l_j(eta) * f_j
// l_j(eta) = ( l(eta) * w_j ) / (eta - zeta_j)
// w = Barycentric weights, ZETA = Barycentric abscissae
inline double l_j(const size_t j,
                  const std::vector<double>& ZETA_barycentric,
                  const std::vector<double>& w_barycentric,
                  const double eta)
{
    const double eps = std::numeric_limits<double>::epsilon();
    // evaluation at nodal point = 1.0
    if (fabs(eta - ZETA_barycentric[j]) < eps)
        return 1.0;

    size_t N = ZETA_barycentric.size() - 1;
    // evaluation at an arbitrary point
    double l_eta = 1.0;
    for (size_t i = 0; i <= N; ++i)
        if (i != j)
            l_eta *= (eta - ZETA_barycentric[i]);

    return (w_barycentric[j] * l_eta);
}
// verified Aug 20 2026

// at zeta_i
inline double l_j_prime_at_zeta_i(const size_t j,
                                  const size_t i,
                                  const std::vector<double>& ZETA_barycentric,
                                  const std::vector<double>& w_barycentric)
{
    size_t N = ZETA_barycentric.size() - 1;
    //  if(i != j)
    //      return (w[j] / w[i]) / (ZETA[i] - ZETA[j]);

    if (i == j)
    {
        double t = 0.0;
        for (size_t k = 0; k <= N; k++)
            if (k != i)
                t -= (w_barycentric[k] / w_barycentric[i]) /
                     (ZETA_barycentric[i] - ZETA_barycentric[k]);

        return t;
    }

    return (w_barycentric[j] / w_barycentric[i]) /
           (ZETA_barycentric[i] - ZETA_barycentric[j]);
}

// computing l_prime on arbitrary zeta value
inline double l_prime(const std::vector<double>& ZETA_barycentric,
                      const double eta)
{
    double val = 0.0;
    size_t N = ZETA_barycentric.size() - 1;

    for (size_t m = 0; m <= N; m++)
    {
        double val_k = 1.0;

        for (size_t k = 0; k <= N; k++)
            if (m != k)
                val_k = val_k * (eta - ZETA_barycentric[k]);

        val += val_k;
    }

    return val;
}

// generic call --> diverts to l_j_prime_at_zeta_i if eta == zeta
inline double l_j_prime(const size_t j,
                        const std::vector<double>& ZETA_barycentric,
                        const std::vector<double>& w_barycentric,
                        const double eta)
{
    const double eps = std::numeric_limits<double>::epsilon();

    size_t N = ZETA_barycentric.size() - 1;

    // if eta is one of the Barycentric coordinate
    for (size_t i = 0; i <= N; i++)
        if (fabs(eta - ZETA_barycentric[i]) < eps)
            return l_j_prime_at_zeta_i(
                j, i, ZETA_barycentric, w_barycentric);

    double l_prime_val = l_prime(ZETA_barycentric, eta);

    double l_j_val = l_j(j, ZETA_barycentric, w_barycentric, eta);

    double val = (w_barycentric[j] * l_prime_val - l_j_val) / (eta - ZETA_barycentric[j]);

    return val;
}

// ==========================================================
// eval. Barycentric Lagrange interpolated function p(eta) from f:[-1, +1] -> R
// N = dim(p(eta))
inline double
evaluate_interpolated_function_value(const size_t N,
                                     const std::vector<double>& ZETA_barycentric,
                                     const std::vector<double>& w_barycentric,
                                     const std::vector<double>& f_hat,
                                     const double eta)
{
    double num = 0.0, den = 0.0, t = 0.0;
    const double eps = std::numeric_limits<double>::epsilon();

    // for p(eta) of degree n-1, we need n-1 basis functions hence
    for (size_t j = 0; j <= N; ++j)
    {
        if (fabs(eta - ZETA_barycentric[j]) < eps)
            return f_hat[j];

        t = (w_barycentric[j] / (eta - ZETA_barycentric[j]));
        num += t * f_hat[j];
        den += t;
    }

    return num / den;
}

// dim(p_prime(eta)) = N + 1
// ZETA[i] is known
inline double evaluate_interpolated_function_derivative_at_given_node(
    const size_t i,
    const std::vector<double>& ZETA_barycentric,
    const std::vector<double>& w_barycentric,
    const std::vector<double>& f)
{
    const size_t Nq = ZETA_barycentric.size();
    std::vector<double> lj_prime(Nq);
    for (size_t j = 0; j < Nq; j++)
        lj_prime[j] =
            l_j_prime_at_zeta_i(j, i, ZETA_barycentric, w_barycentric);

    double t = 0.0;
    for (size_t j = 0; j < Nq; ++j)
        t += lj_prime[j] * f[j];

    return t;
}

// Computing derivatives of constructed interpolating polynomial
inline double
evaluate_interpolated_function_derivative(const size_t N,
                                          const std::vector<double>& ZETA,
                                          const std::vector<double>& w,
                                          const std::vector<double>& f,
                                          const double eta)
{
    const double eps = std::numeric_limits<double>::epsilon();

    for (size_t j = 0; j <= N; ++j)
        if (fabs(eta - ZETA[j]) < eps)
            return evaluate_interpolated_function_derivative_at_given_node(
                j, ZETA, w, f);

    double num = 0.0, den = 0.0, t = 0.0, p_eta = 0.0;
    p_eta = evaluate_interpolated_function_value(N, ZETA, w, f, eta);
    for (size_t j = 0; j <= N; j++)
    {
        t = (w[j] / (eta - ZETA[j]));
        num += (t * ((p_eta - f[j]) / (eta - ZETA[j])));
        den += t;
    }

    return num / den;
}

// return p_prime = [D] {f} at all nodal points
// this function is doubtful
inline void evaluate_interpolated_function_derivative_at_all_nodes(
    const size_t N,
    const std::vector<double>& ZETA,
    const std::vector<double>& w,
    const std::vector<double>& f,
    std::vector<double>& p_prime)
{
    double t = 0.0;
    std::vector<double> D((N + 1) * (N + 1)); // differentiation matrix

    // Compute [D] = D_ij
    for (size_t i = 0; i <= N; ++i)
    {
        t = 0.0;
        for (size_t j = 0; j <= N; ++j)
        {
            if (i != j)
            {
                D[i * (N + 1) + j] = (w[j] / w[i]) / (ZETA[i] - ZETA[j]);
                t += D[i * (N + 1) + j];
            }
        }
        // diagonal entries
        D[i * (N + 1) + i] = -t;
    }

    // Do [D] {f}
    Ax(N + 1, N + 1, D, f, p_prime);
}

#endif
