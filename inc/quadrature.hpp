// quadrature.hpp
/*
In case of Numerical quadrature->
Nq = Number of quadrature points in [-1.0, 1.0]
Nq point rule is exact for any polynomial of degree [2*Nq - 3]
evaluate exact quadrature for polynomial degree of 2*K + 1 = 2*Nq - 3
*/
#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

// legendre polynomial and its derivatives are required to compute GLL
// quadrature weights and roots legendre polynomial is not basis function!
// Atmost we will be required to integrate a polynomial degree upto 2*N_max + 1
// N_max = max polynomial degree to be approximated for solution variable u
// 2*N_max + 1 corresponds to
// Evaluate Legendre polynomial (P_n(zeta)), P'_l(zeta) and P''_l(zeta)
// uisng recurrence relation pg 102 & 103 C.Praveen's notes [17.11]
inline void legendre(size_t l, double zeta, double& P, double& dP, double& ddP)
{
    constexpr double epsilon = std::numeric_limits<double>::epsilon();
    constexpr double tol = 10.0 * epsilon;
    if (l == 0)
    {
        P = 1.0;
        dP = 0.0;
        ddP = 0.0;
        return;
    }
    if (std::fabs(zeta) >= 1.0 - tol)
    {
        P = (zeta > 0.0) ? 1.0 : ((l % 2) ? -1.0 : 1.0);
        dP = (zeta > 0.0) ? static_cast<double>(l * (l + 1)) / 2.0
                          : -P * static_cast<double>(l * (l + 1)) / 2.0;
        ddP = P * static_cast<double>((l - 1) * l * (l + 1) * (l + 2)) * 0.125;
        return;
    }

    double Pnm1 = 1.0; // P_0
    double Pn = zeta;  // P_1

    for (size_t k = 1; k < l; k++)
    {
        const double Pnp1 =
            ((2.0 * k + 1.0) * zeta * Pn - k * Pnm1) / (k + 1.0);

        Pnm1 = Pn;
        Pn = Pnp1;
    }

    P = Pn;
    dP = l * (Pnm1 - zeta * Pn) / (1.0 - zeta * zeta);
    ddP = (2 * zeta * dP - l * (l + 1.0) * Pn) / (1 - zeta * zeta);
}

/**
 * @brief Gauss-Lobatto-Legendre quadrature nodes
 *
 * @param Nq Number of quadrature points in [-1.0, 1.0]
 * @param w quadrature weights
 * @param zeta quadrature abscissae
 *
 * -1.0 & 1.0 always included
 *
 * Nq point rule ==> exact for polynomial of degree upto [2*Nq - 3]
 * @return * void
 */
inline void compute_GLL_quadrature_weights_and_roots(const size_t Nq,
                                                     std::vector<double>& w,
                                                     std::vector<double>& zeta)
{

    constexpr double PI = 3.141592653589793238462643383279502884;
    constexpr double epsilon = std::numeric_limits<double>::epsilon();
    constexpr double tol = 10.0 * epsilon;
    constexpr size_t max_iter = 200;

    if (Nq < 1 || Nq > 18)
        throw std::invalid_argument(
            "QUADRATURE::compute_GLL_quadrature_weights_and_"
            "roots::Nq must satisfy 1 <= Nq <= 18");

    w.resize(Nq);
    zeta.resize(Nq);

    if (Nq == 1)
    {
        w[0] = 2.0;
        zeta[0] = 0.0;
        return;
    }

    zeta[0] = -1.0;
    zeta[Nq - 1] = 1.0;

    const double endpt_wt = 2.0 / (Nq * (Nq - 1.0));

    w[0] = endpt_wt;
    w[Nq - 1] = endpt_wt;

    if (Nq == 2)
        return;

    double P;
    double dP;
    double ddP;

    // Nq - 2 roots in Interior of P'_{Nq - 1}(x) = 0
    for (size_t q = 1; q < Nq - 1; ++q)
    {
        // Chebyshev-Lobatto initial guess
        double zeta_0 =
            -cos(PI * static_cast<double>(q) / static_cast<double>(Nq - 1));

        bool converged = false;

        // qth root of polynomial P'_{Nq - 1}(zeta) = 0
        for (size_t iter = 0; iter < max_iter; iter++)
        {
            legendre(Nq - 1, zeta_0, P, dP, ddP);
            // interior nodes have to satisfy P'_{Nq - 1}(zeta) = 0 for GLL
            // quadrature
            const double dzeta = dP / ddP;
            zeta_0 -= dzeta;
            if (fabs(dzeta) <= tol)
            {
                converged = true;
                break;
            }
        }

        if (!converged)
            throw std::runtime_error("GLL Newton iteration failed to converge");

        zeta[q] = zeta_0;

        legendre(Nq - 1, zeta_0, P, dP, ddP);
        w[q] = 2.0 / (Nq * (Nq - 1.0) * P * P);
    }

    // enforcing symmetry
    for (size_t i = 0; i < Nq / 2; ++i)
    {
        const size_t j = Nq - 1 - i;

        const double zeta_0 = 0.5 * (std::fabs(zeta[i]) + std::fabs(zeta[j]));
        zeta[i] = -zeta_0;
        zeta[j] = zeta_0;

        const double W = 0.5 * (w[i] + w[j]);
        w[i] = W;
        w[j] = W;
    }

    // mid value must be zero
    if (Nq % 2 != 0)
        zeta[Nq / 2] = 0.0;
}

/* GL Accurate upto polynomial of degree 2*Nq - 1
 *
 * zeta in [-1.0, 1.0]
 *
 * Nq = number of quadrature points
 */
inline void compute_GL_quadrature_weights_and_roots(const size_t Nq,
                                                    std::vector<double>& w,
                                                    std::vector<double>& zeta)
{
    if (Nq < 1 || Nq > 18)
        throw std::invalid_argument(
            "QUADRATURE::compute_GL_quadrature_weights_and_"
            "roots::Nq must satisfy 1 <= Nq <= 18");

    constexpr double PI = 3.141592653589793238462643383279502884;
    constexpr double epsilon = std::numeric_limits<double>::epsilon();
    constexpr double tol = 10.0 * epsilon;
    constexpr size_t max_iter = 200;

    w.resize(Nq);
    zeta.resize(Nq);

    if (Nq == 1)
    {
        w[0] = 2.0;
        zeta[0] = 0.0;
        return;
    }

    if (Nq == 2)
    {
        const double z = 1.0 / sqrt(3.0);
        zeta[0] = -z;
        zeta[1] = z;

        w[0] = 1.0;
        w[1] = 1.0;
        return;
    }

    double P;
    double dP;
    double ddP;

    // Nq roots in P_{Nq}(x) = 0 in (-1.0, 1.0)
    for (size_t q = 0; q < Nq; ++q)
    {
        // Tricomi  Approximation initial guess
        double zeta_0 = -cos(PI * static_cast<double>(4.0 * q + 3.0) /
                             static_cast<double>(4.0 * Nq + 2.0));
        bool converged = false;
        // qth root of polynomial P_{Nq}(zeta) = 0
        for (size_t iter = 0; iter < max_iter; iter++)
        {
            legendre(Nq, zeta_0, P, dP, ddP);
            // Interior nodes have to satisfy P_{Nq}(zeta) = 0 for GL quadrature
            const double dzeta = P / dP;
            zeta_0 -= dzeta;
            if (fabs(dzeta) <= tol)
            {
                converged = true;
                break;
            }
        }
        if (!converged)
            throw std::runtime_error("GL Newton iteration failed to converge");

        zeta[q] = zeta_0;
        legendre(Nq, zeta_0, P, dP, ddP);
        w[q] = 2.0 / ((1.0 - zeta_0 * zeta_0) * dP * dP);
    }

    // enforcing symmetry
    for (size_t i = 0; i < Nq / 2; ++i)
    {
        const size_t j = Nq - 1 - i;

        const double zeta_0 = 0.5 * (std::fabs(zeta[i]) + std::fabs(zeta[j]));
        zeta[i] = -zeta_0;
        zeta[j] = zeta_0;

        const double W = 0.5 * (w[i] + w[j]);
        w[i] = W;
        w[j] = W;
    }
    // mid value must be zero
    if (Nq % 2 != 0)
        zeta[Nq / 2] = 0.0;
}

#endif // !QUADRATURE_HPP
// verified on Aug 19 for N upto 17
