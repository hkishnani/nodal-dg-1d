// module_test_L2_projection.cpp
#include "basis.hpp"
#include "file_handling.hpp"
#include "initialization.hpp"
#include "quadrature.hpp"
#include <iomanip>
#include <ios>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#define N_SAMPLE 50

constexpr double PI = 3.141592653589793238462643383279502884;

const std::string basis = "GLL", quadrature = "GLL";
// GLL or GL

// Problem specific
const double X_0 = 0.0;
const double X_L = 2 * PI;

// wavenumber ==> lambda
constexpr inline double U0_func(double x)
{
    double lambda = 2.0 * PI;
    return sin(x * 2.0 * PI / lambda);
}

const std::string L2_projection_folder = "/home/himanshu/himanshu/nodal-dg-1d/"
                                         "test/initialization/l2_projection/";
int main(int argc, char* argv[])
{
    // runtime checks
    if (argc < 3)
    {
        std::cout << "Usage: " << argv[0]
                  << " <K [number of elements]> <N [Polynomial degree for "
                     "approximation of u_h]>"
                  << std::endl;
        return 1;
    }

    size_t K = static_cast<size_t>(std::stoull(argv[1]));
    size_t N = static_cast<size_t>(std::stoull(argv[2]));
    size_t Nq = (quadrature == "GLL") ? (N + 2) : (N + 1);
    // Number of quadrature pts for M, V, zeta_q etc.

    if (N < 0 || N > 11 || K < 1)
    {
        std::cout << "Try with N >= 0 and K >= 1  and N <= 11" << std::endl;
        return 1;
    }

    std::vector<double> w_barycentric, zeta_barycentric;
    if (basis == "GL")
        compute_GL_quadrature_weights_and_roots(
            N + 1, w_barycentric, zeta_barycentric);

    if (basis == "GLL")
        compute_GLL_quadrature_weights_and_roots(
            N + 1, w_barycentric, zeta_barycentric);

    w_barycentric.clear();
    compute_weights_for_barycentric_lagrange_polynomial(zeta_barycentric,
                                                        w_barycentric);

    std::vector<double> X, X_e;
    linspace_eta_vector(X_0, X_L, K + 1, X);

    // full domain discretization

    std::vector<double> u_hat; // Coefficient vector

    // writing error and its L2 norm
    size_t Nq_L2 = 18;
    std::vector<double> zeta_quadrature, w_quadrature;
    double error_l2_norm = 0.0;
    compute_GLL_quadrature_weights_and_roots(
        Nq_L2, w_quadrature, zeta_quadrature);

    std::cout << "i, x, uh, u0" << std::endl;

    for (size_t i = 0; i < K; i++)
    {
        double xl = X[i], xr = X[i + 1];
        double delta_x = xr - xl;
        // left and right end of element k
        l2_projection(u_hat, U0_func, xl, xr, N, Nq, basis, quadrature);

        // for plotting
        linspace_eta_vector(xl, xr, N_SAMPLE, X_e);
        for (size_t n = 0; n < N_SAMPLE; n++)
        {
            double x = X_e[n];
            double eta = x_to_zeta(x, xl, xr);
            double uh = evaluate_interpolated_function_value(
                N, zeta_barycentric, w_barycentric, u_hat, eta);
            std::cout << i << "," << x << "," << uh << "," << U0_func(x)
                      << std::endl;
        }

        for (size_t nq = 0; nq < Nq_L2; nq++)
        {
            // Piecewise Polynomial function evaluation
            double zeta_q = zeta_quadrature[nq];
            double uh_q = evaluate_interpolated_function_value(
                N, zeta_barycentric, w_barycentric, u_hat, zeta_q);

            double xq = zeta_to_x(zeta_q, xl, xr);
            double u0_q = U0_func(xq);

            error_l2_norm += w_quadrature[nq] * (u0_q - uh_q) * (u0_q - uh_q);
        }
        error_l2_norm = error_l2_norm * (delta_x / 2.0);
    }
    std::string err_line =
        std::to_string(K) + "," + std::to_string(N) + "," + std::to_string(error_l2_norm);

    write_line_to_file(L2_projection_folder + "l2_norm_vs_h.csv",
                       err_line);
    return 0;
}
