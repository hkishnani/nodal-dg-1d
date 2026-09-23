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
#define N_SAMPLE 10

constexpr double PI = 3.141592653589793238462643383279502884;

const std::string basis = "GLL";      // GLL or GL
const std::string quadrature = "GLL"; // GLL or GL

// Problem specific
const double X_0 = 0.0;
const double X_L = 2 * PI;

// wavenumber ==> lambda
constexpr inline double U0_func(double x)
{
    double lambda = 2.0;
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
                     "approximation of u_h]> \n";
        return 1;
    }

    size_t K = static_cast<size_t>(std::stoull(argv[1]));
    size_t N = static_cast<size_t>(std::stoull(argv[2]));
    size_t Nq = (quadrature == "GLL") ? (N + 2) : (N + 1);
    // Number of quadrature pts for M, V, zeta_q etc.

    if (N < 0 || N > 11 || K < 1)
    {
        std::cout << "Try with N >= 0 and K >= 1  and N <= 11";
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

    std::vector<double> X;
    linspace_eta_vector(X_0, X_L, K, X);

    std::cout << "i, x, uh, u0" << std::endl;
    for (size_t i = 0; i < K; i++)
    {
        std::vector<double> u_hat(N + 1, 0.0);
        // Coefficient vector

        double xl = X[i], xr = X[i + 1];
        // left and right end of element k

        l2_projection(u_hat, U0_func, xl, xr, N, Nq, basis, quadrature);

        std::vector<double> X_e(N_SAMPLE, 0.0);
        // element values at sample locations X_e

        std::vector<double> Uh_values(N_SAMPLE, 0.0);
        // Uh == Approx. Piecewise Poly evaluated at X_e

        std::vector<double> U0_values(N_SAMPLE, 0.0);
        // U0 == analytical function

        linspace_eta_vector(xl, xr, N_SAMPLE, X_e);

        for (size_t n = 0; n < N_SAMPLE; n++)
        {
            double eta = x_to_zeta(X_e[n], xl, xr);

            Uh_values[n] = evaluate_interpolated_function_value(
                N, zeta_barycentric, w_barycentric, u_hat, eta);

            U0_values[n] = U0_func(X_e[n]);

            std::cout << i << "," << X_e[n] << "," << Uh_values[n] << ","
                      << U0_values[n] << std::endl;
        }
    }

    return 0;
}
