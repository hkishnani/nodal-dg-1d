// linear_advection_steady_state.cpp
// f = au
// df / dx = g(x)
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <cmath>
#include <iterator>
#include <memory>
#include <string>
#include <vector>
#include <format>
#include <iomanip>
#include <ios>
#include "quadrature.hpp"
#include "basis.hpp"
#include "file_handling.hpp"
#include "initialization.hpp"
#include "local_matrix_assembly.hpp"
#include "numerical_flux.hpp"

#define N_SAMPLE 50

// Constants
constexpr double PI = 3.141592653589793238462643383279502884;

const std::string basis = "GLL", quadrature = "GLL"; // GL/GLL

// Problem specific
constexpr double X_0 = 0.0, X_L = 2 * PI;
constexpr double a = 1.0; // constant coefficient

// define initialization function ∿
constexpr double lambda = 2.0;
constexpr double l = 2.0 * PI / lambda;

// U_analytic solution
constexpr inline double U_analytic(double x)
{
    return sin(l * x);
}

// Right hand side function
constexpr inline double g_RHS(double x)
{
    return a * l * cos(l * x);
}

const std::string steady_state_linear_advection_folder =
    std::filesystem::current_path().string();

int main(int argc, char* argv[])
{
    std::cout << std::setprecision(15) << std::scientific;

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
    // Number of quadrature points for M, V, zeta_q etc.

    if (N < 0 || N > 11 || K < 1)
    {
        std::cout << "Try with 0 <= N <= 11 and K >= 1" << std::endl;
        return 1;
    }
    std::string K_N_filename =
        "/K" + std::to_string(K) + "_N" + std::to_string(N) + ".csv";

    // BARYCENTRIC ZETA
    std::vector<double> w_barycentric, zeta_barycentric;
    if (basis == "GL")
        compute_GL_quadrature_weights_and_roots(
            N + 1, w_barycentric, zeta_barycentric);

    if (basis == "GLL")
        compute_GLL_quadrature_weights_and_roots(
            N + 1, w_barycentric, zeta_barycentric);

    // BARYCENTRIC WEIGHTS
    w_barycentric.clear();
    compute_weights_for_barycentric_lagrange_polynomial(zeta_barycentric,
                                                        w_barycentric);

    // DOMAIN DISCRETIZATION
    std::vector<double> X, Xe;
    linspace_eta_vector(X_0, X_L, K + 1, X);

    // ELEMENTWISE COEFFICIENT VECTOR
    std::vector<double> u_hat(N+1, 0.0);

    const double h = (X_L - X_0) / static_cast<double>(K);

    // Writing error and its L2 norm
    size_t Nq_L2 = 18;
    std::vector<double> zeta_quadrature, w_quadrature;
    double error_l2_norm = 0.0;
    compute_GLL_quadrature_weights_and_roots(
        Nq_L2, w_quadrature, zeta_quadrature);

    // iteration on elements
    for (size_t i = 0; i < K; i++)
    {
        double xl = X[i], xr = X[i + 1];
        double delta_x = xr - xl;
        // left and right end of element k
        // Solving system of linear equation here:
        // [S_jq]{f(u_iq)} + flux-terms = [V_qj]^T [W_qq] {g(𝛇_q)}
        // solve for u_hat

        // -------------------
        // for plotting
        linspace_eta_vector(xl, xr, N_SAMPLE, Xe);
        for (size_t n = 0; n < N_SAMPLE; n++)
        {
            double x = Xe[n];
            double eta = x_to_zeta(x, xl, xr);
            double uh = evaluate_interpolated_function_value(
                N, zeta_barycentric, w_barycentric, u_hat, eta);
            double u_val = U_analytic(x);
            std::string i_x_uh_u = std::to_string(i) + "," +
                                   std::format("{:.15f}", x) + "," +
                                   std::format("{:.15f}", uh) + "," +
                                   std::format("{:.15f}", u_val);

            write_line_to_file(
                steady_state_linear_advection_folder + K_N_filename, i_x_uh_u);
        }

        // for L2 norm calculation within each element
        double error_l2_norm_i = 0.0;
        for (size_t nq = 0; nq < Nq_L2; nq++)
        {
            // Piecewise Polynomial function evaluation
            double zeta_q = zeta_quadrature[nq],
                   uh_q = evaluate_interpolated_function_value(
                       N, zeta_barycentric, w_barycentric, u_hat, zeta_q),
                   xq = zeta_to_x(zeta_q, xl, xr),
                   u_q = U_analytic(xq);

            error_l2_norm_i += w_quadrature[nq] * (u_q - uh_q) * (u_q - uh_q);
        }
        error_l2_norm += error_l2_norm_i * (delta_x / 2.0);
    }
    std::string err_line = std::to_string(K) + "," + std::to_string(N) + "," +
                           std::format("{:.15f}", error_l2_norm);

    write_line_to_file(
        steady_state_linear_advection_folder + "/l2_norm_vs_h.csv", err_line);
    return 0;
}