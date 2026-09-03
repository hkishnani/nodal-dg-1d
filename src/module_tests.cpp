// module_tests.cpp
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include "quadrature.hpp"
#include "basis.hpp"
#include "file_handling.hpp"

const std::vector<size_t> N_vector = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
// polynomial degree approximated for solution variable u

int main()
{
    std::cout << std::setprecision(15) << std::scientific;

    // Compute Lagrange weights and its values at GLL roots for inspection for
    // polynomial degree upto 2*N_max + 1
    std::vector<size_t> Nq;
    // const size_t N_max = *std::max_element(N_vector.begin(), N_vector.end());

    // check for polynomial degree upto 2*N_max + 1 ==> Nq = N + 2
    for (auto&& n : N_vector)
        Nq.push_back(n + 2);
    // 2Nq - 3 = 2N + 1 polynomial degree to be checked for quadrature

    std::vector<double> w_GLL, zeta_GLL;
    std::vector<double> w_GL, zeta_GL;
    std::vector<double> w_Lagrange_GL, w_Lagrange_GLL;
    std::vector<double> eta_vector;
    double l_j_value_GL = 0.0, l_j_value_GLL = 0.0;

    std::string post_process_folder =
        "/home/himanshu/himanshu/nodal-dg-1d/test/";

    for (auto&& Nq_i : Nq)
    {
        // assume Nq_i is the number of solution points available
        // we need to find integration for entries of Mass matrix of degree 2*N
        // integration must be exact for polynomial of degree 2*N + 1
        std::cout << "Nq = " << Nq_i << std::endl;

        w_GL.clear();
        zeta_GL.clear();
        w_GLL.clear();
        zeta_GLL.clear();
        eta_vector.clear();
        w_Lagrange_GL.clear();
        w_Lagrange_GLL.clear();

        // pass number of quadrature points Nq_i
        compute_GL_quadrature_weights_and_roots(Nq_i, w_GL, zeta_GL);
        compute_GLL_quadrature_weights_and_roots(Nq_i, w_GLL, zeta_GLL);

        // write quadrature weights and abscissae
        write_csv(post_process_folder +
                      "quadrature_weights_abscissae/GL_abscissae_weights" +
                      std::to_string(Nq_i) + ".csv",
                  zeta_GL,
                  w_GL,
                  "zeta_GL",
                  "w_GL");

        write_csv(post_process_folder +
                      "quadrature_weights_abscissae/GLL_abscissae_weights" +
                      std::to_string(Nq_i) + ".csv",
                  zeta_GLL,
                  w_GLL,
                  "zeta_GLL",
                  "w_GLL");

        std::cout << "=========================" << std::endl;

        // compute Lagrange weights at GL and GLL roots
        w_Lagrange_GL.clear();
        w_Lagrange_GLL.clear();
        compute_weights(zeta_GL, w_Lagrange_GL);
        compute_weights(zeta_GLL, w_Lagrange_GLL);

        linspace_eta_vector(-1.0, 1.0, 50, eta_vector);

        std::string GL_filename = post_process_folder +
                                  "/plot_basis/GL_basis_Nq_" +
                                  std::to_string(Nq_i) + ".csv";

        std::string GLL_filename = post_process_folder +
                                   "/plot_basis/GLL_basis_Nq_" +
                                   std::to_string(Nq_i) + ".csv";

        // calculate l_j(eta) for j = 0, 1, ..., zeta_*.size()
        // eta in [-1.0, 1.0]

        // write basis function values according to GL nodes
        write_line_to_file(GL_filename, "n,j,eta,l_j_GL");
        for (size_t j = 0; j < Nq_i; ++j)
            for (double eta : eta_vector)
            {
                l_j_value_GL = l_j(j, Nq_i, zeta_GL, w_Lagrange_GL, eta);
                write_line_to_file(GL_filename,
                                   std::to_string(Nq_i) + "," +
                                       std::to_string(j) + "," +
                                       std::to_string(eta) + "," +
                                       std::to_string(l_j_value_GL));
            }

        // writing basis function values according to GLL nodes
        write_line_to_file(GLL_filename, "n,j,eta,l_j_GLL");
        for (size_t j = 0; j < Nq_i; ++j)
            for (double eta : eta_vector)
            {
                l_j_value_GLL =
                    l_j(j, zeta_GLL.size(), zeta_GLL, w_Lagrange_GLL, eta);

                write_line_to_file(GLL_filename,
                                   std::to_string(Nq_i) + "," +
                                       std::to_string(j) + "," +
                                       std::to_string(eta) + "," +
                                       std::to_string(l_j_value_GLL));
            }
    }

    return 0;
}
