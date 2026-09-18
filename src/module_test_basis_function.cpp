// module_test_basis_function.cpp

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include "quadrature.hpp"
#include "basis.hpp"
#include "file_handling.hpp"

const std::vector<size_t> N_vector = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
// polynomial degree approximated for solution variable u

std::string post_process_folder = "/home/himanshu/himanshu/nodal-dg-1d/test/";

std::string PLOT_BASIS_FOLDER = post_process_folder + "plot_basis/";

int main()
{
    std::cout << std::setprecision(15) << std::scientific;

    std::vector<size_t> Nq;

    for (auto&& n : N_vector)
        Nq.push_back(n + 2);

    std::vector<double> zeta_GLL, zeta_GL;
    std::vector<double> w_Lagrange_GL, w_Lagrange_GLL;
    std::vector<double> eta_vector, _;
    double l_j_value_GL = 0.0, l_j_value_GLL = 0.0;

    for (auto&& Nq_i : Nq)
    {
        std::string GL_filename =
            PLOT_BASIS_FOLDER + "GL_basis_Nq_" + std::to_string(Nq_i) + ".csv";

        std::string GLL_filename =
            PLOT_BASIS_FOLDER + "GLL_basis_Nq_" + std::to_string(Nq_i) + ".csv";

        std::cout << "Nq = " << Nq_i << std::endl;

        _.clear();
        zeta_GL.clear();  // GL quadrature roots
        zeta_GLL.clear(); // GLL quadrature roots
        w_Lagrange_GL.clear();
        w_Lagrange_GLL.clear();
        eta_vector.clear();

        // ======compute Lagrange weights at GL and GLL roots======
        compute_GL_quadrature_weights_and_roots(Nq_i, _, zeta_GL);
        compute_GLL_quadrature_weights_and_roots(Nq_i, _, zeta_GLL);

        compute_weights_for_barycentric_lagrange_polynomial(zeta_GL,
                                                            w_Lagrange_GL);

        compute_weights_for_barycentric_lagrange_polynomial(zeta_GLL,
                                                            w_Lagrange_GLL);

        linspace_eta_vector(-1.0, 1.0, 50, eta_vector);
        // ======compute Lagrange weights at GL and GLL roots======

        // write basis function values according to GL or GLL nodes
        write_line_to_file(GL_filename, "n,j,eta,l_j_GL");
        write_line_to_file(GLL_filename, "n,j,eta,l_j_GLL");

        for (size_t j = 0; j < Nq_i; ++j)
            for (double eta : eta_vector)
            {
                l_j_value_GL = l_j(j, zeta_GL, w_Lagrange_GL, eta);
                l_j_value_GLL = l_j(j, zeta_GLL, w_Lagrange_GLL, eta);

                write_line_to_file(GL_filename,
                                   std::to_string(Nq_i) + "," +
                                       std::to_string(j) + "," +
                                       std::to_string(eta) + "," +
                                       std::to_string(l_j_value_GL));

                write_line_to_file(GLL_filename,
                                   std::to_string(Nq_i) + "," +
                                       std::to_string(j) + "," +
                                       std::to_string(eta) + "," +
                                       std::to_string(l_j_value_GLL));
            }

        std::cout << "=========================" << std::endl;
    }

    return 0;
}

/**
 * @file module_test_basis_function.cpp
 * @author himanshu kishnani (himanshukish@iisc.ac.in)
 *
 * @brief compute GL and GLL quadrature weights and roots for Nq_i quadrature
 * points assume Nq_i is the number of solution points available we need to find
 * integration for entries of Mass matrix of degree 2*N integration must be
 * exact for polynomial of degree 2*N + 1
 *
 * @version 1.0
 * @date 2026-09-05
 *
 * @copyright Copyright (c) 2026
 *
 */
