// module_tests.cpp

/**
 * ======compute GL and GLL quadrature weights and roots for Nq_i=======
 * assume Nq_i is the number of solution points available we need to find
 * integration for entries of Mass * matrix of degree 2*N integration must be
 * exact for polynomial of degree 2*N + 1
 */

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include "quadrature.hpp"
#include "basis.hpp"
#include "file_handling.hpp"

const std::vector<size_t> N_vector = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
// polynomial degree approximated for solution variable u

std::string QUADRATURE_WEIGHTS_ABSCISSAE_FOLDER =
    "/home/himanshu/himanshu/nodal-dg-1d/test/quadrature_weights_abscissae/";

int main()
{
    std::cout << std::setprecision(15) << std::scientific;

    std::vector<size_t> Nq;

    for (auto&& n : N_vector)
        Nq.push_back(n + 2);

    std::vector<double> w_GLL, zeta_GLL;
    std::vector<double> w_GL, zeta_GL;

    double l_j_value_GL = 0.0, l_j_value_GLL = 0.0;

    for (auto&& Nq_i : Nq)
    {
        std::cout << "Nq = " << Nq_i << std::endl;

        // GL quadrature weights and roots
        w_GL.clear();
        zeta_GL.clear();
        compute_GL_quadrature_weights_and_roots(Nq_i, w_GL, zeta_GL);
        write_csv(QUADRATURE_WEIGHTS_ABSCISSAE_FOLDER + "GL_abscissae_weights" +
                      std::to_string(Nq_i) + ".csv",
                  zeta_GL,
                  w_GL,
                  "zeta_GL",
                  "w_GL"); // WRITE

        // GLL quadrature weights and roots
        w_GLL.clear();
        zeta_GLL.clear();
        compute_GLL_quadrature_weights_and_roots(Nq_i, w_GLL, zeta_GLL);
        write_csv(QUADRATURE_WEIGHTS_ABSCISSAE_FOLDER +
                      "GLL_abscissae_weights" + std::to_string(Nq_i) + ".csv",
                  zeta_GLL,
                  w_GLL,
                  "zeta_GLL",
                  "w_GLL"); // WRITE

        std::cout << "=========================" << std::endl;
    }
    //=========================

    return 0;
}