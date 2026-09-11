// module_tests.cpp
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include "quadrature.hpp"
#include "basis.hpp"
#include "file_handling.hpp"
#include "local_matrix_assembly.hpp"

const std::vector<size_t> N_vector = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
// polynomial degree approximated for solution variable u

std::string MASS_MATRIX_FOLDER =
    "/home/himanshu/himanshu/nodal-dg-1d/test/mass_matrix/";

int main()
{
    std::cout << std::setprecision(15) << std::scientific;

    std::vector<double> M;
    // Mass matrix

    std::string mass_matrix_file_name;

    //=========================
    // degree of polynomial basis function
    for (auto&& N : N_vector)
    {
        // Nq quad. pts for integration of Mass matrix entries
        for (auto&& Nq: {N + 1, N + 2})
        {
            // Lagrange basis functions are defined on GLL or GL point
            for (auto&& basis : {"GLL", "GL"})
            {
                // quadrature rule for integration of Mass matrix entries
                for (auto&& quadrature : {"GLL", "GL"})
                {
                    mass_matrix_file_name =
                        MASS_MATRIX_FOLDER + "mass_matrix_N_" +
                        std::to_string(N) + "_Nq_" + std::to_string(Nq) +
                        "_basis_" + basis + "_quadrature_" + quadrature +
                        ".csv";

                    LOCAL_MASS_MATRIX(M, N, Nq, basis, quadrature);

                    write_double_precision_matrix_to_csv(
                        mass_matrix_file_name, M, N + 1, N + 1);
                }
            }
        }
    }
    //=========================
    return 0;
}

/**
 * ======compute GL and GLL quadrature weights and roots for Nq_i=======
 * assume Nq_i is the number of solution points available we need to find
 * integration for entries of Mass * matrix of degree 2*N integration must be
 * exact for polynomial of degree 2*N + 1
 */
