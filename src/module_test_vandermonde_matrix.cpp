// module_tests_vandermonde_matrix.cpp

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include "file_handling.hpp"
#include "local_matrix_assembly.hpp"

const std::vector<size_t> N_vector = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
// polynomial degree approximated for solution variable u

std::string Vqj_MATRIX_FOLDER =
    "/home/himanshu/himanshu/nodal-dg-1d/test/vandermonde_matrix/";

int main()
{
    std::cout << std::setprecision(15) << std::scientific;

    std::vector<double> V;
    // R_mq matrix

    std::string vqj_matrix_file_name;

    //=========================
    for (auto&& N : N_vector)
        for (auto&& Nq : {N + 1, N + 2})
            for (auto&& basis : {"GLL", "GL"})
                for (auto&& quadrature : {"GLL", "GL"})
                {
                    vqj_matrix_file_name = Vqj_MATRIX_FOLDER + "Vqj_matrix_N_" +
                                           std::to_string(N) + "_Nq_" +
                                           std::to_string(Nq) + "_basis_" +
                                           basis + "_quadrature_" + quadrature +
                                           ".csv";

                    LOCAL_Vqj_MATRIX(V, Nq, N, basis, quadrature);

                    write_double_precision_matrix_to_csv(
                        vqj_matrix_file_name, V, Nq, N + 1);
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
