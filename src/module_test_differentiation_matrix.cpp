// module_test_differentiation.cpp

#include <iomanip>
#include <ios>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include "file_handling.hpp"
#include "local_matrix_assembly.hpp"

// polynomial degree approximation for solution variable u
const std::vector<size_t> N_vector = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

std::string Dqj_MATRIX_FOLDER =
    "/home/himanshu/himanshu/nodal-dg-1d/test/differentiation_matrix/";

int main()
{
    std::cout << std::setprecision(15) << std::scientific;

    std::vector<double> D;
    // D_qj matrix

    std::string Dqj_matrix_filename;

    //====================
    for (auto&& N : N_vector)
        for (auto&& Nq : {N + 1, N + 2})
            for (auto&& basis : {"GLL", "GL"})
                for (auto&& quadrature : {"GLL", "GL"})
                {
                    Dqj_matrix_filename = Dqj_MATRIX_FOLDER + "Dqj_matrix_N_" +
                                          std::to_string(N) + "_Nq_" +
                                          std::to_string(Nq) + "_basis_" +
                                          basis + "_quadrature_" + quadrature +
                                          ".csv";

                    LOCAL_Dqj_MATRIX(D, Nq, N, basis, quadrature);

                    write_double_precision_matrix_to_csv(
                        Dqj_matrix_filename, D, Nq, N + 1);
                }
    //====================
    return 0;
}
