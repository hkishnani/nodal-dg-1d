// linear_advection.cpp
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include "quadrature.hpp"
#include "basis.hpp"
#include "file_handling.hpp"

// Problem specific
constexpr double T0 = 0.0;
constexpr double T = PI;
// [PI, 10 * PI, 100 * PI, 1000 * PI, 2000 * PI]

constexpr double X_L = 0.0;
constexpr double X_R = 2 * PI;

constexpr double a = -2 * PI;
constexpr double lambda = 2.0;
constexpr double l = 2.0 * PI / lambda;

const std::vector<size_t> N_vector = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
// polynomial degree approximated for solution variable u

const std::vector<size_t> K_vector = {2, 4, 8, 16, 32, 64, 128};
// num elements

const size_t K = K_vector[0];
const size_t N = N_vector[0];

const double h = (X_R - X_L) / static_cast<double>(K);
const double delta_t = 1e-5;

// define initialization function
constexpr double U0(double x)
{
    return sin(l * x);
}

// define nodal basis functions evaluation zeta in [-1.0, 1.0]
// Analytical solution U(x, t) --> MOC
constexpr double U(const double x, const double t)
{
    return U0(x - a * t);
}

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


    std::string post_process_folder =
        "/home/himanshu/himanshu/nodal-dg-1d/test/";



    return 0;
}
