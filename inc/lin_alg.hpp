// lin_alg.hpp
#ifndef LIN_ALG_HPP
#define LIN_ALG_HPP

#include <vector>
#include <stdexcept>
#include <eigen3/Eigen/Dense>
// #include <iostream>

// Matrix-Vector Product-> A [m rows n columns] x[n rows 1 column]
// always write as flattened matrix -->
// read as column major if transpose_A == true
inline void Ax(const size_t m,
               const size_t n,
               const std::vector<double>& A,
               const std::vector<double>& x,
               std::vector<double>& b)
{

#ifndef __OPTIMIZE__
    if (x.size() != n)
        throw std::invalid_argument(
            "Ax error: x size (" + std::to_string(x.size()) +
            ") must match columns n (" + std::to_string(n) + ").");

    if (A.size() != m * n)
        throw std::invalid_argument(
            "Ax error: A size (" + std::to_string(A.size()) +
            ") must equal m * n (" + std::to_string(m * n) + ").");
#endif

    // read as row major
    b.resize(m);
    // i^th row and j^th column
    for (size_t i = 0; i < m; ++i)
    {
        b[i] = 0.0;
        for (size_t j = 0; j < n; ++j)
            b[i] += A[i * n + j] * x[j];
    }
}

inline void ATx(const size_t m,
                const size_t n,
                const std::vector<double>& A,
                const std::vector<double>& x,
                std::vector<double>& b)
{

#ifndef __OPTIMIZE__
    if (x.size() != m)
        throw std::invalid_argument(
            "ATx error: x size (" + std::to_string(x.size()) +
            ") must match rows m (" + std::to_string(m) + ").");

    if (A.size() != m * n)
        throw std::invalid_argument(
            "ATx error: A size (" + std::to_string(A.size()) +
            ") must equal m * n (" + std::to_string(m * n) + ").");
#endif

    b.resize(n);
    // A = m rows x n cols
    for (size_t i = 0; i < n; i++)
    {
        b[i] = 0.0;
        for (size_t j = 0; j < m; j++)
            b[i] += A[j * n + i] * x[j];
    }
}

#endif // LIN_ALG_HPP

// keep it bare bones
