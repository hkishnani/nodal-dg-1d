// lin_alg.hpp
#ifndef LIN_ALG_HPP
#define LIN_ALG_HPP

#include <vector>

// Matrix-Vector Product-> A [m rows n columns] x[n rows 1 column]
inline void Ax(const size_t m,
        const size_t n,
        const std::vector<double>& A,
        const std::vector<double>& x,
        std::vector<double>& b)
{
    b.resize(m);
    // i^th row and j^th column
    for (size_t i = 0; i < m; ++i)
    {
        b[i] = 0.0;
        for (size_t j = 0; j < n; ++j)
            b[i] += A[i * n + j] * x[j];
    }
}
#endif // LIN_ALG_HPP
