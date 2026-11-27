// numerical_flux.hpp
// f( uhl(xi + 0.5), uhr(xi + 0.5) )
#include "basis.hpp"
#include <cmath>
#include <string>
#include <vector>

// f(uh) = a * uh for linear advection flux
double upwind_flux(const double uh_L, const double uh_R, const double a, auto f)
{
    return 0.5 * (f(uh_L) + f(uh_R)) - 0.5 * std::fabs(a) * (uh_R - uh_L);
}

double lax_friedrich_flux(const double uh_L,
                          const double uh_R,
                          const double delta_x,
                          const double delta_t,
                          auto f)
{
    return 0.5 * (f(uh_L) + f(uh_R)) -
           0.5 * std::fabs(delta_x / delta_t) * (uh_R - uh_L);
    // for steady state, we can choose delta_t = delta_x / a
    // for steady state, this is equal to upwind_flux
}

// interface flux at f face in F_vector
void interface_flux(std::vector<double>& F_vector,
                    auto f_func,
                    const std::vector<double>& u_hat_L,
                    const std::vector<double>& u_hat_R,
                    const std::vector<double>& w_barycentric,
                    const std::vector<double>& zeta_barycentric,
                    const std::string& basis,
                    const double delta_x,
                    const double delta_t,
                    const double a,
                    const std::string& flux_name)
{
    F_fvector.clear();
    F_fvector.resize(zeta_barycentric.size());
    double uh_L = 0.0;
    double uh_R = 0.0;

    for (size_t j = 0; j < zeta_barycentric.size(); j++)
    {
        F_fvector[j] = l_j(j, zeta_barycentric, w_barycentric, )
    }
}
