#include "PROPOSAL/scattering/stochastic_deflection/bremsstrahlung/BremsGinneken.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/Logging.h"
#include <cmath> 
#include <stdexcept>
#include <algorithm>

using namespace PROPOSAL;
using namespace std;

double stochastic_deflection::BremsGinneken::f_nu_g(double nu, double n, double k_4) const {
    return pow(nu, (1/n + 1)) + pow((0.2/k_4), (1/n)) * nu - pow((0.2/k_4), (1/n));
}

double stochastic_deflection::BremsGinneken::get_nu_g(double E, double n, double k_4) const {
    // auto f = [this, &n, &k_4](double nu) {
    //     return f_nu_g(nu, n, k_4);
    // };

    auto f = [&, n, k_4](double nu) {
        return f_nu_g(nu, n, k_4);
    };

    if (f_nu_g(0.5, n, k_4) * f_nu_g(1., n, k_4) > 0) {
        Logging::Get("proposal.deflection")->warn("Boundaries of Bisection are: {} and {}, E = {}", f_nu_g(0.5, n, k_4), f_nu_g(1., n, k_4), E);
    }

    Logging::Get("proposal.deflection")->warn("Use Bisection: {}, {}", f_nu_g(0.5, n, k_4), f_nu_g(1., n, k_4));
    auto nu_g_x = Bisection(f, 0.5, 1., 1e-5, 100);
    
    return (nu_g_x.first + nu_g_x.second) / 2;
}

double stochastic_deflection::BremsGinneken::get_rms_theta(double e_i, double e_f, double mass, double Z) const {
    auto nu = (e_i - e_f) / (e_i - mass); 
    if (nu <= 0.5) {
        auto k_1 = 0.092 * pow(e_i, -1.0/3.0);
        auto k_2 = 0.052 / e_i * pow(Z, -0.25);
        auto k_3 = 0.22 * pow(e_i, -0.92);
        auto rms_theta = max(min(k_1 * sqrt(nu), k_2), k_3 * nu);
        return rms_theta; 
    } else {
        auto k_4 = 0.26 * pow(e_i, -0.91);
        auto m = 0.5;
        auto d = 1.8;
        auto n = 0.81 * pow(e_i, m) / (pow(e_i, m) + d);
        auto rms_theta = k_4 * pow(nu, 1 + n) * pow(1 - nu, -n);
        if (rms_theta < 0.2) {
            return rms_theta; 
        } else {
            double nu_g = 0.5;
            if (e_i > 0.7) {    // limit of parametrization
                nu_g = get_nu_g(e_i, n, k_4);
            };

            cout << "final nu_g = " << nu_g << ", E = " << e_i << endl;
            auto k_5 = k_4 * pow(nu_g, 1 + n) * pow(1 - nu_g, 0.5 - n);
            rms_theta = k_5 * pow(1 - nu, -0.5);
            if (rms_theta < 0.2) {
                // If no case matches
                Logging::Get("proposal.deflection")->warn("BremsGinneken deflection: rms_theta = {} < 0.2, nu_g = {}", rms_theta, nu_g);
            }
            return rms_theta;
        }
    }
}

UnitSphericalVector 
stochastic_deflection::BremsGinneken::CalculateStochasticDeflection(
    double e_i, double e_f, std::vector<double> const& rnd, size_t component) const
{
    // All energies should be in units of GeV
    e_i = e_i / 1000.0;
    e_f = e_f / 1000.0;
    auto muon_mass = mass / 1000.0;
    // Use Z of medium
    auto Z = Component::GetComponentForHash(component).GetNucCharge(); 

    auto rms_theta = get_rms_theta(e_i, e_f, muon_mass, Z);
    
    auto lambda = 1 / (rms_theta * rms_theta);
    // Need sqrt because of sampling in theta^2
    auto theta_muon = std::sqrt(SampleFromExponential(rnd[0], lambda));
    
    return UnitSphericalVector {theta_muon, 2 * PI * rnd[1]};
}
