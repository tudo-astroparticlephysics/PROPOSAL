#include "PROPOSAL/scattering/stochastic_deflection/bremsstrahlung/TsaiParametrizationBremsstrahlungGinneken.h"
#include "PROPOSAL/Constants.h"
#include <cmath> 
#include<stdexcept>
using namespace PROPOSAL;

double get_nu_g(double E, double x_min = 0.00005) {
    auto n = 0.81 * pow(E, 0.5) / (pow(E, 0.5) + 1.8);
    auto k_4 = 0.26 * pow(E, -0.91);
    auto v = logspace(log10(x_min), log10(0.5), 20);
    auto nu_lin = reduce_array(v, 1);
    auto f_nus = f_nu_g(nu_lin, n, k_4);
    auto nu_g = nu_lin[distance(f_nus.begin(), min_element(f_nus.begin(), f_nus.end()))];
    
    return nu_g;
}

std::vector<double> f_nu_g(std::vector<double> nu, double n, double k_4) {
    std::vector<double> nus(nu.size());
    for (int i = 0; i < nu.size(); i++) {
        nus[i] = abs(pow(nu[i], (1/n + 1)) + pow((0.2/k_4), (1/n)) * nu[i] - pow((0.2/k_4), (1/n)));
    }
    return nus;
}

std::vector<double> reduce_array(std::vector<double> array, double number) {
    for (int i = 0; i < array.size(); i++) {
        array[i] = number - array[i];
    }
    return array;
}

template<typename T>
std::vector<double> logspace(T start_in, T end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(pow(10, start + delta * i));
    }
  linspaced.push_back(pow(10, end)); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}

double get_rms_theta(double e_i, double e_f, double mass, double Z) {
    auto nu = (e_i - e_f) / (e_i - mass); 
    if (nu <= 0.5) {
        auto k_1 = 0.092 * pow(e_i, -1.0/3.0);
        auto k_2 = 0.052 / e_i * pow(Z, -0.25);
        auto k_3 = 0.22 * pow(e_i, -0.92);
        auto rms_theta = ;
        return rms_theta; 
    }
    else if (nu > 0.5) {
        auto k_4 = 0.26 * pow(e_i, -0.91);
        auto m = 0.5;
        auto d = 1.8;
        auto n = 0.81 * pow(e_i, m) / (pow(e_i, m) + d);
        auto rms_theta = k_4 * pow(nu, 1 + n) * pow(1 - nu, -n);
        if (rms_theta < 0.2) {
            return rms_theta; 
        } else 

        auto nu_g = get_nu_g(e_i);
        auto k_5 = k_4 * pow(nu_g, 1 + n) * pow(1 - nu_g, 0.5 - n);
        rms_theta = k_5 * pow(1 - nu, -0.5);
        if (rms_theta >= 0.2) {
            return rms_theta;
        } else {
            // If no case matches
            throw invalid_argument("bremsstrahlung deflection: no matching nu/theta");
        }
    }
}

UnitSphericalVector 
stochastic_deflection::TsaiParametrizationBremsstrahlungGinneken::CalculateStochasticDeflection(
    double e_i, double e_f, std::vector<double> const& rnd) const 
{

    auto Z = 1; // use Z of medium
    auto rms_theta = get_rms_theta(e_i, e_f, mass, Z)
    
    auto lambda = 1 / (rms_theta * rms_theta);
    // Need sqrt because of sampling in theta^2
    auto theta_muon = std::sqrt(SampleFromExponential(rnd[0], lambda));
    
    return UnitSphericalVector {theta_muon, 2 * PI * rnd[1]};
}

