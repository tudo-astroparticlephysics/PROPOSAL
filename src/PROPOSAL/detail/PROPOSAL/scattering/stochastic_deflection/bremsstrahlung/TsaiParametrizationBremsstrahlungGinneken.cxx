#include "PROPOSAL/scattering/stochastic_deflection/bremsstrahlung/TsaiParametrizationBremsstrahlungGinneken.h"
#include "PROPOSAL/Constants.h"
#include <cmath> 
using namespace PROPOSAL;

UnitSphericalVector 
stochastic_deflection::TsaiParametrizationBremsstrahlungGinneken::CalculateStochasticDeflection(
    double e_i, double e_f, std::vector<double> const& rnd) const 
{

    return UnitSphericalVector {theta_muon, 2 * PI * rnd[1]};
}