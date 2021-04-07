
#include <cmath>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/crosssection/parametrization/Photonuclear.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;
using std::make_shared;

std::vector<double> crosssection::HardComponent::x = { 3, 4, 5, 6, 7, 8, 9 };

double crosssection::Photonuclear::GetLowerEnergyLim(const ParticleDef& p_def) const noexcept {
    return p_def.mass;
}

crosssection::HardComponent::HardComponent(const ParticleDef& particle_def)
{
    hash_combine(hash, std::string("hard_component"));
    const auto& y = particle_def.hard_component_table;

    if (!y.empty()) {
        for (unsigned int i = 0; i < y.size(); i++) {
            interpolant_.emplace_back(x, y.at(i), 4, false, false);
        }
    } else {
        Logging::Get("proposal.parametrization")
            ->error(
                "No HardComponent tables provided for the given particle %s",
                particle_def.name.c_str());
    }
}

double crosssection::HardComponent::CalculateHardComponent(
    double energy, double v)
{
    //TODO: This hard cutoff in v causes kinks in dNdx for higher energies
    //TODO: And the hard cutoff in energy causes a kink in dEdx for E=1e5
    if (energy < 1.0e5 || v < 1.0e-7) {
        return 0;
    }

    double aux, sum, lov, loe;

    sum = 0;
    aux = 1;
    lov = std::log(v) / LOG10;
    loe = std::log(energy) / LOG10 - 3;

    for (unsigned int i = 0; i < interpolant_.size(); i++) {
        if (i > 0) {
            aux *= lov;
        }

        sum += aux * interpolant_[i].InterpolateArray(loe);
    }
    return sum / v;
}

double crosssection::SoftComponent::CalculateHardComponent(double, double)
{
    return 0;
}

crosssection::KinematicLimits crosssection::Photonuclear::GetKinematicLimits(
    const ParticleDef& p_def, const Component& comp, double energy) const
{
    auto kin_lim = KinematicLimits();
    kin_lim.v_min
        = (MPI + MPI * MPI / (2 * comp.GetAverageNucleonWeight())) / energy;
    kin_lim.v_max = 1.f;
    if (p_def.mass < MPI) {
        auto aux = p_def.mass / comp.GetAverageNucleonWeight();
        kin_lim.v_max
            -= comp.GetAverageNucleonWeight() * (1 + aux * aux) / (2 * energy);
    }
    // vMax calculated above is always smaller than 1-m/E
    // in comparison, the following inequality arise
    // (M-m)^2 >= 0
    // limits.vMax = std::min(limits.vMax, 1 - p_def.mass/energy);
    if (kin_lim.v_max < kin_lim.v_min)
        kin_lim.v_max = kin_lim.v_min;
    return kin_lim;
}
