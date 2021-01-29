
#include <cmath>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
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
