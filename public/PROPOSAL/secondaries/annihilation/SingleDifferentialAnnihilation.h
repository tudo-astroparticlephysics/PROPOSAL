
#pragma once

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crossection/parametrization/Annihilation.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/secondaries/DefaultFactory.h"
#include "PROPOSAL/secondaries/annihilation/Annihilation.h"

#include "PROPOSAL/crossection/CrossSectionDNDX/CrossSectionDNDXBuilder.h"
/* #include <unordered_map> */

using PROPOSAL::Components::Component;
/* using std::unordered_map; */
using std::array;
using std::vector;

namespace PROPOSAL {
namespace secondaries {
    class SingleDifferentialAnnihilation : public secondaries::Annihilation {
        dndx_map_t dndx;

    public:
        static constexpr int n_rnd = 2;

        template <typename Param>
        SingleDifferentialAnnihilation(
            Param&& param, const ParticleDef& p, const Medium& m, bool interpol)
            : dndx(build_cross_section_dndx(std::forward<Param>(param), p, m,
                  std::make_shared<EnergyCutSettings>(
                      0, std::numeric_limits<double>::infinity(), false),
                  interpol))
        {
        }

        double CalculateRho(double, double, const Component&) final;
        tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double) final;
        tuple<double, double> CalculateEnergy(double, double) final;

        size_t RequiredRandomNumbers() { return n_rnd; }
        vector<Loss::secondary_t> CalculateSecondaries(
            double, Loss::secondary_t, const Component&, vector<double>);
    };
} // namespace secondaries
} // namespace PROPOSAL
