#pragma once

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crosssection/parametrization/Annihilation.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/secondaries/RegisteredInDefault.h"
#include "PROPOSAL/secondaries/annihilation/Annihilation.h"

#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXBuilder.h"
/* #include <unordered_map> */

using PROPOSAL::Components::Component;
/* using std::unordered_map; */
using std::array;
using std::vector;

namespace PROPOSAL {
namespace secondaries {
    class SingleDifferentialAnnihilation : public secondaries::Annihilation {
        Medium m;
        dndx_map_t dndx;

    public:
        static constexpr int n_rnd = 2;

        template <typename Param>
        SingleDifferentialAnnihilation(
            Param&& param, const ParticleDef& p, const Medium& medium, bool interpol)
            : m(medium),
              dndx(build_cross_section_dndx(std::forward<Param>(param), p, m,
                      std::make_shared<EnergyCutSettings>(0.f, 1.f, false), interpol))
        {
        }

        double CalculateRho(double, double, const Component&) final;
        tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double) final;
        tuple<double, double> CalculateEnergy(double, double) final;

        size_t RequiredRandomNumbers() const noexcept { return n_rnd; }
        vector<Loss::secondary_t> CalculateSecondaries(
            double, Loss::secondary_t, const Component&, vector<double>);
    };
} // namespace secondaries
} // namespace PROPOSAL
