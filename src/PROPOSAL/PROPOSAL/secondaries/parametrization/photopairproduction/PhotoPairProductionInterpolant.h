#pragma once

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/secondaries/parametrization/photopairproduction/PhotoPairProduction.h"

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXBuilder.h"
#include "PROPOSAL/crosssection/parametrization/PhotoPairProduction.h"

#include <cmath>
#include <tuple>

using std::make_tuple;

namespace PROPOSAL {
namespace secondaries {
    template <class Param>
    class PhotoPairProductionInterpolant : public PhotoPairProduction {

        using dndx_ptr = std::unique_ptr<CrossSectionDNDX>;

        Medium medium;
        std::unique_ptr<
            std::unordered_map<size_t, std::tuple<double, dndx_ptr>>>
            dndx;

    public:
        static constexpr int n_rnd = 5;

        PhotoPairProductionInterpolant() = default;
        PhotoPairProductionInterpolant(ParticleDef p, Medium m)
            : medium(m)
            , dndx(detail::build_dndx(std::true_type {}, true,
                                      Param(), p, medium, nullptr))
        {}

        double CalculateRho(double energy, double rnd,
                            const Component& comp) override {
            if (!dndx)
                throw std::logic_error("dndx Interpolant for PhotoPairProduction not defined.");
            for (auto& it : *dndx) {
                if (comp.GetHash() == it.first) {
                    auto& calc = *std::get<1>(it.second);
                    auto lim = calc.GetIntegrationLimits(energy);
                    auto rate = rnd * calc.Calculate(energy, lim.max);
                    auto rho = calc.GetUpperLimit(energy, rate);
                    return rho;
                }
            }
            std::ostringstream s;
            s << "Component (" << comp.GetName()
            << ") can not be found in the precalculated photopairproduction "
               "tables.";
            throw std::out_of_range(s.str());
        };

        // TODO: two random numbers are always discarded here, this might be improved
        std::tuple<Cartesian3D, Cartesian3D> CalculateDirections(
                const Vector3D& dir, double energy, double rho,
                const Component& comp, double rnd1, double, double) override {
            auto k = energy / ME;

            auto cosphi0 = std::cos(1. / k);
            auto cosphi1 = std::cos(1. / k);

            auto theta0 = rnd1 * 2. * PI;
            auto theta1 = std::fmod(theta0 + PI, 2. * PI);
            if (cosphi0 == -1.)
                cosphi0 *= (-1);
            if (cosphi1 == -1.)
                cosphi1 *= (-1);
            auto dir_0 = Cartesian3D(dir);
            dir_0.deflect(cosphi0, theta0);
            auto dir_1 = Cartesian3D(dir);
            dir_1.deflect(cosphi1, theta1);
            return std::make_tuple(dir_0, dir_1);
        };

        std::tuple<double, double> CalculateEnergy(double energy, double rho,
                                                   double rnd) override {
            if (rnd > 0.5)
                return make_tuple(energy * (1 - rho), energy * rho);
            return make_tuple(energy * rho, energy * (1 - rho));
        };

        size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }

        std::vector<ParticleState> CalculateSecondaries(
                StochasticLoss loss, const Component& comp,
                std::vector<double>& rnd) final {
            auto rho = CalculateRho(loss.parent_particle_energy, rnd[0], comp);
            auto secondary_energies = CalculateEnergy(loss.parent_particle_energy, rho, rnd[1]);
            auto secondary_dir = CalculateDirections(
                    loss.direction, loss.parent_particle_energy, rho, comp,
                    rnd[2], rnd[3], rnd[4]);

            auto sec = std::vector<ParticleState>();
            sec.emplace_back(ParticleType::EMinus, loss.position,
                             std::get<0>(secondary_dir),
                             std::get<0>(secondary_energies), loss.time, 0.);
            sec.emplace_back(ParticleType::EPlus, loss.position,
                             std::get<1>(secondary_dir),
                             std::get<1>(secondary_energies), loss.time, 0.);
            return sec;
        };
    };
} // namespace secondaries
} // namespace PROPOSAL
