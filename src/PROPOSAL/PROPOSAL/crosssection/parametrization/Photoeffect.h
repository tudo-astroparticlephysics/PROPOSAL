#pragma once

#include "PROPOSAL/crosssection/parametrization/ParametrizationDirect.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"

namespace PROPOSAL {
    class Component;

    namespace crosssection {
        class Photoeffect : public ParametrizationDirect {
            virtual double PhotoeffectKshellCrossSection(double, const Component&) = 0;
        protected:
            double GetCutOff(const Component& comp) const;
        public:
            Photoeffect() = default;

            virtual double PhotonAtomCrossSection(double, const Component&);

            double CalculatedNdx(double, const ParticleDef&, const Medium&, cut_ptr) override;
            double CalculatedNdx(double, size_t, const ParticleDef&, const Medium&, cut_ptr) override;
            std::vector<std::pair<size_t, double>> CalculatedNdx_PerTarget(
                    double, const ParticleDef&, const Medium&, cut_ptr) override;

            // no continuous losses
            double CalculatedEdx(double, const ParticleDef&, const Medium&, cut_ptr) override { return 0.; };
            double CalculatedE2dx(double, const ParticleDef&, const Medium&, cut_ptr) override { return 0.; };

            double CalculateCumulativeCrosssection(double, size_t, double, const ParticleDef&, const Medium&, cut_ptr) override {
                throw std::logic_error("No cumulative crosssection defined for Photoproduction.");
            };

            double CalculateStochasticLoss(size_t, double, double, const ParticleDef&, const Medium&, cut_ptr) override {
                return 1.; // all energy is always lost, i.e. v=1
            };

            double GetLowerEnergyLim(const ParticleDef&, const Medium&, cut_ptr) const override;

            size_t GetHash(const ParticleDef&, const Medium& m, cut_ptr) const noexcept override;

            InteractionType GetInteractionType() const noexcept override;
        };

        class PhotoeffectSauter : public Photoeffect {
        public:
            PhotoeffectSauter();
            std::unique_ptr<ParametrizationDirect> clone() const final;
            double PhotoeffectKshellCrossSection(double, const Component&) override;
        };

        template <> struct ParametrizationName<PhotoeffectSauter> {
            static constexpr auto value = "Sauter";
        };
    }
}
