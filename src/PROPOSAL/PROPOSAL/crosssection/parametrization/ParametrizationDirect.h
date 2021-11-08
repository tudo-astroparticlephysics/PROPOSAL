#pragma once

#include <memory>
#include "PROPOSAL/particle/Particle.h"

namespace PROPOSAL{
    class EnergyCutSettings;
    struct ParticleDef;
    class Medium;
}
namespace PROPOSAL {
    namespace crosssection {
        class ParametrizationDirect {
        protected:
            using cut_ptr = std::shared_ptr<const EnergyCutSettings>;
            size_t hash;
        public:
            ParametrizationDirect() : hash(0) {};
            virtual ~ParametrizationDirect() = default;

            virtual std::unique_ptr<ParametrizationDirect> clone() const = 0;

            virtual double CalculatedEdx(
                    double, const ParticleDef&, const Medium&, cut_ptr) = 0;
            virtual double CalculatedE2dx(
                    double, const ParticleDef&, const Medium&, cut_ptr) = 0;
            virtual double CalculatedNdx(
                    double, const ParticleDef&, const Medium&, cut_ptr) = 0;
            virtual double CalculatedNdx(
                    double, size_t, const ParticleDef&, const Medium&, cut_ptr) = 0;
            virtual double CalculateCumulativeCrosssection(
                    double, size_t, double, const ParticleDef&, const Medium&, cut_ptr) = 0;
            virtual std::vector<std::pair<size_t, double>> CalculatedNdx_PerTarget(
                    double, const ParticleDef&, const Medium&, cut_ptr) = 0;
            virtual double CalculateStochasticLoss(
                    size_t, double, double, const ParticleDef&, const Medium&, cut_ptr) = 0;
            virtual double GetLowerEnergyLim(
                    const ParticleDef&, const Medium&, cut_ptr) const = 0;
            virtual size_t GetHash(
                    const ParticleDef&, const Medium&, cut_ptr) const noexcept = 0;
            virtual InteractionType GetInteractionType() const noexcept = 0;
        };
    }
}