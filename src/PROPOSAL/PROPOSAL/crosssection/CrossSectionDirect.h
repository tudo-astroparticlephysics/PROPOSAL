#pragma once

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/parametrization/ParametrizationDirect.h"

namespace PROPOSAL {
    class CrossSectionDirect : public CrossSectionBase {
    public:
        template <typename Param,
                typename _name = crosssection::ParametrizationName<Param>>
        CrossSectionDirect(Param param_, ParticleDef p,
                           Medium m, std::shared_ptr<const EnergyCutSettings> cut,
                           bool interpol) : param(param_.clone()),
                                            p(p),
                                            m(m),
                                            cut(cut),
                                            interpol(interpol),
                                            param_name(_name::value){};
        double CalculatedEdx(double energy) {
            return param->CalculatedEdx(energy, p, m, cut); };
        double CalculatedE2dx(double energy) {
            return param->CalculatedE2dx(energy, p, m, cut); };
        double CalculatedNdx(double energy) {
            return param->CalculatedNdx(energy, p, m, cut);
        };
        double CalculatedNdx(double energy, size_t hash) {
            return param->CalculatedNdx(energy, hash, p, m, cut);
        };
        double CalculateCumulativeCrosssection(double energy, size_t hash, double v) {
            return param->CalculateCumulativeCrosssection(energy, hash, v, p, m, cut);
        };
        std::vector<std::pair<size_t, double>> CalculatedNdx_PerTarget(
                double energy) {
            return param->CalculatedNdx_PerTarget(energy, p, m, cut);
        };
        double CalculateStochasticLoss(size_t hash, double energy, double rate) {
            return param->CalculateStochasticLoss(hash, energy, rate, p, m, cut);
        };
        double GetLowerEnergyLim() const {
            return param->GetLowerEnergyLim(p, m, cut);
        };
        size_t GetHash() const noexcept {
            return param->GetHash(p, m, cut);
        };
        InteractionType GetInteractionType() const noexcept {
            return param->GetInteractionType();
        };
        std::string GetParametrizationName() const noexcept {
            return param_name;
        };

    private:
        std::unique_ptr<crosssection::ParametrizationDirect> param;
        ParticleDef p;
        Medium m;
        std::shared_ptr<const EnergyCutSettings> cut;
        bool interpol;
        std::string param_name;
    };
}

