#pragma once
#include "PROPOSAL/crosssection/CrossSection.h"
namespace PROPOSAL {
    struct CrossSectionMultiplier : public CrossSectionBase {
        CrossSectionMultiplier(std::shared_ptr<CrossSectionBase> cross, double multiplier) : CrossSectionBase(), cross_(cross), multiplier_(multiplier) {}

        double CalculatedEdx(double energy) override {
            return multiplier_ * cross_->CalculatedEdx(energy);
        }

        double CalculatedE2dx(double energy) override {
            return multiplier_ * cross_->CalculatedE2dx(energy);
        }

        double CalculatedNdx(double energy) override {
            return multiplier_ * cross_->CalculatedNdx(energy);
        }

        double CalculatedNdx(double energy, size_t comp_hash) override {
            return multiplier_ * cross_->CalculatedNdx(energy, comp_hash);
        };

        double CalculateCumulativeCrosssection(double energy, size_t comp_hash, double v) override {
            return multiplier_ * cross_->CalculateCumulativeCrosssection(energy, comp_hash, v);
        }

        std::vector<std::pair<size_t, double>> CalculatedNdx_PerTarget(double energy) override {
            auto rates = cross_->CalculatedNdx_PerTarget(energy);
            for (auto &val : rates)
                val.second *= multiplier_;
            return rates;
        }

        double CalculateStochasticLoss(size_t comp_hash, double energy, double rate) override {
            return cross_->CalculateStochasticLoss(comp_hash, energy, rate/multiplier_);
        };

        double GetLowerEnergyLim() const override {
            return cross_->GetLowerEnergyLim();
        }

        size_t GetHash() const noexcept override {
            auto hash = cross_->GetHash();
            hash_combine(hash, multiplier_);
            return hash;
        }

        InteractionType GetInteractionType() const noexcept override {
            return cross_->GetInteractionType();
        }

        std::string GetParametrizationName() const noexcept override {
            return cross_->GetParametrizationName();
        }

    private:
        std::shared_ptr<CrossSectionBase> cross_;
        double multiplier_;
    };

}

namespace PROPOSAL {
    inline auto make_crosssection_multiplier(std::shared_ptr<CrossSectionBase> cross, double multiplier) {
        return std::unique_ptr<CrossSectionBase>(new CrossSectionMultiplier(cross, multiplier));
    }
}