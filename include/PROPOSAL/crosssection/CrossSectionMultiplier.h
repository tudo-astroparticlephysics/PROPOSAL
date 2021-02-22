#pragma once

namespace PROPOSAL {
    template <typename T>
    class CrossSectionMultiplier : public T {
    public:
        CrossSectionMultiplier(std::shared_ptr<CrossSectionBase> cross, double multiplier) : cross_(cross), multiplier_(multiplier) {}

        double CalculatedEdx(double energy) override {
            return multiplier_ * cross_->CalculatedEdx(energy);
        }

        double CalculatedE2dx(double energy) override {
            return multiplier_ * cross_->CalculatedE2dx(energy);
        }

        double CalculatedNdx(double energy, std::shared_ptr<const Component> comp = nullptr) override {
            return multiplier_ * cross_->CalculatedNdx(energy, comp);
        };

        double CalculateStochasticLoss(std::shared_ptr<const Component> const& comp, double energy, double rate) override {
            return cross_->CalculateStochasticLoss(comp, energy, rate/multiplier_);
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

        std::vector<std::shared_ptr<const Component>> GetTargets() const noexcept override {
            return cross_->GetTargets();
        }
    private:
        std::shared_ptr<CrossSectionBase> cross_;
        double multiplier_;

    };

    template <typename T>
    auto make_crosssection_multiplier(std::shared_ptr<T> cross, double multiplier) {
        return std::unique_ptr<T>(new CrossSectionMultiplier<T>(cross, multiplier));
    }

}
