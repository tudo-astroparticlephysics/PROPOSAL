#pragma once

#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/propagation_utility/TimeBuilder.h"

namespace PROPOSAL {
    class PropagationUtilityStochastic : public PropagationUtility {
    public:
        struct Collection {

            bool operator==(const Collection& lhs);

            // obligatory pointers
            std::shared_ptr<Interaction> interaction_calc; // replace with more fitting utility?
            std::shared_ptr<ApproximateTimeBuilder> time_calc;

        };

        PropagationUtilityStochastic(Collection const& collection);

        Interaction::Loss EnergyStochasticloss(double, double) const override;
        double EnergyDecay(double, std::function<double()>, double) const override;
        std::pair<double, double> EnergyDistanceStochasticInteraction(double, std::function<double()>) const override;
        double EnergyRandomize(double, double, std::function<double()>) const override;
        double EnergyDistance(double, double) const override;
        double LengthContinuous(double, double) const override;
        double TimeElapsed(double, double, double, double) const override;

        // TODO: return value doesn't tell what it include. Maybe it would be better
        // to give a tuple of two directions back. One is the mean over the
        // displacement and the other is the actual direction. With a get method
        // there could be a possible access with the position of the object stored
        // in an enum.

        std::tuple<Cartesian3D, Cartesian3D> DirectionsScatter(
                double, double, double, const Vector3D&, std::function<double()>) const override;
        Cartesian3D DirectionDeflect(InteractionType, double, double,
                                     const Vector3D&, std::function<double()>) const override;

        double GetLowerPropagationLim() const override;


        Collection collection;
    };
} // namespace PROPOSAL
