#include "PROPOSAL/propagation_utility/PropagationUtilityStochastic.h"

using namespace PROPOSAL;

bool PropagationUtilityStochastic::Collection::operator==(const Collection& lhs)
{
    if (interaction_calc != lhs.interaction_calc)
        return false;
    return true;
}

PropagationUtilityStochastic::PropagationUtilityStochastic(
        PropagationUtilityStochastic::Collection const& collect)
        : collection(collect)
{
    if (collect.interaction_calc == nullptr) {
        throw std::invalid_argument("Interaction calculator need to be defined.");
    }
}

Interaction::Loss PropagationUtilityStochastic::EnergyStochasticloss(double energy,
                                                                     double rnd) const
                                                                     {
    auto rates = collection.interaction_calc->Rates(energy);
    auto loss = collection.interaction_calc->SampleLoss(energy, rates, rnd);

    return loss;
}

double PropagationUtilityStochastic::EnergyDecay(double, std::function<double()>,
        double) const
{
    return 0; // no decay, e.g. particle is stable
}

std::pair<double, double> PropagationUtilityStochastic::EnergyDistanceStochasticInteraction(double E_i,
                                                       std::function<double()> rnd) const
{
    auto rndi = -std::log(rnd());
    auto grammage_to_interaction = rndi * collection.interaction_calc->MeanFreePath(E_i);
    return std::make_pair(E_i, grammage_to_interaction); // no continuous losses
}

double PropagationUtilityStochastic::EnergyRandomize(double, double final_energy,
                                                     std::function<double()>) const
{
    return final_energy; // no randomization
}

double PropagationUtilityStochastic::EnergyDistance(double initial_energy,
                                                    double) const
{
    return initial_energy; // no continuous losses
}

double PropagationUtilityStochastic::TimeElapsed(double initial_energy, double final_energy, double distance, double density) const
{
    return collection.time_calc->TimeElapsed(
            initial_energy, final_energy, distance, density);
}

std::tuple<Cartesian3D, Cartesian3D> PropagationUtilityStochastic::DirectionsScatter(
        double, double, double, const Vector3D& direction, std::function<double()>) const
{
    auto dir = Cartesian3D(direction.GetCartesianCoordinates());
    return std::make_tuple(dir, dir); // no scattering
}

Cartesian3D PropagationUtilityStochastic::DirectionDeflect(
        InteractionType, double, double, const Vector3D& direction,
        std::function<double()>) const
{
    return direction; // no scattering as of now
}

double PropagationUtilityStochastic::LengthContinuous(
        double initial_energy, double final_energy) const
{
    assert(initial_energy == final_energy);
    return initial_energy; // no continuous losses
}

double PropagationUtilityStochastic::GetLowerPropagationLim() const {
    return 0.; // TODO : better fitting value? -> replace with LowerLim
}