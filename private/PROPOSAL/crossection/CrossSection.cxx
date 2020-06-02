#include <cassert>

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// CrossSection
// ------------------------------------------------------------------------- //

CrossSection::CrossSection(
    const Parametrization& param, std::shared_ptr<const EnergyCutSettings> cuts)
    : parametrization_(param.clone())
    , prob_for_component_(param.GetMedium()->GetNumComponents(), 0)
    , sum_of_rates_(0)
    , components_(parametrization_->GetMedium()->GetComponents())
    , rnd_(0)
    , cuts_(cuts)
    , type_id_(GetParametrization().GetInteractionType())
{
}

CrossSection::CrossSection(const CrossSection& cross_section)
    : parametrization_(cross_section.parametrization_->clone())
    , prob_for_component_(cross_section.prob_for_component_)
    , sum_of_rates_(cross_section.sum_of_rates_)
    , components_(parametrization_->GetMedium()->GetComponents())
    , rnd_(cross_section.rnd_)
    , cuts_(cross_section.cuts_)
    , type_id_(cross_section.type_id_)
{
}

CrossSection::CrossSection()
        : parametrization_(nullptr)
        , prob_for_component_(0, 0)
        , sum_of_rates_(0)
        , components_(components_empty_)
        , rnd_(0)
        , cuts_(nullptr)
{
}

CrossSection::~CrossSection() { delete parametrization_; }

bool CrossSection::operator==(const CrossSection& cross_section) const
{

    if (typeid(*this) != typeid(cross_section))
        return false;
    else if (*parametrization_ != *cross_section.parametrization_)
        return false;
    else if (prob_for_component_ != cross_section.prob_for_component_)
        return false;
    else if (sum_of_rates_ != cross_section.sum_of_rates_)
        return false;
    else if (rnd_ != cross_section.rnd_)
        return false;
    else if (cuts_ != cross_section.cuts_)
        return false;
    else if (type_id_ != cross_section.type_id_)
        return false;
    else
        return this->compare(cross_section);
}

bool CrossSection::operator!=(const CrossSection& cross_section) const
{
    return !(*this == cross_section);
}

namespace PROPOSAL {
std::ostream& operator<<(std::ostream& os, CrossSection const& cross)
{
    std::stringstream ss;
    ss << " CrossSection (" << &cross << ") ";
    os << Helper::Centered(80, ss.str()) << '\n';

    os << *cross.parametrization_ << '\n';

    os << Helper::Centered(80, "");
    return os;
}
} // namespace PROPOSAL

double CrossSection::GetEnergyCut(double energy)
{
    auto physical_limits = parametrization_->GetKinematicLimits(energy);
    auto energy_cut = physical_limits.vMin;

    assert(physical_limits.vMin <= physical_limits.vMax);

    if (cuts_)
        return std::min(std::max(physical_limits.vMin, cuts_->GetCut(energy)),
            physical_limits.vMax);
    else
        return physical_limits.vMin;
}

std::pair<double, double> CrossSection::StochasticDeflection(
    double energy, double energy_loss)
{
    // per default the particle is not deflected
    (void)energy;
    (void)energy_loss;
    return std::make_pair(1, 0);
}

size_t CrossSection::GetHash() const{
    std::size_t hash = 0;
    hash_combine(hash, GetParametrization().GetHash(), GetParametrization().GetMultiplier(), cuts_->GetHash());
    return hash;
}

const std::vector<Components::Component> CrossSection::components_empty_{};
