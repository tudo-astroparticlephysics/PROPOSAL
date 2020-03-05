
#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// CrossSection
// ------------------------------------------------------------------------- //

CrossSection::CrossSection(const Parametrization& param, std::shared_ptr<const EnergyCutSettings> cuts)
    : parametrization_(param.clone())
    , prob_for_component_(param.GetMedium()->GetNumComponents(), 0)
    , sum_of_rates_(0)
    , components_(parametrization_->GetMedium()->GetComponents())
    , rnd_(0)
    , cuts_(cuts)
{
}

CrossSection::CrossSection(const CrossSection& cross_section)
    : parametrization_(cross_section.parametrization_->clone())
    , prob_for_component_(cross_section.prob_for_component_)
    , sum_of_rates_(cross_section.sum_of_rates_)
    , components_(parametrization_->GetMedium()->GetComponents())
    , rnd_(cross_section.rnd_)
    , cuts_(cross_section.cuts_)
{
}

CrossSection::~CrossSection()
{
    delete parametrization_;
}

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
    else if(cuts_ != cross_section.cuts_)
        return false;
    else
        return this->compare(cross_section);
}

bool CrossSection::operator!=(const CrossSection& cross_section) const
{
    return !(*this == cross_section);
}

// ------------------------------------------------------------------------- //
std::ostream& PROPOSAL::operator<<(std::ostream& os, CrossSection const& cross)
{
    std::stringstream ss;
    ss << " CrossSection (" << &cross << ") ";
    os << Helper::Centered(80, ss.str()) << '\n';

    os << *cross.parametrization_ << '\n';

    os << Helper::Centered(80, "");
    return os;
}

double CrossSection::GetEnergyCut(double energy){
    Parametrization::KinematicLimits physical_limits = parametrization_->GetKinematicLimits(energy);

    if(cuts_ == nullptr) throw std::logic_error("CrossSection is only stochastic, therefore there is no energy cut defined.");

    double energy_cut = std::min(physical_limits.vMax, cuts_->GetCut(energy));

    if (energy_cut < physical_limits.vMin)
    {
        energy_cut = physical_limits.vMin;
    }

    return energy_cut;
}


std::pair<double, double> CrossSection::StochasticDeflection(double energy, double energy_loss){
    // per default the particle is not deflected
    (void) energy;
    (void) energy_loss;
    return std::make_pair(1, 0);
}
