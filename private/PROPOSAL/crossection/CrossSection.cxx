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
    hash_combine(hash, GetParametrization().GetHash(), GetParametrization().GetMultiplier());
    return hash;
}

const std::vector<Components::Component> CrossSection::components_empty_{};

CrossSectionBuilder::CrossSectionBuilder(std::string name) : CrossSection(), name(name){
    // Set zero return functions as default
    dEdx_function = [](double x)->double {
        (void) x; throw std::logic_error("dEdx_function must be set first");
    };
    dE2dx_function = [](double x)->double {
        (void)x; throw std::logic_error("dE2dx_function must be set first");
    };
    dNdx_function = [](double x)->double {
        (void)x; throw std::logic_error("dNdx_function must be set first");
    };
    dNdx_rnd_function = [](double x1, double x2)->double {
        (void)x1; (void)x2; throw std::logic_error("dNdx_rnd_function must be set first");
    };
    StochasticLoss_function = [](double x1, double x2, double x3)->double {
        (void)x1; (void)x2; (void)x3; throw std::logic_error("StochastcLoss_function must be set first");
    };
    CumulativeCrossSection_function = [](double x1, double x2, double x3)->double {
        (void)x1; (void)x2; (void)x3; throw std::logic_error("CumulativeCrossSection_function must be set first");
    };
}

bool CrossSectionBuilder::compare(const CrossSection& cross) const{
    const CrossSectionBuilder* cross_compare = static_cast<const CrossSectionBuilder*>(&cross);
    if(&dEdx_function != &cross_compare->dEdx_function)
        return false;
    if(&dE2dx_function != &cross_compare->dE2dx_function)
        return false;
    if(&dNdx_function != &cross_compare->dNdx_function)
        return false;
    if(&dNdx_rnd_function != &cross_compare->dNdx_rnd_function)
        return false;
    if(&StochasticLoss_function != &cross_compare->StochasticLoss_function)
        return false;
    if(&CumulativeCrossSection_function != &cross_compare->CumulativeCrossSection_function)
        return false;
}

double CrossSectionBuilder::CalculatedEdx(double energy){
    return dEdx_function(energy);
}

double CrossSectionBuilder::CalculatedE2dx(double energy){
    return dE2dx_function(energy);
}

double CrossSectionBuilder::CalculatedNdx(double energy){
    return dNdx_function(energy);
}

double CrossSectionBuilder::CalculatedNdx(double energy, double rnd){
    return dNdx_rnd_function(energy, rnd);
}

double CrossSectionBuilder::CalculateStochasticLoss(double energy, double rnd1, double rnd2){
    return StochasticLoss_function(energy, rnd1, rnd2);
}

double CrossSectionBuilder::CalculateCumulativeCrossSection(double energy, int component, double v){
    return CumulativeCrossSection_function(energy, component, v);
}

void CrossSectionBuilder::SetdEdx_function(Function func){
    dEdx_function = func;
}

void CrossSectionBuilder::SetdE2dx_function(Function func){
    dE2dx_function = func;
}

void CrossSectionBuilder::SetdNdx_function(Function func){
    dNdx_function = func;
}

void CrossSectionBuilder::SetdNdx_rnd_function(std::function<double(double, double)> func){
    dNdx_rnd_function = func;
}

void CrossSectionBuilder::SetStochasticLoss_function(std::function<double(double, double, double)> func){
    StochasticLoss_function = func;
}

void CrossSectionBuilder::SetCumulativeCrossSection_function(std::function<double(double, double, double)> func){
    CumulativeCrossSection_function = func;
}

double CrossSectionBuilder::CalculateStochasticLoss(double energy, double rnd1){
    (void) energy; (void) rnd1;
    return 0;
}

size_t CrossSectionBuilder::GetHash() const{
    std::size_t seed = 0;
    hash_combine(seed, name);

    return seed;
}

Parametrization* CrossSectionBuilder::param;

