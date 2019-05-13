
#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// CrossSection
// ------------------------------------------------------------------------- //

CrossSection::CrossSection(const DynamicData::Type& type, const Parametrization& param)
    : type_id_(type)
    , parametrization_(param.clone())
    , prob_for_component_(param.GetMedium().GetNumComponents(), 0)
    , sum_of_rates_(0)
    , components_(parametrization_->GetMedium().GetComponents())
    , rnd_(0)
{
}

CrossSection::CrossSection(const CrossSection& cross_section)
    : type_id_(cross_section.type_id_)
    , parametrization_(cross_section.parametrization_->clone())
    , prob_for_component_(cross_section.prob_for_component_)
    , sum_of_rates_(cross_section.sum_of_rates_)
    , components_(parametrization_->GetMedium().GetComponents())
    , rnd_(cross_section.rnd_)
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
    if (type_id_ != cross_section.type_id_)
        return false;
    else if (*parametrization_ != *cross_section.parametrization_)
        return false;
    else if (prob_for_component_ != cross_section.prob_for_component_)
        return false;
    else if (sum_of_rates_ != cross_section.sum_of_rates_)
        return false;
    else if (rnd_ != cross_section.rnd_)
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
    std::string name;
    switch (cross.type_id_)
    {
        case DynamicData::Brems:
            name = "Bremsstrahlung";
            break;
        case DynamicData::Epair:
            name = "EPairProduction";
            break;
        case DynamicData::DeltaE:
            name = "Ionization";
            break;
        case DynamicData::NuclInt:
            name = "PhotoNuclear";
            break;
        case DynamicData::MuPair:
            name = "MuPairProduction";
            break;
        case DynamicData::WeakInt:
            name = "WeakInteraction";
            break;
        default:
            break;
    }
    std::stringstream ss;
    ss << " CrossSection (" << &cross << ") ";
    os << Helper::Centered(80, ss.str()) << '\n';

    os << "type: " << name << '\n';
    os << *cross.parametrization_ << '\n';

    os << Helper::Centered(80, "");
    return os;
}
