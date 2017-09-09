
#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// CrossSection
// ------------------------------------------------------------------------- //

CrossSection::CrossSection(const Parametrization& param)
    : parametrization_(param.clone())
    , prob_for_component_(param.GetMedium().GetNumComponents(), 0)
    , sum_of_rates_(0)
    , components_(param.GetMedium().GetComponents())
    , rnd_(0)
{

    // int number_of_components = param.GetMedium().GetNumComponents();
    //
    // dndx_integral_.resize(number_of_components);
    //
    // for(IntegralVec::iterator it = dndx_integral_.begin(); it != dndx_integral_.end(); ++it)
    // {
    //         *it =  new Integral(IROMB, IMAXS, IPREC);
    // }

    // prob_for_component_.resize(number_of_components);
}

CrossSection::CrossSection(const CrossSection& cross_section)
    :parametrization_(cross_section.parametrization_->clone())
    , prob_for_component_(cross_section.prob_for_component_)
    , sum_of_rates_(cross_section.sum_of_rates_)
    , components_(parametrization_->GetMedium().GetComponents())
    , rnd_(cross_section.rnd_)
{
}

CrossSection::~CrossSection()
{
    delete parametrization_;
    // for(IntegralVec::iterator it = dndx_integral_.begin(); it != dndx_integral_.end(); ++it)
    // {
    //         delete *it;
    // }
}
