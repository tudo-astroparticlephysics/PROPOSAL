
#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// CrossSection
// ------------------------------------------------------------------------- //

CrossSection::CrossSection(Parametrization& param)
    : parametrization_(param)
    , prob_for_component_(param.GetMedium().GetNumComponents(), 0)
    , sum_of_rates_(0)
    , dedx_integral_(IROMB, IMAXS, IPREC)
    , dndx_integral_(param.GetMedium().GetNumComponents(), Integral(IROMB, IMAXS, IPREC))
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
    :parametrization_(cross_section.parametrization_)
    , prob_for_component_(cross_section.prob_for_component_)
    , sum_of_rates_(cross_section.sum_of_rates_)
    , dedx_integral_(cross_section.dedx_integral_)
    , dndx_integral_(cross_section.dndx_integral_)
    , components_(cross_section.components_)
    , rnd_(cross_section.rnd_)
{
}

CrossSection::~CrossSection()
{
    // for(IntegralVec::iterator it = dndx_integral_.begin(); it != dndx_integral_.end(); ++it)
    // {
    //         delete *it;
    // }
}
