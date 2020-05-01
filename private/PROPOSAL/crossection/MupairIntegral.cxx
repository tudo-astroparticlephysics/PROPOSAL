
#include <functional>

#include "PROPOSAL/crossection/MupairIntegral.h"
#include "PROPOSAL/crossection/parametrization/MupairProduction.h"

using namespace PROPOSAL;

MupairIntegral::MupairIntegral(
    unique_ptr<MupairProduction>&& param, shared_ptr<const EnergyCutSettings> cut)
    : CrossSectionIntegral(forward<unique_ptr<MupairProduction>>(param), cut)
{
}
