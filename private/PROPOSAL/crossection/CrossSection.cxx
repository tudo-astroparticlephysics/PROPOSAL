#include <cassert>

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/CrossSectionIntegral.h"
#include "PROPOSAL/crossection/CrossSectionInterpolant.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

CrossSection::CrossSection(shared_ptr<const EnergyCutSettings> cut)
    : cuts_(cut)
{
}

namespace PROPOSAL {

} // namespace PROPOSAL
