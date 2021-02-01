#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDX.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;

CrossSectionDEDX::CrossSectionDEDX(crosssection::Parametrization<Medium> const&,
    ParticleDef const&, Medium const&, EnergyCutSettings const&) {};

CrossSectionDEDX::CrossSectionDEDX(
    crosssection::Parametrization<Component> const&, ParticleDef const&,
    Component const&, EnergyCutSettings const&) {};
