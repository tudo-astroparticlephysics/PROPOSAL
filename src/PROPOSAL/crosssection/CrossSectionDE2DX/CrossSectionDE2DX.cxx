#include "PROPOSAL/crosssection/CrossSectionDE2DX/CrossSectionDE2DX.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;

CrossSectionDE2DX::CrossSectionDE2DX(
    crosssection::Parametrization<Component> const& param, ParticleDef const& p,
    Component const& comp, EnergyCutSettings const& cut)
{
    hash_combine(
        hash, param.GetHash(), p.GetHash(), comp.GetHash(), cut.GetHash());
}

CrossSectionDE2DX::CrossSectionDE2DX(
    crosssection::Parametrization<Medium> const& param, ParticleDef const& p,
    Medium const& medium, EnergyCutSettings const& cut)
{
    hash_combine(
        hash, param.GetHash(), p.GetHash(), medium.GetHash(), cut.GetHash());
}
