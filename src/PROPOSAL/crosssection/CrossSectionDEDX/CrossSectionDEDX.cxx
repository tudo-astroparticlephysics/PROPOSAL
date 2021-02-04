#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDX.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;


CrossSectionDEDX::CrossSectionDEDX(crosssection::Parametrization<Medium> const& param,
    ParticleDef const& p, Medium const& m, EnergyCutSettings const& cut) {
    hash_combine(hash, param.GetHash(), p.GetHash(), m.GetHash(), cut.GetHash());
};

CrossSectionDEDX::CrossSectionDEDX(
    crosssection::Parametrization<Component> const& param, ParticleDef const& p,
    Component const& c, EnergyCutSettings const& cut) {
    hash_combine(hash, param.GetHash(), p.GetHash(), c.GetHash(), cut.GetHash());
};
