
#include "PROPOSAL/propagation_utility/PropagationUtilityFactory.h"
#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include "PROPOSAL/Output.h"

using namespace PROPOSAL;

PropagationUtilityFactory::Definition::Definition()
    : e_cut(500)
    , v_cut(0.05)
    , do_interpolation(true)
    , medium(MediumFactory::Water)
    , density_correction(1.0)
{
}

PropagationUtilityFactory::Definition::~Definition()
{
}

PropagationUtility* PropagationUtilityFactory::CreatePropagationUtility(const ParticleDef& particle_def, const Definition& def)
{
    Medium* med = MediumFactory::Get().CreateMedium(def.medium, def.density_correction);
    EnergyCutSettings cuts(def.e_cut, def.v_cut);

    PropagationUtility* sec;

    if (def.do_interpolation)
    {
        sec = new PropagationUtilityInterpolant(particle_def, *med, cuts, def);
    }
    else
    {
        sec = new PropagationUtilityIntegral(particle_def, *med, cuts, def);
    }

    delete med;

    return sec;
}
