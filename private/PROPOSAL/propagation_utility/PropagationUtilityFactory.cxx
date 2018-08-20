
#include "PROPOSAL/propagation_utility/PropagationUtilityFactory.h"
#include "PROPOSAL/medium/Medium.h"

using namespace PROPOSAL;

// UtilityFactory::Definition::Definition()
//     : e_cut(500)
//     , v_cut(0.05)
//     , medium(MediumFactory::Water)
//     , density_correction(1.0)
// {
// }
//
// UtilityFactory::Definition::~Definition() {}
//
// Utility* UtilityFactory::CreateUtility(const ParticleDef& particle_def, const Definition& def)
// {
//     Medium* med = MediumFactory::Get().CreateMedium(def.medium, def.density_correction);
//     EnergyCutSettings cuts(def.e_cut, def.v_cut);
//
//     Utility* sec = new Utility(particle_def, *med, cuts, def);
//
//     delete med;
//
//     return sec;
// }
