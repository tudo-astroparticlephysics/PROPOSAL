
#include <boost/algorithm/string.hpp> // case insensitive string compare for configuration file

#include "PROPOSAL/crossection/IonizIntegral.h"
#include "PROPOSAL/crossection/IonizInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"
#include "PROPOSAL/crossection/factories/IonizationFactory.h"

#include "PROPOSAL/Output.h"

using namespace PROPOSAL;

IonizationFactory::IonizationFactory()
{
}

IonizationFactory::~IonizationFactory()
{
}

// ------------------------------------------------------------------------- //
// PhotoQ2
// ------------------------------------------------------------------------- //

CrossSection* IonizationFactory::CreateIonizationIntegral(const ParticleDef& particle_def,
                                                          const Medium& medium,
                                                          const EnergyCutSettings& cuts,
                                                          double multiplier) const
{
    return new IonizIntegral(Ionization(particle_def, medium, cuts, multiplier));
}

CrossSection* IonizationFactory::CreateIonizationInterpolant(const ParticleDef& particle_def,
                                                             const Medium& medium,
                                                             const EnergyCutSettings& cuts,
                                                             double multiplier,
                                                             InterpolationDef def) const
{
    return new IonizInterpolant(Ionization(particle_def, medium, cuts, multiplier), def);
}


// ------------------------------------------------------------------------- //
// Most general creator
// ------------------------------------------------------------------------- //

CrossSection* IonizationFactory::CreateIonization(const ParticleDef& particle_def,
                                                  const Medium& medium,
                                                  const EnergyCutSettings& cuts,
                                                  double multiplier,
                                                  bool interpolate,
                                                  InterpolationDef def) const
{
    if (interpolate)
    {
       return CreateIonizationInterpolant(particle_def, medium, cuts, multiplier, def);
    }
    else
    {
       return CreateIonizationIntegral(particle_def, medium, cuts, multiplier);
    }
}

CrossSection* IonizationFactory::CreateIonization(const ParticleDef& particle_def,
                                                  const Medium& medium,
                                                  const EnergyCutSettings& cuts,
                                                  const Definition& def,
                                                  bool interpolate,
                                                  InterpolationDef interpolation_def) const
{
    if (interpolate)
    {
       return CreateIonizationInterpolant(particle_def, medium, cuts, def.multiplier, interpolation_def);
    }
    else
    {
       return CreateIonizationIntegral(particle_def, medium, cuts, def.multiplier);
    }
}
