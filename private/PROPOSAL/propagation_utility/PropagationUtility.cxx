#include <boost/bind.hpp>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/propagation_utility/PropagationUtility.h"

#include "PROPOSAL/crossection/CrossSection.h"

#include "PROPOSAL/crossection/BremsIntegral.h"
#include "PROPOSAL/crossection/PhotoIntegral.h"
#include "PROPOSAL/crossection/IonizIntegral.h"
#include "PROPOSAL/crossection/EpairIntegral.h"

#include "PROPOSAL/crossection/BremsInterpolant.h"
#include "PROPOSAL/crossection/PhotoInterpolant.h"
#include "PROPOSAL/crossection/IonizInterpolant.h"
#include "PROPOSAL/crossection/EpairInterpolant.h"

#include "PROPOSAL/crossection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/crossection/parametrization/Photonuclear.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"
#include "PROPOSAL/crossection/parametrization/EpairProduction.h"


using namespace std;
using namespace PROPOSAL;

/******************************************************************************
*                            Propagation utility                              *
******************************************************************************/

Utility::Definition::Definition()
    // : do_interpolation(true)
    // , interpolation_def()
    : brems_def()
    , photo_def()
    , epair_def()
    , ioniz_def()
{
}

Utility::Definition::~Definition()
{
}

// ------------------------------------------------------------------------- //
// Constructors
// ------------------------------------------------------------------------- //

// Standard constructor
// Utility::Utility(const ParticleDef& particle_def)
//       //TODO(mario): init different Fri 2017/09/01
//     : particle_def_(particle_def)
//     , medium_(new Water())
//     , cut_settings_()
//     , crosssections_()
// {
// }

Utility::Utility(const ParticleDef& particle_def,
                 const Medium& medium,
                 const EnergyCutSettings& cut_settings,
                 Definition utility_def)
    : particle_def_(particle_def)
    , medium_(medium.clone())
    , cut_settings_(cut_settings)
    , crosssections_()
    , utility_def_(utility_def)
{
    crosssections_.push_back(BremsstrahlungFactory::Get().CreateBremsstrahlung(particle_def_,
                                                                               *medium_,
                                                                               cut_settings_,
                                                                               utility_def.brems_def,
                                                                               false));

    crosssections_.push_back(PhotonuclearFactory::Get().CreatePhotonuclear(particle_def_,
                                                                           *medium_,
                                                                           cut_settings_,
                                                                           utility_def.photo_def,
                                                                           false));

    crosssections_.push_back(EpairProductionFactory::Get().CreateEpairProduction(particle_def_,
                                                                                 *medium_,
                                                                                 cut_settings_,
                                                                                 utility_def.epair_def));

    crosssections_.push_back(IonizationFactory::Get().CreateIonization(particle_def_,
                                                                       *medium_,
                                                                       cut_settings_,
                                                                       utility_def.ioniz_def,
                                                                       false));
}

Utility::Utility(const ParticleDef& particle_def,
                 const Medium& medium,
                 const EnergyCutSettings& cut_settings,
                 Definition utility_def,
                 const InterpolationDef& interpolation_def)
    : particle_def_(particle_def)
    , medium_(medium.clone())
    , cut_settings_(cut_settings)
    , crosssections_()
    , utility_def_(utility_def)
{
    crosssections_.push_back(BremsstrahlungFactory::Get().CreateBremsstrahlung(particle_def_,
                                                                               *medium_,
                                                                               cut_settings_,
                                                                               utility_def.brems_def,
                                                                               true,
                                                                               interpolation_def));

    crosssections_.push_back(PhotonuclearFactory::Get().CreatePhotonuclear(particle_def_,
                                                                           *medium_,
                                                                           cut_settings_,
                                                                           utility_def.photo_def,
                                                                           true,
                                                                           interpolation_def));

    crosssections_.push_back(EpairProductionFactory::Get().CreateEpairProduction(particle_def_,
                                                                                 *medium_,
                                                                                 cut_settings_,
                                                                                 utility_def.epair_def,
                                                                                 interpolation_def));

    crosssections_.push_back(IonizationFactory::Get().CreateIonization(particle_def_,
                                                                       *medium_,
                                                                       cut_settings_,
                                                                       utility_def.ioniz_def,
                                                                       true,
                                                                       interpolation_def));
}

Utility::Utility(const std::vector<CrossSection*>& crosssections)
{
    try
    {
        const ParticleDef& particle_def = crosssections.at(0)->GetParametrization().GetParticleDef();
        const Medium& medium = crosssections.at(0)->GetParametrization().GetMedium();
        const EnergyCutSettings& cuts = crosssections.at(0)->GetParametrization().GetEnergyCuts();

        for (std::vector<CrossSection*>::const_iterator it = crosssections.begin(); it != crosssections.end(); ++it)
        {
            if ((*it)->GetParametrization().GetParticleDef() != particle_def)
            {
                log_fatal("Particle definition of the cross section must be equal!");
            }

            if ((*it)->GetParametrization().GetMedium() != medium)
            {
                log_fatal("Medium of the cross section must be equal!");
            }

            if ((*it)->GetParametrization().GetEnergyCuts() != cuts)
            {
                log_fatal("Energy cuts of the cross section must be equal!");
            }

            crosssections_.push_back((*it)->clone());
        }

        particle_def_ = particle_def;
        medium_ = medium.clone();
        cut_settings_ = cuts;
    }
    catch (const std::out_of_range& e)
    {
        log_fatal("At least one cross section is needed for initializing utility.");
    }
}


Utility::Utility(const Utility& collection)
    :particle_def_(collection.particle_def_)
    ,medium_(collection.medium_->clone())
    ,cut_settings_(collection.cut_settings_)
    ,crosssections_(collection.crosssections_.size(), NULL)
    ,utility_def_(collection.utility_def_)
{
    for (unsigned int i = 0; i < crosssections_.size(); ++i)
    {
        crosssections_[i] = collection.crosssections_[i]->clone();
    }
}

Utility::~Utility()
{
    delete medium_;

    for(std::vector<CrossSection*>::const_iterator iter = crosssections_.begin(); iter != crosssections_.end(); ++iter)
    {
        delete *iter;
    }

    crosssections_.clear();
}

/******************************************************************************
*                               UtilitIntegral                               *
******************************************************************************/

// UtilityIntegral::UtilityIntegral(const ParticleDef& particle_def)
//     : Utility(particle_def)
// {
// }

// UtilityIntegral::UtilityIntegral(const ParticleDef& particle_def,
//                                  const Medium& medium,
//                                  const EnergyCutSettings& cuts,
//                                  const Bremsstrahlung& brems,
//                                  const Photonuclear& photo,
//                                  const Ionization& ioniz,
//                                  const EpairProduction& epair)
//     : Utility(particle_def, medium, cuts)
// {
//     crosssections_.push_back(new BremsIntegral(brems));
//     crosssections_.push_back(new PhotoIntegral(photo));
//     crosssections_.push_back(new IonizIntegral(ioniz));
//     crosssections_.push_back(new EpairIntegral(epair));
// }
//
// UtilityIntegral::UtilityIntegral(const Utility& collection)
//     :Utility(collection)
// {
// }
//
// UtilityIntegral::~UtilityIntegral()
// {
// }

/******************************************************************************
*                             UtilityInterpolant                             *
******************************************************************************/

// UtilityInterpolant::UtilityInterpolant(const ParticleDef& particle_def)
//     : Utility(particle_def)
// {
// }

// UtilityInterpolant::UtilityInterpolant(const ParticleDef& particle_def,
//                                        const Medium& medium,
//                                        const EnergyCutSettings& cuts,
//                                        const Bremsstrahlung& brems,
//                                        const Photonuclear& photo,
//                                        const Ionization& ioniz,
//                                        const EpairProduction& epair,
//                                        InterpolationDef def)
//     : Utility(particle_def, medium, cuts)
//     , interpolation_def(def)
// {
//     crosssections_.push_back(new BremsInterpolant(brems, def));
//     crosssections_.push_back(new PhotoInterpolant(photo, def));
//     crosssections_.push_back(new IonizInterpolant(ioniz, def));
//     crosssections_.push_back(new EpairInterpolant(epair, def));
// }
//
// UtilityInterpolant::UtilityInterpolant(const Utility& collection)
//     :Utility(collection)
// {
// }
//
// UtilityInterpolant::~UtilityInterpolant()
// {
// }

/******************************************************************************
*                            Utility Decorator                            *
******************************************************************************/

UtilityDecorator::UtilityDecorator(const Utility& utility)
    : utility_(utility)
{
}

UtilityDecorator::UtilityDecorator(const UtilityDecorator& decorator)
    : utility_(decorator.utility_)
{

}

UtilityDecorator::~UtilityDecorator()
{
}

// ------------------------------------------------------------------------- //
double UtilityDecorator::FunctionToIntegral( double energy)
{
    double result = 0.0;

    const std::vector<CrossSection*> crosssections = utility_.GetCrosssections();

    for(std::vector<CrossSection*>::const_iterator iter = crosssections.begin(); iter != crosssections.end(); ++iter)
    {
        result += (*iter)->CalculatedEdx(energy);
    }

    return -1.0 / result;
}
