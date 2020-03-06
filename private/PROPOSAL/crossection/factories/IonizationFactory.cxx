#include <algorithm>
#include "PROPOSAL/crossection/IonizIntegral.h"
#include "PROPOSAL/crossection/IonizInterpolant.h"
#include "PROPOSAL/crossection/factories/IonizationFactory.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"

#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

IonizationFactory::IonizationFactory()
    : ioniz_map_str_()
    , ioniz_map_enum_()
    , string_enum_()
{
    Register("ionizbetheblochrossi", BetheBlochRossi, &IonizBetheBlochRossi::create);
    Register("ionizbergerseltzerbhabha", IonizBergerSeltzerBhabha, &IonizBergerSeltzerBhabha::create);
    Register("ionizbergerseltzermoller", IonizBergerSeltzerMoller, &IonizBergerSeltzerMoller::create);
    Register("none", None, nullptr);
}

// ------------------------------------------------------------------------- //
IonizationFactory::~IonizationFactory()
{
    ioniz_map_str_.clear();
    ioniz_map_enum_.clear();
    string_enum_.clear();
}

// ------------------------------------------------------------------------- //
void IonizationFactory::Register(const std::string& name, Enum enum_t, RegisterFunction create)
{
    ioniz_map_str_[name]    = create;
    ioniz_map_enum_[enum_t] = create;
    string_enum_.insert(name, enum_t);
}

// ------------------------------------------------------------------------- //
// Most general creator
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
CrossSection* IonizationFactory::CreateIonization(const ParticleDef& particle_def,
                                                  std::shared_ptr<const Medium> medium,
                                                  const EnergyCutSettings& cuts,
                                                  const Definition& def) const
{
    if(def.parametrization == IonizationFactory::Enum::None){
        log_fatal("Can't return Ionization Crosssection if parametrization is None");
        return NULL;
    }

    IonizMapEnum::const_iterator it = ioniz_map_enum_.find(def.parametrization);

    if (it != ioniz_map_enum_.end())
    {
        return new IonizIntegral(*it->second(particle_def, medium, cuts, def.multiplier));
    } else
    {
        log_fatal("Ionization %s not registered!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
CrossSection* IonizationFactory::CreateIonization(const ParticleDef& particle_def,
                                                  std::shared_ptr<const Medium> medium,
                                                  const EnergyCutSettings& cuts,
                                                  const Definition& def,
                                                  InterpolationDef interpolation_def) const
{
    if(def.parametrization == IonizationFactory::Enum::None){
        log_fatal("Can't return Ionization Crosssection if parametrization is None");
        return NULL;
    }

    IonizMapEnum::const_iterator it = ioniz_map_enum_.find(def.parametrization);

    if (it != ioniz_map_enum_.end())
    {
        return new IonizInterpolant(*it->second(particle_def, medium, cuts, def.multiplier), interpolation_def);
    } else
    {
        log_fatal("Ionization %s not registered!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
IonizationFactory::Enum IonizationFactory::GetEnumFromString(const std::string& name)
{
    std::string name_lower = name;
    std::transform(name.begin(), name.end(), name_lower.begin(), ::tolower);

    auto& left = string_enum_.GetLeft();
    auto it = left.find(name_lower);
    if (it != left.end())
    {
        return it->second;
    } else
    {
        log_fatal("Ionization %s not registered!", name.c_str());
        return IonizationFactory::Fail; // Just to prevent warnings

    }
}

// ------------------------------------------------------------------------- //
std::string IonizationFactory::GetStringFromEnum(const IonizationFactory::Enum& enum_t)
{
    auto& right = string_enum_.GetRight();
    auto it = right.find(enum_t);
    if (it != right.end())
    {
        return it->second;
    } else
    {
        log_fatal("Ionization %s not registered!", typeid(enum_t).name());
        return ""; // Just to prevent warnings
    }
}
