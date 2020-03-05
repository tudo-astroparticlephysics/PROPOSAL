
#include <algorithm>

#include "PROPOSAL/crossection/MupairIntegral.h"
#include "PROPOSAL/crossection/MupairInterpolant.h"
#include "PROPOSAL/crossection/factories/MupairProductionFactory.h"
#include "PROPOSAL/crossection/parametrization/MupairProduction.h"

#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

MupairProductionFactory::MupairProductionFactory()
    : mupair_map_str_()
    , mupair_map_enum_()
    , string_enum_()
{
    Register("mupairkelnerkokoulinpetrukhin", KelnerKokoulinPetrukhin, std::make_pair(&MupairKelnerKokoulinPetrukhin::create, &MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin>::create));
    Register("none", None, std::make_pair(nullptr, nullptr));
}

MupairProductionFactory::~MupairProductionFactory()
{
    string_enum_.clear();
    mupair_map_str_.clear();
    mupair_map_enum_.clear();
}

// ------------------------------------------------------------------------- //
// Most general creator
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
CrossSection* MupairProductionFactory::CreateMupairProduction(const ParticleDef& particle_def,
                                                            std::shared_ptr<const Medium> medium,
                                                            std::shared_ptr<const EnergyCutSettings> cuts,
                                                            const Definition& def,
                                                            std::shared_ptr<const InterpolationDef> interpolation_def = nullptr) const
{
    if(def.parametrization == MupairProductionFactory::Enum::None){
        log_fatal("Can't return MuPairproduction Crosssection if parametrization is None");
        return NULL;
    }

    MupairMapEnum::const_iterator it = mupair_map_enum_.find(def.parametrization);

    if (it != mupair_map_enum_.end())
    {
        if(interpolation_def==nullptr){
            return new MupairIntegral(*it->second.first(particle_def, medium, def.multiplier), cuts);
        }
        else{
            return new MupairInterpolant(*it->second.second(particle_def, medium, def.multiplier, interpolation_def), cuts, interpolation_def);
        }
    } else
    {
        log_fatal("MupairProduction %s not registerd!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warnings
    }
}

CrossSection* MupairProductionFactory::CreateMupairProduction(const MupairProduction& parametrization,
                                                              std::shared_ptr<const EnergyCutSettings> cuts,
                                                              std::shared_ptr<const InterpolationDef> interpolation_def=nullptr) const
{
    if(interpolation_def==nullptr){
        return new MupairIntegral(parametrization, cuts);
    }
    else{
        return new MupairInterpolant(parametrization, cuts, interpolation_def);
    }
}

// ------------------------------------------------------------------------- //
void MupairProductionFactory::Register(const std::string& name,
                                     Enum enum_t,
                                     std::pair<RegisterFunction, RegisterFunctionInterpolant> create)
{
    mupair_map_str_[name]    = create;
    mupair_map_enum_[enum_t] = create;
    string_enum_.insert(name, enum_t);
}

// ------------------------------------------------------------------------- //
MupairProductionFactory::Enum MupairProductionFactory::GetEnumFromString(const std::string& name)
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
        log_fatal("MupairProduction %s not registerd!", name.c_str());
        return Fail; // Just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
std::string MupairProductionFactory::GetStringFromEnum(const MupairProductionFactory::Enum& enum_t)
{
    auto& right = string_enum_.GetRight();
    auto it = right.find(enum_t);
    if (it != right.end())
    {
        return it->second;
    } else
    {
        log_fatal("MupairProduction %s not registerd!", typeid(enum_t).name());
        return ""; // Just to prevent warnings
    }
}
