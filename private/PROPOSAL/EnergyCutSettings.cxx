#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/methods.h"

#include <sstream>
#include <stdexcept>

using namespace PROPOSAL;

//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//

double EnergyCutSettings::GetCut(double energy) const
{
    // if both ecut and vcut are setted (in their bounds), the minimum is used
    return std::min(ecut_ / energy, vcut_);
}

//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//


EnergyCutSettings::EnergyCutSettings(
    const double ecut, const double vcut, const bool continuous_randomization)
    : ecut_(ecut)
    , vcut_(vcut)
    , continuous_randomization_(continuous_randomization)
{
    if (ecut_ < 0)
        throw std::invalid_argument(
            "Absolut EnergyCut (ecut) must be larger > 0! \nINFO: If a ecut = "
            "std::numeric_limits<double>::infinity() is set the particle will "
            "only propagated continously. ");

    if (vcut_ <= 0 || vcut_ >= 1)
        throw std::invalid_argument(
            "Relative EnergyCut (vcut) must be between [0, 1]. \nINF0: If a vc"
            "ut = 1 is set, the particle is only propagate continous.");

    log_debug(
        "EnergyCut set to ecut(%f), vcut(%f), continuous_randomization(%b)",
        ecut_, vcut_, continuous_randomization_);
};

EnergyCutSettings::EnergyCutSettings(const nlohmann::json& config)
{
    if(!config.contains("e_cut"))
        throw std::invalid_argument("A e_cut must be defined. e_cut = [0-inf]");
    if(!config.contains("v_cut"))
        throw std::invalid_argument("A v_cut must be defined. v_cut = [0-1]");
    if(!config.contains("cont_rand"))
        throw std::invalid_argument("cont_rand must be activated/deactivated. cont_rand = true/false");

    double ecut, vcut;
    bool cont_rand;
    config.at("e_cut").get_to(ecut);
    config.at("v_cut").get_to(vcut);
    config.at("cont_rand").get_to(cont_rand);

    *this = EnergyCutSettings(ecut, vcut, cont_rand);
};

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

bool EnergyCutSettings::operator==(
    const EnergyCutSettings& energyCutSettings) const
{
    if (ecut_ != energyCutSettings.ecut_)
        return false;
    if (vcut_ != energyCutSettings.vcut_)
        return false;
    if (continuous_randomization_
        != energyCutSettings.continuous_randomization_)
        return false;

    return true;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

std::ostream& PROPOSAL::operator<<(
    std::ostream& os, PROPOSAL::EnergyCutSettings const& cut_settings)
{
    std::stringstream ss;
    ss << " EnergyCutSettings (" << &cut_settings << ") ";
    os << Helper::Centered(60, ss.str()) << '\n';

    os << "Ecut: " << cut_settings.ecut_ << std::endl;
    os << "Vcut: " << cut_settings.vcut_ << std::endl;
    os << "Continuous randomization: " << cut_settings.continuous_randomization_
       << std::endl;

    os << Helper::Centered(60, "");
    return os;
}
