
#include <sstream>
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//

double EnergyCutSettings::GetCut(double energy) const
{
    // if both ecut and vcut are setted (in their bounds), the minimum is used
    //
    // if just one of them is setted, this is used
    //
    // if both are setted out of their bound (usually -1, -1),
    // then 1 (which is the max. energy transfer) is returned,
    // which is always greater than v_max
    // that means, just continuous losses are calculated, no stochastic ones

    if (ecut_ > 0)
    {
        if (vcut_ > 0 && vcut_ <= 1)
        {
            return std::min(ecut_ / energy, vcut_);
        } else
        {
            return ecut_ / energy;
        }
    } else
    {
        if (vcut_ > 0 && vcut_ <= 1)
        {
            return vcut_;
        } else
        {
            return 1.;
        }
    }
}

//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//

EnergyCutSettings::EnergyCutSettings()
    : ecut_(500.)
    , vcut_(0.05)
{
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

EnergyCutSettings::EnergyCutSettings(const EnergyCutSettings& cuts)
    : ecut_(cuts.ecut_)
    , vcut_(cuts.vcut_)
{
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

EnergyCutSettings::EnergyCutSettings(const double ecut, const double vcut)
    : ecut_(ecut)
    , vcut_(vcut)
{
}

//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//

EnergyCutSettings& EnergyCutSettings::operator=(const EnergyCutSettings& energyCutSettings)
{
    if (this != &energyCutSettings)
    {
        EnergyCutSettings tmp(energyCutSettings);
        swap(tmp);
    }
    return *this;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

bool EnergyCutSettings::operator==(const EnergyCutSettings& energyCutSettings) const
{
    if (ecut_ != energyCutSettings.ecut_)
        return false;
    if (vcut_ != energyCutSettings.vcut_)
        return false;

    // else
    return true;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

bool EnergyCutSettings::operator!=(const EnergyCutSettings& energyCutSettings) const
{
    return !(*this == energyCutSettings);
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

std::ostream& PROPOSAL::operator<<(std::ostream& os, PROPOSAL::EnergyCutSettings const& cut_settings)
{
    std::stringstream ss;
    ss << " EnergyCutSettings (" << &cut_settings << ") ";
    os << Helper::Centered(60, ss.str()) << '\n';

    os << "Ecut: " << cut_settings.ecut_ << std::endl;
    os << "Vcut: " << cut_settings.vcut_ << std::endl;

    os << Helper::Centered(60, "");
    return os;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

void EnergyCutSettings::swap(EnergyCutSettings& energyCutSettings)
{
    std::swap(ecut_, energyCutSettings.ecut_);
    std::swap(vcut_, energyCutSettings.vcut_);
}
