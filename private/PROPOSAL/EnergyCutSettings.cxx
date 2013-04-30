#include "PROPOSAL/EnergyCutSettings.h"
#include <iostream>

using namespace std;


//----------------------------------------------------------------------------//

EnergyCutSettings::EnergyCutSettings()
    :ecut_(500.)
    ,vcut_(0.05)
{

}
//----------------------------------------------------------------------------//

EnergyCutSettings::EnergyCutSettings(const EnergyCutSettings &cuts)
    :ecut_(cuts.ecut_)
    ,vcut_(cuts.vcut_)
{

}
//----------------------------------------------------------------------------//

EnergyCutSettings::EnergyCutSettings(double ecut, double vcut)
    :ecut_(ecut)
    ,vcut_(vcut)
{

}
//----------------------------------------------------------------------------//

EnergyCutSettings& EnergyCutSettings::operator=(const EnergyCutSettings &energyCutSettings){
    if (this != &energyCutSettings)
    {
      EnergyCutSettings tmp(energyCutSettings);
      swap(tmp);
    }
    return *this;
}

//----------------------------------------------------------------------------//
bool EnergyCutSettings::operator==(const EnergyCutSettings &nergyCutSettings) const
{
    if( ecut_   != nergyCutSettings.ecut_)  return false;
    if( vcut_   != nergyCutSettings.vcut_)  return false;

    //else
    return true;
}
//----------------------------------------------------------------------------//

bool EnergyCutSettings::operator!=(const EnergyCutSettings &energyCutSettings) const
{
  return !(*this == energyCutSettings);
}


//----------------------------------------------------------------------------//

double EnergyCutSettings::GetCut(double energy){

    double aux;
    double vCut;

    aux =   ecut_/energy;

    if(ecut_>0)
    {
        if(vcut_>0 && vcut_<=1)
        {
            if(aux<vcut_)
            {
                vCut   =   aux;
            }
            else
            {
                vCut   =   vcut_;
            }
        }
        else
        {
            vCut   =   aux;
        }
    }
    else
    {
        if(vcut_>0 && vcut_<=1)
        {
            vCut   =   vcut_;
        }
        else
        {
            vCut   =   1.;
        }
    }
    return vCut;
}

//----------------------------------------------------------------------------//

void EnergyCutSettings::swap(EnergyCutSettings &energyCutSettings)
{
    using std::swap;

    swap( ecut_   , energyCutSettings.ecut_);
    swap( vcut_   , energyCutSettings.vcut_);
}

//----------------------------------------------------------------------------//

void EnergyCutSettings::SetEcut(double ecut){
    ecut_ = ecut;
}

//----------------------------------------------------------------------------//

void EnergyCutSettings::SetVcut(double vcut){
    vcut_ = vcut;
}
