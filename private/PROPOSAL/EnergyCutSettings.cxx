#include "PROPOSAL/EnergyCutSettings.h"
#include <iostream>

using namespace std;


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double EnergyCutSettings::GetCut(double energy)
{

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
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


EnergyCutSettings::EnergyCutSettings()
    :ecut_  ( 500. )
    ,vcut_  ( 0.05 )
{

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


EnergyCutSettings::EnergyCutSettings(const EnergyCutSettings &cuts)
    :ecut_  ( cuts.ecut_ )
    ,vcut_  ( cuts.vcut_ )
{

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


EnergyCutSettings::EnergyCutSettings(double ecut, double vcut)
    :ecut_  ( ecut )
    ,vcut_  ( vcut )
{

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



EnergyCutSettings& EnergyCutSettings::operator=(const EnergyCutSettings &energyCutSettings)
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


bool EnergyCutSettings::operator==(const EnergyCutSettings &energyCutSettings) const
{
    if( ecut_   != energyCutSettings.ecut_)  return false;
    if( vcut_   != energyCutSettings.vcut_)  return false;

    //else
    return true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool EnergyCutSettings::operator!=(const EnergyCutSettings &energyCutSettings) const
{
  return !(*this == energyCutSettings);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


ostream& operator<<(ostream& os, EnergyCutSettings const& cut_settings)
{
    os<<"--------EnergyCutSettings( "<<&cut_settings<<" )--------"<<endl;
    os<<"\tEcut: "<<cut_settings.ecut_<<endl;
    os<<"\tVcut: "<<cut_settings.vcut_<<endl;
    os<<"------------------------------------";
    return os;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void EnergyCutSettings::swap(EnergyCutSettings &energyCutSettings)
{
    using std::swap;

    swap( ecut_   , energyCutSettings.ecut_);
    swap( vcut_   , energyCutSettings.vcut_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void EnergyCutSettings::SetEcut(double ecut){
    ecut_ = ecut;
}

void EnergyCutSettings::SetVcut(double vcut){
    vcut_ = vcut;
}
