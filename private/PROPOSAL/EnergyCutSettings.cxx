#include "PROPOSAL/EnergyCutSettings.h"

using namespace std;


//----------------------------------------------------------------------------//

EnergyCutSettings::EnergyCutSettings()
    :ecut_(500.)
    ,vcut_(0.05)
{

}
//----------------------------------------------------------------------------//

EnergyCutSettings::EnergyCutSettings(const EnergyCutSettings &cuts)
{
    *this = cuts;
}
//----------------------------------------------------------------------------//

EnergyCutSettings& EnergyCutSettings::operator=(const EnergyCutSettings &cuts){
    return *this;
}
//----------------------------------------------------------------------------//

EnergyCutSettings::EnergyCutSettings(double ecut, double vcut)
    :ecut_(ecut)
    ,vcut_(vcut)
{

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

void EnergyCutSettings::SetEcut(double ecut){
    ecut_ = ecut;
}

//----------------------------------------------------------------------------//

void EnergyCutSettings::SetVcut(double vcut){
    vcut_ = vcut;
}
