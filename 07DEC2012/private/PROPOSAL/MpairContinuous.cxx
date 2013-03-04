/*! \file   MpairContinuous.cxx
*   \brief  Source file for the continuous pairproduction routines.
*
*   \author Jan-Hendrik Koehne
*/

#include "PROPOSAL/MpairContinuous.h"
#include "PROPOSAL/Medium.h"

using namespace std;


MpairContinuous:: MpairContinuous(Mpairproduction* MPAIR)
:Mpairproduction(*MPAIR)
{
    jt_         =   false;
    integral_   =   new Integral(IROMB, IMAXS, IPREC);
    cm_         =   1;
}

//----------------------------------------------------------------------------//

double MpairContinuous::dEdx()
{
    if(cm_ <= 0)
    {
        return 0;
    }

    int i;
    double sum  =   0;

    for(i=0; i<medium_->get_numCompontents(); i++)
    {
        if(get_jt())
        {
            return max(interpolateJ_->interpolate(particle_->e), 0.0);
        }

        else
        {
            setEnergy(i);
        }
        sum +=  integral_->integrateWithLog(vMin_, vUp_, this);
    }

    return particle_->get_energy()*sum;
}

//----------------------------------------------------------------------------//

double MpairContinuous::function(double v)
{
    double aux;

    aux =   v*Mpair_->mpair(v,get_component());
    return aux;
}

//----------------------------------------------------------------------------//

void MpairContinuous::activate_interpolation()
{
        set_jt(false);
        interpolateJ_   =   new Interpolate(NUM1, e_low, e_hi, this, g, true, false, true, g, false, false, false);
        set_jt(true);
}

//----------------------------------------------------------------------------//

double  MpairContinuous::functionInt(double e)
{
    particle_->setEnergy(e);
    return dEdx();
}
