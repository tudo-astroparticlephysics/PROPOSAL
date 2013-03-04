/*
 * EpairContinuous.cxx
 *
 *  Created on: 22.06.2010
 *      Author: koehne
 */

#include "PROPOSAL/EpairContinuous.h"
#include <algorithm>
#include <cmath>
#include "PROPOSAL/Medium.h"


using namespace std;

EpairContinuous::EpairContinuous(Epairproduction *cros)
:Epairproduction(*cros)
{
    reverse_    =   false;
    jt_         =   false;
    integral_   =   new Integral(IROMB, IMAXS, IPREC);
}

//----------------------------------------------------------------------------//

double EpairContinuous::function(double v)
{
    if(reverse_)
    {
        v   =   1-v;
    }

    return v*epair->ePair(v, cros->get_component());
}

//----------------------------------------------------------------------------//



double EpairContinuous::dEdx()
{
    if(cros->get_ce()<=0)
    {
        return 0;
    }

    if(jt_)
    {
        return max(interpolateJ_->interpolate(particle_->e), 0.0);
    }
    int i;

    double sum  =   0;

    for(i=0; i<medium_->get_numCompontents(); i++)
    {
        setEnergy(i);

        double r1   =   0.8;
        double rUp  =   vUp*(1-HALF_PRECISION);
        bool rflag  =   false;

        if(r1<rUp)
        {
            if(2*function(r1)<function(rUp))
            {
                rflag   =   true;
            }
        }

        if(rflag)
        {
            if(r1>vUp)
            {
                r1  =   vUp;
            }

            if(r1<vMin)
            {
                r1  =   vMin;
            }

            sum         +=  integral_->integrateWithLog(vMin, r1, this);
            reverse_    =   true;
            double r2   =   max(1-vUp, COMPUTER_PRECISION);

            if(r2>1-r1)
            {
                r2  =   1-r1;
            }

            sum         +=  integral_->integrateOpened(1-vUp, r2, this)+integral_->integrateWithLog(r2, 1-r1, this);
            reverse_    =   false;
        }

        else
        {
            sum +=  integral_->integrateWithLog(vMin, vUp, this);
        }
    }

    return cros->get_ce()*particle_->e*sum;
}

//----------------------------------------------------------------------------//

double  EpairContinuous::functionInt(double e)
{
    particle_->setEnergy(e);
    return dEdx();
}


