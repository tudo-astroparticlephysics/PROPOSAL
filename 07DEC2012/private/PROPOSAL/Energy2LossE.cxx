/*
 * Energy2LossE.cxx
 *
 *  Created on: 24.06.2010
 *      Author: koehne
 */
#include "PROPOSAL/Energy2LossE.h"
#include <cmath>
#include "algorithm"
#include "PROPOSAL/Energy2LossX.h"

using namespace std;

Energy2LossE::Energy2LossE(Energy2Loss *e2loss)
:Energy2Loss(e2loss)
{
    df          =   false;
    jt          =   false;
    integral_   =   new Integral(IROMB, IMAXS, IPREC2);
}

//----------------------------------------------------------------------------//

double Energy2LossE::function(double E)
{
    double aux;
    aux =   cros->function(E);
    return aux*e2loss->e2lx->dE2dx();
}

//----------------------------------------------------------------------------//

double Energy2LossE::dE2de(double ei, double ef)
{
    if(jt)
    {
        if(abs(ei-ef)>abs(ei)*HALF_PRECISION)
            {
            double aux  =   interpolateJ_->interpolate(ei);
            double aux2 =   aux-interpolateJ_->interpolate(ef);

            if(abs(aux2)>abs(aux)*HALF_PRECISION)
            {
                return max(aux2, 0.0);
            }
        }

        else
        {
            return max(interpolateJdf_->interpolate((ei+ef)/2)*(ef-ei), 0.0);
        }

    }
    return integral_->integrateWithLog(ei, ef, this);

}

//----------------------------------------------------------------------------//

double Energy2LossE::functionInt(double e)
{
    if(df)
    {
        return function(e);
    }
    else
    {
        return integral_->integrateWithLog(e, particle_->low, this);
    }
}

//----------------------------------------------------------------------------//

//Setter

void Energy2LossE::set_jt(bool newJT)
{
    jt  =   newJT;
}

void Energy2LossE::set_df(bool newDF)
{
    df  =   newDF;
}




