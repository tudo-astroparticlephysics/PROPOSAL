/*! \file   Scattering.cxx
*   \brief  Source file for the scattering routines.
*
*   For more details see the class documentation.
*
*   \date   29.06.2010
*   \author Jan-Hendrik Koehne
*/


#include "PROPOSAL/Scattering.h"
#include <cmath>
#include "algorithm"

#include "PROPOSAL/Propagate.h"
#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Interpolate.h"

#include "PROPOSAL/Bremsstrahlung.h"

using namespace std;

double Scattering::cutoff   =   1;


Scattering::Scattering(PROPOSALParticle *p)
{
    df              =   false;
    jt              =   false;
    this->particle_ =   p;
    propagate_      =   p->propagate_;
    integral_       =   new Integral(IROMB, IMAXS, IPREC2);
}

//----------------------------------------------------------------------------//


double Scattering::function(double E)
{
    const bool DEBUG    =   false;
    double aux, aux2;

    aux     =   propagate_->get_cros()->function(E);
    aux2    =   RY*particle_->e/(particle_->p*particle_->p);
    aux     *=  aux2*aux2;

    if(DEBUG)
    {
        cout <<" $" << o->f(particle_->e) <<" \t "<< o->f(aux);
    }

    return aux;
}

//----------------------------------------------------------------------------//

double Scattering::gettho(double dr, double ei, double ef)
{

    double aux;

    if(jt)
    {
        if(fabs(ei-ef)>fabs(ei)*HALF_PRECISION)
        {
            aux         =   interpolateJ_->interpolate(ei);
            double aux2 =   aux - interpolateJ_->interpolate(ef);

            if(fabs(aux2)>fabs(aux)*HALF_PRECISION)
            {
                aux =   aux2;
            }
            else
            {
                aux =   interpolateJdf_->interpolate((ei+ef)/2)*(ef-ei);
            }
        }
        else
        {
            aux =   interpolateJdf_->interpolate((ei+ef)/2)*(ef-ei);
        }
    }
    else
    {
        aux =   integral_->integrateWithLog(ei, ef, this);
    }

    // setLpm sets xo_ as well
    propagate_->get_cros()->get_bremsstrahlung()->setLpm();

    aux =   sqrt(max(aux, 0.0)/propagate_->get_cros()->get_bremsstrahlung()->get_xo())*particle_->c*max(1 + 0.038*log(dr/propagate_->get_cros()->get_bremsstrahlung()->get_xo()), 0.0);


    return min(aux, cutoff);
}

//----------------------------------------------------------------------------//

double Scattering::functionInt(double e)
{
    if(df)
    {
        return function(e);
    }
    else
    {
        return integral_->integrateWithLog(e, PhysicsModel::ebig_, this);
    }
}

//----------------------------------------------------------------------------//

// Setter

void Scattering::set_jt(bool newJT)
{
    jt  =   newJT;
}

void Scattering::set_df(bool newDF)
{
    df  =   newDF;
}
