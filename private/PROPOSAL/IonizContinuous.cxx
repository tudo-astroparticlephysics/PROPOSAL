/*! \file   IonizContinuous.cxx
*   \brief  Source file for the continuous ionizationloss routines.
*
*   For more details see the class documentation.
*
*   \date   24.06.2010
*   \author Jan-Hendrik Koehne
*/


#include "PROPOSAL/IonizContinuous.h"
#include "PROPOSAL/IonizStochastic.h"
#include <cmath>
#include "algorithm"
#include "PROPOSAL/Medium.h"

using namespace std;


IonizContinuous::IonizContinuous(Ionizationloss *cros)
:Ionizationloss(*cros)
{
    jt_         =   false;
    integral_   =   new Integral(IROMB, IMAXS, IPREC);
}

//----------------------------------------------------------------------------//


double  IonizContinuous::function(double v)
{
    return v*ioniz->d2Ndvdx(v)*ioniz->inelCorrection(v);
}

//----------------------------------------------------------------------------//


double IonizContinuous::delta()
{
    double X;

    X   =   log(beta*gamma)/log(10);

    if(X<medium_->get_X0())
    {
        return medium_->get_d0()*pow(10 , 2*(X-medium_->get_X0()));
    }
    else if(X<medium_->get_X1())
    {
        return medium_->get_C1()*X + medium_->get_C() + medium_->get_a()*pow(medium_->get_X1() - X , medium_->get_m());
    }
    else
    {
        return medium_->get_C1()*X + medium_->get_C();
    }
}

//----------------------------------------------------------------------------//



double IonizContinuous::dEdx()
{
    if(cros->get_ci()<=0)
    {
        return 0;
    }

    if(jt_)
    {
        return max(interpolateJ_->interpolate(particle_->e), 0.);
    }

    double result, aux;
    setEnergy();

    aux     =   beta*gamma/(1.e-6*medium_->get_I());
    result  =   log(vUp*(2*ME*particle_->e))+2*log(aux);
    aux     =   vUp/(2*(1 + 1/gamma));
    result  +=  aux*aux;
    aux     =   beta*beta;
    result  -=  aux*(1 + vUp/vMax) + delta();

    if(result>0)
    {
        result*=IONK*particle_->c*particle_->c*medium_->get_ZA()/(2*aux);
    }
    else
    {
        result=0;
    }
    return cros->get_ci()*(medium_->get_massDensity()*(1 + ionizerror)*result + particle_->e*(integral_->integrateWithLog(vMin, vUp, this)));
}

//----------------------------------------------------------------------------//


double IonizContinuous::functionInt(double e)
{
    particle_->setEnergy(e);
    return dEdx();
}

//----------------------------------------------------------------------------//
// Setter

void IonizContinuous::set_jt(bool newJT)
{
    jt_ = newJT;
}

