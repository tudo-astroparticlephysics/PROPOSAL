/*
 * BremsContinuous.cxx
 *
 *  Created on: 21.06.2010
 *      Author: koehne
 */

#include "PROPOSAL/BremsContinuous.h"
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/PROPOSALParticle.h"
#include <algorithm>

using namespace std;

// constructors

BremsContinuous::BremsContinuous()
{
    set_jt(false);
}

//----------------------------------------------------------------------------//

BremsContinuous::BremsContinuous(Bremsstrahlung *cros)
: Bremsstrahlung (*cros)
{
    set_jt(false);
    integral_   =   new Integral(IROMB, IMAXS, IPREC);
}

//----------------------------------------------------------------------------//

// destructors

BremsContinuous::~BremsContinuous()
{

}

//----------------------------------------------------------------------------//


	// Memberfunctions

double BremsContinuous::function(double v)
{
    return v*brems_->Sel(v, cros->get_component());
}

//----------------------------------------------------------------------------//


double BremsContinuous::dEdx()
{
    double sum  =   0;

    if((cros->get_cb())<=0)
    {
        return 0;
    }

    if(jt_)
    {
        return max(interpolateJ_->interpolate(particle_->get_energy()), 0.0);
    }

    for(int i=0; i<(medium_->get_numCompontents()); i++)
    {
        setEnergy(i);
        sum +=  integral_->integrateOpened(0, vUp_, this);
    }

    return (cros->get_cb())*(particle_->get_energy())*sum;

}

//----------------------------------------------------------------------------//

double BremsContinuous::functionInt(double e)
{
    particle_->setEnergy(e);
    return dEdx();
}

//----------------------------------------------------------------------------//


// Setter

void BremsContinuous::set_integral(Integral *integral)
{
    integral_   =   integral;
}

void BremsContinuous::set_interpolateJ(Interpolate* interpolateJ)
{
    interpolateJ_   =   interpolateJ;
}

void BremsContinuous::set_jt(bool newjt)
{
    jt_ =   newjt;
}

