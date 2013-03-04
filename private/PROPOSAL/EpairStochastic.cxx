/*
 * EpairStochastic.cxx
 *
 *  Created on: 29.09.2010
 *      Author: koehne
 */

#include "PROPOSAL/EpairStochastic.h"
#include "PROPOSAL/Medium.h"
#include <cmath>
#include <algorithm>


using namespace std;


EpairStochastic::EpairStochastic(Epairproduction *cros)
:Epairproduction(*cros)
{

    jt_ =   false;

    int i;
    integral_.resize(medium_->get_numCompontents());

    for(i=0; i<medium_->get_numCompontents(); i++)
    {
        integral_.at(i) =   new Integral(IROMB, IMAXS, IPREC);
    }

    H_  =   (double *)calloc(medium_->get_numCompontents(),sizeof(double));


}

//----------------------------------------------------------------------------//


double EpairStochastic::function(double v)
{
    return cros->get_ce()*epair->ePair(v, cros->get_component());
}

//----------------------------------------------------------------------------//

double EpairStochastic::functionInt(double e, double v)
{
    particle_->setEnergy(e);
    setEnergy(cros->get_component());

    if(vUp==vMax)
    {
        return 0.;
    }

    v   =   vUp*exp(v*log(vMax/vUp));

    return integral_.at(cros->get_component())->integrateWithLog(vUp, v, this);

}

//----------------------------------------------------------------------------//

double EpairStochastic::functionInt(double e)
{
    return interpolateJ_[cros->get_component()].interpolate(e, 1.);
}

//----------------------------------------------------------------------------//

double EpairStochastic::dNdx()
{
    if(cros->get_ce()<=0)
    {
        return 0;
    }

    int i;
    sum =   0;

    for(i=0; i<medium_->get_numCompontents(); i++)
    {
        if(jt_)
        {
            sum +=  max(interpolateJo_[i].interpolate(particle_->e), 0.0);
        }
        else
        {
            setEnergy(i);
            sum +=  integral_.at(i)->integrateWithLog(vUp, vMax, this);
        }
    }

    return sum;
}

//----------------------------------------------------------------------------//


double EpairStochastic::dNdx(double rnd)
{

    if(cros->get_ce()<=0)
    {
        return 0.;
    }

    if(jt_)
    {
        this->rnd   =   rnd;
    }

    int i;
    sum =   0;

    for(i=0; i<medium_->get_numCompontents(); i++)
    {
        if(jt_)
        {
            H_[i]   =   max(interpolateJo_[i].interpolate(particle_->e), 0.0);
        }
        else
        {
            setEnergy(i);
            H_[i]   =   integral_.at(i)->integrateWithLog(vUp, vMax, this, rnd);
        }

        sum +=  H_[i];
    }
    return sum;
}

//----------------------------------------------------------------------------//



double EpairStochastic::e(double rnd)
{
    int i;
    double rand, rsum;

    rand    =   rnd*sum;
    rsum    =   0;

    for(i=0; i<medium_->get_numCompontents(); i++)
    {
        rsum    +=  H_[i];
        if(rsum>rand)
        {
            if(jt_)
            {
                setEnergy(i);

                if(vUp==vMax)
                {
                    return particle_->e*vUp;
                }

                return particle_->e*(vUp*exp(interpolateJ_[i].findLimit(particle_->e, this->rnd*H_[i])*log(vMax/vUp)));
            }

            else
            {
                cros->set_component(i);
                return particle_->e*integral_.at(i)->getUpperLimit();
            }
        }
    }

    cerr<<"Error (in EpairStochastic/e): sum was not initialized correctly"<<endl;

    return 0.;
}


