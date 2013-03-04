/*! \file   PhotoStochastic.cxx
*   \brief  Source file for the stochastic photonuclear routines.
*
*   For more details see the class documentation.
*
*   \date   23.06.2010
*   \author Jan-Hendrik Koehne
*/


#include "PROPOSAL/PhotoStochastic.h"
#include "PROPOSAL/Medium.h"

using namespace std;

PhotoStochastic::PhotoStochastic(Photonuclear *cros)
:Photonuclear(*cros)
{

    jt_     =   false;

    integral_.resize(medium_->get_numCompontents());

    for(int i=0; i<medium_->get_numCompontents(); i++)
    {
        integral_.at(i) =   new Integral(IROMB, IMAXS, IPREC);
    }

    H_  =   (double *)calloc(medium_->get_numCompontents(),sizeof(double));

}

//----------------------------------------------------------------------------//

double PhotoStochastic::function(double v)
{
    return cros->get_cp()*photo->photoN(v, cros->get_component());
}

//----------------------------------------------------------------------------//

double PhotoStochastic::dNdx()
{

    if(cros->get_cp()<=0)
    {
        return 0;
    }

    sum_    =   0;

    for(int i=0; i<medium_->get_numCompontents(); i++)
    {

        if(jt_)
        {
            sum_    +=  max(interpolateJo_[i].interpolate(particle_->e), 0.0);
        }
        else
        {
            setEnergy(i);
            sum_    +=  integral_.at(i)->integrateWithLog(vUp, vMax, this);
        }
    }

    return sum_;
}

//----------------------------------------------------------------------------------------------------//

double PhotoStochastic::dNdx(double rnd)
{

    if(cros->get_cp()<=0)
    {
        return 0.;
    }

    if(jt_)
    {
        this->rnd   =   rnd;
    }

    sum_    =   0;

    for(int i=0; i<medium_->get_numCompontents(); i++)
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

        sum_    +=  H_[i];
    }

    return sum_;
}

//----------------------------------------------------------------------------------------------------//

double PhotoStochastic::e(double rnd)
{

    double rand, rsum;

    rand    =   rnd*sum_;
    rsum    =   0;

    for(int i=0; i<medium_->get_numCompontents(); i++)
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

    cerr<<"Error (in PhotoStochastic/e): sum was not initialized correctly";
    return 0;
}

//----------------------------------------------------------------------------//

double PhotoStochastic::functionInt(double e, double v)
{
    particle_->setEnergy(e);
    setEnergy(cros->get_component());

    if(vUp==vMax)
    {
        return 0;
    }

    v   =   vUp*exp(v*log(vMax/vUp));

    return integral_.at(cros->get_component())->integrateWithLog(vUp, v, this);
}

//----------------------------------------------------------------------------//

double PhotoStochastic::functionInt(double e)
{
    return interpolateJ_[cros->get_component()].interpolate(e, 1.);
}


