/*! \file   MpairStochastic.cxx
*   \brief  Source file for the stochastic pairproduction routines.
*
*   \author Jan-Hendrik Koehne
*/

#include "PROPOSAL/MpairStochastic.h"
#include "PROPOSAL/Medium.h"

using namespace std;



MpairStochastic:: MpairStochastic(Mpairproduction* MPAIR)
:Mpairproduction(*MPAIR)
{
    jt_ =   false;
    cm_ =   1;

    integral_.resize(medium_->get_numCompontents());

    for(int i=0; i<medium_->get_numCompontents(); i++)
    {
        integral_.at(i) = new Integral(IROMB, IMAXS, IPREC);
    }

    H_.resize(medium_->get_numCompontents());
}

//----------------------------------------------------------------------------//

double MpairStochastic::function(double v)
{
    return cm_*Mpair_->mpair(v, get_component());
}

//----------------------------------------------------------------------------//

double MpairStochastic::dNdx()
{
    if(get_cm()<=0)
    {
        return 0;
    }

    sum =   0;

    for(int i=0; i<medium_->get_numCompontents(); i++)
    {
        if(get_jt())
        {
            sum +=  max(interpolateJo_.at(i)->interpolate(particle_->e), 0.0);
        }
        else
        {
            setEnergy(i);
            sum +=  integral_.at(i)->integrateWithLog(vUp_, vMax_, this);
        }
    }
    return sum;
}

//----------------------------------------------------------------------------//

double MpairStochastic::dNdx(double rnd)
{
    if(get_cm()<=0)
    {
        return 0.;
    }
    if(get_jt())
    {
        this->rnd   =   rnd;
    }

    sum =   0;

    for(int i=0; i<medium_->get_numCompontents(); i++)
    {
        if(jt_)
        {
            H_.at(i)    =   max(interpolateJo_.at(i)->interpolate(particle_->e), 0.0);
        }
        else
        {
            setEnergy(i);
            H_.at(i)    =   integral_.at(i)->integrateWithLog(vUp_, vMax_, this, rnd);
        }

        sum +=  H_.at(i);
    }
    return sum;
}

//----------------------------------------------------------------------------//

double MpairStochastic::e(double rnd)
{
    double rand, rsum;

    rand    =   rnd*sum;
    rsum    =   0;

    for(int i=0; i<medium_->get_numCompontents(); i++)
    {
        rsum    +=  H_.at(i);

        if(rsum>rand)
        {
            if(get_jt())
            {
                setEnergy(i);

                if(vUp_==vMax_)
                {
                    return particle_->e*vUp_;
                }
                    return particle_->e*(vUp_*exp(interpolateJ_.at(i)->findLimit(particle_->e, this->rnd*H_.at(i))*log(vMax_/vUp_)));
            }
            else
            {
                set_component(i);
                return particle_->e*integral_.at(i)->getUpperLimit();
            }
        }
    }

    cerr<<"Error (in MpairStochastic/e): sum was not initialized correctly"<<endl;

    return 0.;
}

//----------------------------------------------------------------------------//

double MpairStochastic::functionInt(double e, double v)
{
    particle_->setEnergy(e);
    setEnergy(get_component());

    if(vUp_==vMax_)
    {
        return 0.;
    }

    v   =   vUp_*exp(v*log(vMax_/vUp_));

    return integral_.at(get_component())->integrateWithLog(vUp_, v, this);

}

//----------------------------------------------------------------------------//

double MpairStochastic::functionInt(double e)
{
    return interpolateJ_.at(get_component())->interpolate(e, 1.);
}

//----------------------------------------------------------------------------//

void MpairStochastic::activate_interpolation()
{
    set_jt(false);

    interpolateJ_.resize(medium_->get_numCompontents());
    interpolateJo_.resize(medium_->get_numCompontents());

    for(int i=0; i<medium_->get_numCompontents(); i++)
    {
        set_component(i);
        interpolateJ_.at(i)     =   new Interpolate(NUM1, e_low, e_hi, NUM1, 0, 1, this, g, false, false, true, g, false, false, false, g, true, false, false);
        interpolateJo_.at(i)    =   new Interpolate(NUM1, e_low, e_hi, this, g, false, false, true, g, true, false, false);
    }

    set_jt(true);
}
