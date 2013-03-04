/*
 * BremsStochastic.cxx
 *
 *  Created on: 21.06.2010
 *      Author: koehne
 */

#include "PROPOSAL/BremsStochastic.h"
#include "PROPOSAL/Medium.h"

using namespace std;

// constructors

BremsStochastic::BremsStochastic()
{
    set_jt(false);
}

//----------------------------------------------------------------------------//

BremsStochastic::BremsStochastic(Bremsstrahlung *cros)
:Bremsstrahlung(*cros)
{

    set_jt(false);
    VectorOfIntegral_.resize(medium_->get_numCompontents());
    h_  =   (double *)calloc(medium_->get_numCompontents(),sizeof(double));

    for(int i =0; i<(medium_->get_numCompontents()); i++)
    {
            VectorOfIntegral_.at(i) =   new Integral(IROMB, IMAXS, IPREC);
    }

}

//----------------------------------------------------------------------------//

// destructors

BremsStochastic::~BremsStochastic()
{

}

//----------------------------------------------------------------------------//

// Memberfunctions


double BremsStochastic::function(double v)
{
    return (cros->get_cb()) * (brems_->Sel(v,cros->get_component()));
}

//----------------------------------------------------------------------------//


double BremsStochastic::dNdx()
{

    if((cros->get_cb())<=0)
    {
        return 0;
    }

    sum_    =   0;

    for(int i=0; i<(medium_->get_numCompontents()); i++)
    {

        if(jt_)
        {
            sum_    +=  max((interpolateJo_[i]).interpolate(particle_->get_energy()), 0.0);
        }
        else
        {
            setEnergy(i);
            sum_    +=  VectorOfIntegral_.at(i)->integrateWithLog(vUp_, vMax_, this);
        }

    }

    return sum_;

}

//----------------------------------------------------------------------------//

double BremsStochastic::dNdx(double rnd)
{

    if((cros->get_cb())<=0)
    {
        return 0.;
    }
    if(jt_)
    {
        rnd_    =   rnd;
    }

    sum_    =   0;

    for(int i =0; i<(medium_->get_numCompontents()); i++)
    {
        if(jt_)
        {
            h_[i]   =   max(interpolateJo_[i].interpolate(particle_->get_energy()) , 0.0);
        }
        else
        {
            setEnergy(i);
            h_[i]   =   VectorOfIntegral_.at(i)->integrateWithLog(vUp_, vMax_, this, rnd);
        }
        sum_    +=  h_[i];
    }

    return sum_;
}

//----------------------------------------------------------------------------//

double BremsStochastic::e(double rnd)
{

    double rand;
    double rsum;

    rand    =   rnd*sum_;
    rsum    =   0;

    for(int i=0; i<(medium_->get_numCompontents()); i++)
    {
        rsum    +=  h_[i];
        if(rsum > rand)
        {

            if(jt_)
            {
                setEnergy(i);

                if(vUp_==vMax_)
                {
                    return (particle_->get_energy())*vUp_;
                }

                return (particle_->get_energy())*(vUp_*exp(interpolateJ_[i].findLimit((particle_->get_energy()), (this->rnd_)*h_[i])*log(vMax_/vUp_)));
            }

            else
            {
                cros->set_component(i);
                return (particle_->get_energy())*VectorOfIntegral_.at(i)->getUpperLimit();
            }

        }
    }

    cerr<<"Error (in BremsStochastic/e): sum was not initialized correctly";
    return 0;

}

//----------------------------------------------------------------------------//

double BremsStochastic::functionInt(double e, double v)
{

    particle_->setEnergy(e);
    setEnergy(cros->get_component());

    if(vUp_==vMax_)
    {
        return 0;
    }

    v   =   vUp_*exp(v*log(vMax_/vUp_));

    return VectorOfIntegral_.at(cros->get_component())->integrateWithLog(vUp_, v, this);


}

//----------------------------------------------------------------------------//

double BremsStochastic::functionInt(double e)
{

    return interpolateJ_[cros->get_component()].interpolate(e, 1.);

}

//----------------------------------------------------------------------------//

// Setter

void BremsStochastic::set_integral(vector<Integral*> integral)
{
    VectorOfIntegral_   =   integral;
}

void BremsStochastic::set_h(double* h)
{
    h_  =   h;
}

void BremsStochastic::set_sum(double sum)
{
    sum_    =   sum;
}

void BremsStochastic::set_interpolateJ(Interpolate* interpolateJ)
{
    interpolateJ_   =   interpolateJ;
}

void BremsStochastic::set_jt(bool newjt)
{
    jt_ =   newjt;
}

void BremsStochastic::set_rnd(double rnd)
{
    rnd_    =   rnd;
}

void BremsStochastic::set_interpolateJo(Interpolate* interpolateJo)
{
    interpolateJo_  =   interpolateJo;
}


