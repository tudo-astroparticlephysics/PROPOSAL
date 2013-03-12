/*
 * CrossSections.cxx
 *
 *  Created on: 2013.03.12
 *      Author: koehne
 */

#include "PROPOSAL/CrossSections.h"

using namespace std;

//----------------------------------------------------------------------------//

//destructors

CrossSections::~CrossSections()
{

}

//----------------------------------------------------------------------------//

//constructors

CrossSections::CrossSections()
{

}

//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//

//Memberfunctions

double CrossSections::functionInt(double e)
{

//    if(df_)
//    {
//        return function(e);
//    }
//    else
//    {
//        return integral_->integrateWithLog(e, particle_->low, this);
//    }
      return 65.3;

}

//----------------------------------------------------------------------------//

double CrossSections::function(double E)
{

//    const bool DEBUG    =   false;
//    double result;
//    double aux;

//    particle_->setEnergy(E);
//    result  =    0;

//    if(DEBUG)
//    {
//        cerr<<" * "<<o->f((particle_->get_energy()));
//    }

//    aux     =   ionization_->get_Continuous()->dEdx();
//    result  +=  aux;

//    if(DEBUG)
//    {
//        cerr<<" \t "<<o->f(aux);
//    }

//    aux     =   bremsstrahlung_->get_Continuous()->dEdx();
//    result  +=  aux;

//    if(DEBUG)
//    {
//        cerr<<" \t "<<o->f(aux);
//    }

//    aux     =   photonuclear_->get_Continuous()->dEdx();
//    result  +=  aux;

//    if(DEBUG)
//    {
//        cerr<<" \t "<<o->f(aux);
//    }

//    aux     =   epairproduction_->get_Continuous()->dEdx();
//    result  +=  aux;

//    if(DEBUG)
//    {
//        cerr<<" \t "<<o->f(aux);
//    }

//    return -1/result;
      return 23.3;
}

//----------------------------------------------------------------------------//


double CrossSections::getdx(double ei, double ef, double dist)
{
//    if(jt_)
//    {
//        if(fabs(ei-ef) > fabs(ei)*HALF_PRECISION)
//        {
//            double aux;

//            ini_    =   interpolateJ_->interpolate(ei);
//            aux     =   ini_ - interpolateJ_->interpolate(ef);

//            if(fabs(aux) > fabs(ini_)*HALF_PRECISION)
//            {
//                return max(aux, 0.0);
//            }

//        }

//        ini_    =   0;
//        return max((interpolateJdf_->interpolate((ei + ef)/2))*(ef - ei), 0.0);

//    }
//    else
//    {
//        return integral_->integrateWithLog(ei, ef, this, -dist);
//    }
      return 2.4;
}

//----------------------------------------------------------------------------//

double CrossSections::getef(double ei, double dist)
{

//    if(jt_)
//    {
//        if(ini_ != 0)
//        {
//            double aux;
//            aux     =   interpolateJ_->findLimit(ini_-dist);

//            if(fabs(aux) > fabs(ei)*HALF_PRECISION)
//            {
//                return min( max(aux, particle_->get_low()), ei);
//            }
//        }

//        return min( max(ei+dist/interpolateJdf_->interpolate(ei + dist/(2*interpolateJdf_->interpolate(ei))), particle_->get_low()), ei);
//    }
//    else
//    {
//        return integral_->getUpperLimit();
//    }
      return 82.;

}

//----------------------------------------------------------------------------//


void CrossSections::SetConstants()
{
//    set_component(0);

//    ci_         =   1.;
//    cb_         =   1.;
//    cp_         =   1.;
//    ce_         =   1.;
//    cd_         =   1.;
//    cm_         =   1.;
//    bremserror  =   0;
//    epairerror  =   0;
//    photoerror  =   0;
//    ionizerror  =   0;
//    df_         =   false;
//    jt_         =   false;
//    lpm_        =   false;
//    e_hi        =   PhysicsModel::ebig_;
//    g           =   5;

}

//----------------------------------------------------------------------------//

// Getter



