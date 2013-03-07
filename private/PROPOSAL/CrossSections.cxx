/*
 * CrossSections.cxx
 *
 *  Created on: 21.06.2010
 *      Author: koehne
 */

#include "PROPOSAL/CrossSections.h"
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/Interpolate.h"
#include "PROPOSAL/Ionizationloss.h"
#include "PROPOSAL/Epairproduction.h"
#include "PROPOSAL/Decay.h"
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Photonuclear.h"
#include "PROPOSAL/IonizContinuous.h"
#include "PROPOSAL/BremsContinuous.h"
#include "PROPOSAL/PhotoContinuous.h"
#include "PROPOSAL/PhotoStochastic.h"
#include "PROPOSAL/EpairContinuous.h"


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
    SetConstants();
}

//----------------------------------------------------------------------------//

CrossSections::CrossSections(CrossSections *cros)
{
    SetConstants();
    this->particle_ =   cros->particle_;
    this->medium_   =   cros->medium_;
    this->cros      =   cros;
}

//----------------------------------------------------------------------------//

CrossSections::CrossSections(PROPOSALParticle *p, Medium *m)
{

    SetConstants();
    this->particle_         =   p;
    e_low                   =   particle_->low;
    this->medium_           =   m;
    this->cros              =   this;
    this->decay_            =   new Decay(this);
    this->ionization_       =   new Ionizationloss(this);
    this->bremsstrahlung_   =   new Bremsstrahlung(this);
    photonuclear_           =   new Photonuclear(this);
    this->epairproduction_  =   new Epairproduction(this);
    this->integral_         =   new Integral(IROMB, IMAXS, IPREC2);

}

//----------------------------------------------------------------------------//

//Memberfunctions

double CrossSections::functionInt(double e)
{

    if(df_)
    {
        return function(e);
    }
    else
    {
        return integral_->integrateWithLog(e, particle_->low, this);
    }

}

//----------------------------------------------------------------------------//

double CrossSections::function(double E)
{

    const bool DEBUG    =   false;
    double result;
    double aux;

    particle_->setEnergy(E);
    result  =    0;

    if(DEBUG)
    {
        cerr<<" * "<<o->f((particle_->get_energy()));
    }

    aux     =   ionization_->get_Continuous()->dEdx();
    result  +=  aux;

    if(DEBUG)
    {
        cerr<<" \t "<<o->f(aux);
    }

    aux     =   bremsstrahlung_->get_Continuous()->dEdx();
    result  +=  aux;

    if(DEBUG)
    {
        cerr<<" \t "<<o->f(aux);
    }

    aux     =   photonuclear_->get_Continuous()->dEdx();
    result  +=  aux;

    if(DEBUG)
    {
        cerr<<" \t "<<o->f(aux);
    }

    aux     =   epairproduction_->get_Continuous()->dEdx();
    result  +=  aux;

    if(DEBUG)
    {
        cerr<<" \t "<<o->f(aux);
    }

    return -1/result;
}

//----------------------------------------------------------------------------//


double CrossSections::getdx(double ei, double ef, double dist)
{
    if(jt_)
    {
        if(fabs(ei-ef) > fabs(ei)*HALF_PRECISION)
        {
            double aux;

            ini_    =   interpolateJ_->interpolate(ei);
            aux     =   ini_ - interpolateJ_->interpolate(ef);

            if(fabs(aux) > fabs(ini_)*HALF_PRECISION)
            {
                return max(aux, 0.0);
            }

        }

        ini_    =   0;
        return max((interpolateJdf_->interpolate((ei + ef)/2))*(ef - ei), 0.0);

    }
    else
    {
        return integral_->integrateWithLog(ei, ef, this, -dist);
    }

}

//----------------------------------------------------------------------------//

double CrossSections::getef(double ei, double dist)
{

    if(jt_)
    {
        if(ini_ != 0)
        {
            double aux;
            aux     =   interpolateJ_->findLimit(ini_-dist);

            if(fabs(aux) > fabs(ei)*HALF_PRECISION)
            {
                return min( max(aux, particle_->get_low()), ei);
            }
        }

        return min( max(ei+dist/interpolateJdf_->interpolate(ei + dist/(2*interpolateJdf_->interpolate(ei))), particle_->get_low()), ei);
    }
    else
    {
        return integral_->getUpperLimit();
    }

}

//----------------------------------------------------------------------------//


void CrossSections::SetConstants()
{
    set_component(0);

    ci_         =   1.;
    cb_         =   1.;
    cp_         =   1.;
    ce_         =   1.;
    cd_         =   1.;
    cm_         =   1.;
    bremserror  =   0;
    epairerror  =   0;
    photoerror  =   0;
    ionizerror  =   0;
    df_         =   false;
    jt_         =   false;
    lpm_        =   false;
    e_hi        =   PhysicsModel::ebig_;
    g           =   5;

}

//----------------------------------------------------------------------------//

// Getter

Decay* CrossSections::get_decay()
{
    return decay_;
}

Ionizationloss* CrossSections::get_ionization()
{
    return ionization_;
}

Bremsstrahlung* CrossSections::get_bremsstrahlung()
{
    return bremsstrahlung_;
}

Photonuclear* CrossSections::get_photonuclear()
{
    return photonuclear_;
}

Epairproduction* CrossSections::get_epairproduction()
{
    return epairproduction_;
}

//----------------------------------------------------------------------------//

//Setter

void CrossSections::set_decay(Decay *decay)
{
    decay_  =   decay;
}

void CrossSections::set_ionization(Ionizationloss *ionization)
{
    ionization_ =   ionization;
}

void CrossSections::set_bremsstrahlung(Bremsstrahlung *bremsstrahlung)
{
    bremsstrahlung_ =   bremsstrahlung;
}

void CrossSections::set_photonuclear(Photonuclear *photonuclear)
{
    photonuclear_   =   photonuclear;
}

void CrossSections::set_epairproduction(Epairproduction *epairproduction)
{
    epairproduction_    =   epairproduction;
}

void CrossSections::set_component(int component)
{
    component_  =   component;
}

void CrossSections::set_ci(double ci)
{
    ci_ =   ci;
}

void CrossSections::set_cb(double cb)
{
    cb_ =   cb;
}

void CrossSections::set_cp(double cp)
{
    cp_ =   cp;
}

void CrossSections::set_ce(double ce)
{
    ce_ =   ce;
}

void CrossSections::set_cd(double cd)
{
    cd_ =   cd;
}

void CrossSections::set_bremserror(double newbrems)
{
    bremserror  =   newbrems;
}

void CrossSections::set_epairerror(double newepair)
{
    epairerror  =   newepair;
}

void CrossSections::set_photoerror(double newphoto)
{
    photoerror  =   newphoto;
}

void CrossSections::set_ionizerror(double newioniz)
{
    ionizerror  =   newioniz;
}

void CrossSections::set_particle(PROPOSALParticle *particle)
{
    particle_   =   particle;
}

void CrossSections::set_medium(Medium *medium)
{
    medium_ =   medium;
}

void CrossSections::set_integral(Integral *integral)
{
    integral_   =   integral;
}

void CrossSections::set_ini(double ini)
{
    ini_    =   ini;
}

void CrossSections::set_df(bool df)
{
    df_ =   df;
}

void CrossSections::set_interpolateJ(Interpolate *interpolateJ)
{
    interpolateJ_   =   interpolateJ;
}

void CrossSections::set_interpolateJdf(Interpolate *interpolateJdf)
{
    interpolateJdf_ =   interpolateJdf;
}

void CrossSections::set_jt(bool jt)
{
    jt_ =   jt;
}

void CrossSections::set_lpm(bool lpm)
{
    lpm_    =   lpm;
}


