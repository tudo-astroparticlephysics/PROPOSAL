#include "PROPOSAL/Ionization.h"
#include <algorithm>

using namespace std;

Ionization::Ionization(){

}
//----------------------------------------------------------------------------//

Ionization::Ionization(const Ionization &ioniz)
    :beta_(0)
    ,gamma_(0)
{
    *this = ioniz;
}
//----------------------------------------------------------------------------//

Ionization& Ionization::operator=(const Ionization &ioniz){
    return *this;
}

//----------------------------------------------------------------------------//

Ionization::Ionization(Particle* particle,
                             Medium* medium,
                             EnergyCutSettings* cut_settings)
    :beta_(0)
    ,gamma_(0)
{
    particle_                   = particle;
    medium_                     = medium;
    cut_settings_               = cut_settings;
    vMax_                       = 0;
    vUp_                        = 0;
    vMin_                       = 0;
    ebig_                       = BIGENERGY;
    do_dedx_Interpolation_      = false;
    do_dndx_Interpolation_  = false;
    multiplier_                 = 1.;


    integral_   = new Integral(IROMB, IMAXS, IPREC);

}

//----------------------------------------------------------------------------//

void Ionization::SetIntegralLimits(int component){

    double aux;

    beta_    =   particle_->GetMomentum()/particle_->GetEnergy();
    gamma_   =   particle_->GetEnergy()/particle_->GetMass();
    vMin_    =   (1.e-6*medium_->GetI())/particle_->GetEnergy();
    aux     =   ME/particle_->GetMass();
    vMax_    =   2*ME*(gamma_*gamma_-1)/((1 + 2*gamma_*aux + aux*aux)*particle_->GetEnergy());
    vMax_    =   min(vMax_, 1. - particle_->GetMass()/particle_->GetEnergy());

    if(vMax_<vMin_)
    {
        vMax_    =   vMin_;
    }

    vUp_ =   min(vMax_, cut_settings_->GetCut(particle_->GetEnergy()));

    if(vUp_<vMin_)
    {
        vUp_    =   vMin_;
    }

}

//----------------------------------------------------------------------------//

double Ionization::CalculatedEdx(){

    if(multiplier_<=0)
    {
        return 0;
    }

//    if(jt_)
//    {
//        return max(interpolateJ_->interpolate(particle_->e), 0.);
//    }

    double result, aux;

    SetIntegralLimits(0);

    aux     =   beta_*gamma_/(1.e-6*medium_->GetI());
    result  =   log(vUp_*(2*ME*particle_->GetEnergy()))+2*log(aux);
    aux     =   vUp_/(2*(1 + 1/gamma_));
    result  +=  aux*aux;
    aux     =   beta_*beta_;
    result  -=  aux*(1 + vUp_/vMax_) + Delta();

    if(result>0)
    {
        result*=IONK*particle_->GetCharge()*particle_->GetCharge()*medium_->GetZA()/(2*aux);
    }
    else
    {
        result=0;
    }
    return multiplier_*(medium_->GetMassDensity()*result
                        + particle_->GetEnergy()*(integral_->IntegrateWithLog(vMin_, vUp_, boost::bind(&Ionization::FunctionToContinuousIntegral, this, _1))));
}
//----------------------------------------------------------------------------//

double Ionization::CalculatedNdx(){
    return 0;
}
//----------------------------------------------------------------------------//

double Ionization::CalculatedNdx(double rnd){
    return 0;
}
//----------------------------------------------------------------------------//

double Ionization::CalculateStochasticLoss(){
    return 0;
}
//----------------------------------------------------------------------------//

void Ionization::EnableDNdxInterpolation(){
    do_dndx_Interpolation_=true;
}
//----------------------------------------------------------------------------//

void Ionization::EnableDEdxInterpolation(){
    do_dedx_Interpolation_=true;
}

//----------------------------------------------------------------------------//

void Ionization::DisableDNdxInterpolation(){
    do_dndx_Interpolation_  =   false;
}

//----------------------------------------------------------------------------//

void Ionization::DisableDEdxInterpolation(){
    do_dedx_Interpolation_  =   false;
}


//----------------------------------------------------------------------------//

double Ionization::FunctionToContinuousIntegral(double variable){
    return variable*D2Ndvdx(variable)*InelCorrection(variable);
}

//----------------------------------------------------------------------------//

double Ionization::FunctionToStochasticalIntegral(double variable){
    return 0;
}

//----------------------------------------------------------------------------//

double Ionization::Delta()
{
    double X;

    X   =   log(beta_*gamma_)/log(10);

    if( X < medium_->GetX0())
    {
        return medium_->GetD0()*pow(10 , 2*(X - medium_->GetX0()));
    }
    else if(X < medium_->GetX1())
    {
        return medium_->GetC1() * X + medium_->GetC()
                + medium_->GetA() * pow(medium_->GetX1() - X , medium_->GetM());
    }
    else
    {
        return medium_->GetC1()*X + medium_->GetC();
    }
}

//----------------------------------------------------------------------------//



double Ionization::D2Ndvdx(double v)
{
    double result, aux, aux2;

    aux     =   beta_*beta_;
    aux2    =   v/(1 + 1/gamma_);
    aux2    *=  0.5*aux2;
    result  =   1 - aux*(v/vMax_) + aux2;
    result  *=  IONK*particle_->GetCharge()*particle_->GetCharge()*medium_->GetZA()/(2*aux*particle_->GetEnergy()*v*v);

    return medium_->GetMassDensity()*result;
}

//----------------------------------------------------------------------------//



double Ionization::InelCorrection(double v)
{
    double result, a, b, c;

    a       =   log(1 + 2*v*particle_->GetEnergy()/ME);
    b       =   log((1 - v/vMax_)/(1 - v));
    c       =   log((2*gamma_*(1 - v)*ME)/(particle_->GetMass()*v));
    result  =   a*(2*b + c) - b*b;

    return (ALPHA/(2*PI))*result;
}
//----------------------------------------------------------------------------//
