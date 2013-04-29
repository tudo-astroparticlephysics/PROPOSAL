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
    do_dndx_Interpolation_      = false;
    multiplier_                 = 1.;


    integral_   = new Integral(IROMB, IMAXS, IPREC);

}

//----------------------------------------------------------------------------//

void Ionization::SetIntegralLimits(int component){

    double aux;

    beta_    =   particle_->GetMomentum()/particle_->GetEnergy();
    gamma_   =   particle_->GetEnergy()/particle_->GetMass();
    vMin_    =   (1.e-6*medium_->GetI())/particle_->GetEnergy();
    aux      =   ME/particle_->GetMass();
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

    if(do_dedx_Interpolation_)
    {
        return max(dedx_interpolant_->interpolate(particle_->GetEnergy()), 0.);
    }

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
                        + particle_->GetEnergy()*(integral_->Integrate(vMin_, vUp_, boost::bind(&Ionization::FunctionToDEdxIntegral, this, _1),4)));
}
//----------------------------------------------------------------------------//

double Ionization::CalculatedNdx(){
    if(multiplier_<=0)
    {
        return 0;
    }

    if(do_dndx_Interpolation_)
    {
        return max(dndx_interpolant_1d_->interpolate(particle_->GetEnergy()), 0.);
    }
    else{
        SetIntegralLimits(0);
        return integral_->Integrate(vUp_,vMax_,boost::bind(&Ionization::FunctionToDNdxIntegral, this, _1),3,1);
    }
}
//----------------------------------------------------------------------------//

double Ionization::CalculatedNdx(double rnd){
    if(multiplier_<=0)
    {
        return 0;
    }

    if(do_dndx_Interpolation_)
    {
        return max(dndx_interpolant_1d_->interpolate(particle_->GetEnergy()), 0.);
    }
    else
    {
        SetIntegralLimits(0);
        return integral_->IntegrateWithSubstitution(vUp_,vMax_,boost::bind(&Ionization::FunctionToDNdxIntegral, this, _1),1,rnd);
    }
}
//----------------------------------------------------------------------------//

double Ionization::CalculateStochasticLoss(){
    return 0;
}

//----------------------------------------------------------------------------//

double Ionization::CalculateStochasticLoss(double rnd1, double rnd2){
    double rand, rsum;

    double sum = this->CalculatedNdx(rnd1);

    rand=medium_->GetSumCharge()*rnd2;
    rsum=0;

    for(int i=0; i<medium_->GetNumCompontents(); i++){
        rsum+=medium_->GetAtomInMolecule().at(i)* medium_->GetNucCharge().at(i);

        if(rsum>rand){

            if(do_dndx_Interpolation_)
            {
                SetIntegralLimits(0);
                if(vUp_==vMax_){
                    return particle_->GetEnergy()*vUp_;
                }
                return particle_->GetEnergy()*(vUp_*exp(dndx_interpolant_2d_->findLimit(particle_->GetEnergy(), rnd1*sum)*log(vMax_/vUp_)));
            }
            else
            {
                return particle_->GetEnergy()*integral_->GetUpperLimit();
            }
        }
    }

    cerr<<"Error (in IonizStochastic/e): m.totZ was not initialized correctly"<<endl;
        return 0;
}

//----------------------------------------------------------------------------//

double Ionization::CalculateStochasticLossNew(double rnd1, double rnd2){
    double rand, rsum;

    double sum = this->CalculatedNdx();

    rand=medium_->GetSumCharge()*rnd2;
    rsum=0;

    Integral* get_upper_integral = new Integral(IROMB, IMAXS, IPREC);


    for(int i=0; i<medium_->GetNumCompontents(); i++){
        rsum+=medium_->GetAtomInMolecule().at(i)* medium_->GetNucCharge().at(i);

        if(rsum>rand){

            if(do_dndx_Interpolation_)
            {
                SetIntegralLimits(0);
                if(vUp_==vMax_){
                    return particle_->GetEnergy()*vUp_;
                }
                return particle_->GetEnergy()*(vUp_*exp(dndx_interpolant_2d_->findLimit(particle_->GetEnergy(), rnd1*sum)*log(vMax_/vUp_)));
            }
            else
            {
                cout<<"GETRANDOMX::::::::::: "<<integral_->GetRandomX()<<endl;
                get_upper_integral->SetRandomX( integral_->GetRandomX());
                SetIntegralLimits(0);
                return particle_->GetEnergy()*get_upper_integral->GetUpperLimit(vUp_,vMax_,sum, rnd1,boost::bind(&Ionization::FunctionToDNdxIntegral, this, _1),3,1);
            }
        }
    }

    cerr<<"Error (in IonizStochastic/e): m.totZ was not initialized correctly"<<endl;
        return 0;
}

//----------------------------------------------------------------------------//


void Ionization::EnableDNdxInterpolation(){
    double energy = particle_->GetEnergy();

    dndx_interpolant_2d_ = new Interpolant(NUM1, particle_->GetLow(), BIGENERGY, NUM1, 0, 1, boost::bind(&Ionization::FunctionToBuildDNdxInterpolant2D, this, _1, _2),order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false, order_of_interpolation_, true, false, false);
    dndx_interpolant_1d_ = new Interpolant(NUM1, particle_->GetLow(), BIGENERGY, boost::bind(&Ionization::FunctionToBuildDNdxInterpolant, this, _1),order_of_interpolation_, false, false, true, order_of_interpolation_, true, false, false);

    do_dndx_Interpolation_=true;
    particle_->SetEnergy(energy);
}
//----------------------------------------------------------------------------//

void Ionization::EnableDEdxInterpolation(){
    double energy = particle_->GetEnergy();
    dedx_interpolant_ = new Interpolant(NUM1, particle_->GetLow(), BIGENERGY, boost::bind(&Ionization::FunctionToBuildDEdxInterpolant, this, _1), order_of_interpolation_, true, false, true, order_of_interpolation_, false, false, true);
    do_dedx_Interpolation_=true;
    particle_->SetEnergy(energy);
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

double Ionization::FunctionToDEdxIntegral(double variable){
    return variable*D2Ndvdx(variable)*InelCorrection(variable);
}

//----------------------------------------------------------------------------//

double Ionization::FunctionToDNdxIntegral(double variable){
    return multiplier_ * D2Ndvdx(variable) * (1+InelCorrection(variable));
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

double Ionization::FunctionToBuildDEdxInterpolant(double energy)
{
    particle_->SetEnergy(energy);
    return CalculatedEdx();
}
//----------------------------------------------------------------------------//

double Ionization::FunctionToBuildDNdxInterpolant(double energy){
    return dndx_interpolant_2d_->interpolate(energy, 1.0);
}



double Ionization::FunctionToBuildDNdxInterpolant2D(double energy , double v){
    particle_->SetEnergy(energy);
    SetIntegralLimits(0);

    if(vUp_==vMax_){
        return 0;
    }
    v=vUp_*exp(v*log(vMax_/vUp_));

    return integral_->Integrate(vUp_, v,boost::bind(&Ionization::FunctionToDNdxIntegral, this, _1),3,1);
}

