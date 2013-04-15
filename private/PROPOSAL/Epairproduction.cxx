#include "PROPOSAL/Epairproduction.h"
#include <algorithm>
#include "boost/bind.hpp"


using namespace std;

Epairproduction::Epairproduction()
    :reverse_(false)
    ,eLpm_(0)
    ,dndx_integral_()
    ,dndx_interpolant_1d_()
    ,dndx_interpolant_2d_()
    ,prob_for_component_()
{
    integral_             = new Integral(IROMB, IMAXS, IPREC);
    integral_for_dEdx_    = new Integral(IROMB, IMAXS, IPREC);
}
//----------------------------------------------------------------------------//

Epairproduction::Epairproduction(const Epairproduction &epair)
{
    *this = epair;
}
//----------------------------------------------------------------------------//

Epairproduction& Epairproduction::operator=(const Epairproduction &epair){
    return *this;
}

//----------------------------------------------------------------------------//

Epairproduction::Epairproduction(Particle* particle,
                             Medium* medium,
                             EnergyCutSettings* cut_settings)
    :reverse_(false)
    ,eLpm_(0)
    ,dndx_integral_()
    ,dndx_interpolant_1d_()
    ,dndx_interpolant_2d_()
    ,prob_for_component_()
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
    parametrization_            = 1;
    lpm_effect_enabled_         = false;
    init_lpm_effect_            = true;
    component_                  = 0;


    integral_             = new Integral(IROMB, IMAXS, IPREC);
    integral_for_dEdx_    = new Integral(IROMB, IMAXS, IPREC);

    dndx_integral_.resize(medium_->GetNumCompontents());


    for(int i =0 ; i<medium_->GetNumCompontents();i++){
        dndx_integral_.at(i) = new Integral(IROMB, IMAXS, IPREC);
    }

    prob_for_component_.resize(medium_->GetNumCompontents());
}

//----------------------------------------------------------------------------//


void Epairproduction::SetIntegralLimits(int component){

    component_ = component;
    double aux;

    vMin_    =   4*ME/particle_->GetEnergy();
    vMax_    =   1 - (3./4)*SQRTE*(particle_->GetMass()/particle_->GetEnergy())
                * pow(medium_->GetNucCharge().at(component) , 1./3);
    aux      =   particle_->GetMass()/particle_->GetEnergy();
    aux      =   1-6*aux*aux;
    vMax_    =   min(vMax_, aux);
    vMax_    =   min(vMax_, 1-particle_->GetMass()/particle_->GetEnergy());

    if(vMax_ < vMin_)
    {
        vMax_    =   vMin_;
    }

    vUp_     =   min(vMax_, cut_settings_->GetCut(particle_->GetEnergy()));

    if(vUp_ < vMin_)
    {
        vUp_     =   vMin_;
    }

}

//----------------------------------------------------------------------------//

double Epairproduction::CalculatedEdx(){
    if(multiplier_<=0)
    {
        return 0;
    }

    if(do_dedx_Interpolation_)
    {
        return max(dedx_interpolant_->interpolate(particle_->GetEnergy()), 0.0);
    }


    double sum  =   0;

    for(int i=0; i<medium_->GetNumCompontents(); i++)
    {
        SetIntegralLimits(i);
        double r1   =   0.8;
        double rUp  =   vUp_*(1-HALF_PRECISION);
        bool rflag  =   false;

        if(r1<rUp)
        {
            if(2*FunctionToDEdxIntegral(r1)<FunctionToDEdxIntegral(rUp))
            {
                rflag   =   true;
            }
        }

        if(rflag)
        {
            if(r1>vUp_)
            {
                r1  =   vUp_;
            }

            if(r1<vMin_)
            {
                r1  =   vMin_;
            }

            sum         +=  integral_for_dEdx_->Integrate(vMin_, r1, boost::bind(&Epairproduction::FunctionToDEdxIntegral, this, _1),4);
            reverse_    =   true;
            double r2   =   max(1-vUp_, COMPUTER_PRECISION);

            if(r2>1-r1)
            {
                r2  =   1-r1;
            }

            sum         +=  integral_for_dEdx_->Integrate(1-vUp_, r2, boost::bind(&Epairproduction::FunctionToDEdxIntegral, this, _1),2)
                        +   integral_for_dEdx_->Integrate(r2, 1-r1, boost::bind(&Epairproduction::FunctionToDEdxIntegral, this, _1),4);

            reverse_    =   false;
        }

        else
        {
            sum +=  integral_for_dEdx_->Integrate(vMin_, vUp_, boost::bind(&Epairproduction::FunctionToDEdxIntegral, this, _1),4);
        }
    }

    return multiplier_*particle_->GetEnergy()*sum;
}
//----------------------------------------------------------------------------//

double Epairproduction::CalculatedNdx(){
    if(multiplier_<=0)
    {
        return 0;
    }


    double sum  =   0;

    for(int i=0; i<medium_->GetNumCompontents(); i++){
        if(do_dndx_Interpolation_)
        {
            sum += max(dndx_interpolant_1d_.at(i)->interpolate(particle_->GetEnergy()), 0.0);
        }
        else
        {
            SetIntegralLimits(i);
            sum += dndx_integral_.at(i)->Integrate(vUp_,vMax_, boost::bind(&Epairproduction::FunctionToDNdxIntegral, this, _1),4);
        }
    }

    return sum;
}
//----------------------------------------------------------------------------//

double Epairproduction::CalculatedNdx(double rnd){
    if(multiplier_<=0)
    {
        return 0;
    }

    double sum  =   0;

    for(int i=0; i<medium_->GetNumCompontents(); i++){
        if(do_dndx_Interpolation_)
        {
            prob_for_component_.at(i) = max(dndx_interpolant_1d_.at(i)->interpolate(particle_->GetEnergy()), 0.0);
        }
        else
        {
            SetIntegralLimits(i);
            prob_for_component_.at(i) = dndx_integral_.at(i)->IntegrateWithLog(vUp_,vMax_, boost::bind(&Epairproduction::FunctionToDNdxIntegral, this, _1),rnd);
        }

        sum += prob_for_component_.at(i);
    }

    return sum;

}
//----------------------------------------------------------------------------//

double Epairproduction::CalculateStochasticLoss(){
    return 0;
}
//----------------------------------------------------------------------------//

double Epairproduction::CalculateStochasticLoss(double rnd1, double rnd2){

    double rand;
    double rsum;

    //double rnd_ = rnd1;
    double sum = this->CalculatedNdx(rnd1);
    rand    =   rnd2*sum;
    rsum    =   0;

    for(int i=0; i<(medium_->GetNumCompontents()); i++)
    {
        rsum    += prob_for_component_.at(i);
        if(rsum > rand)
        {

            if(do_dndx_Interpolation_)
            {
                SetIntegralLimits(i);

                if(vUp_==vMax_)
                {
                    return (particle_->GetEnergy())*vUp_;
                }

                return (particle_->GetEnergy())*(vUp_*exp(dndx_interpolant_2d_.at(i)->findLimit((particle_->GetEnergy()), (rnd1)*prob_for_component_.at(i))*log(vMax_/vUp_)));
            }

            else
            {
                component_ = i;
                return (particle_->GetEnergy())*dndx_integral_.at(i)->GetUpperLimit();
            }
        }
    }

    //TOMASZ sometime everything is fine, just the probability for interaction is zero
    bool prob_for_all_comp_is_zero=true;
    for(int i=0; i<(medium_->GetNumCompontents()); i++)
    {
        SetIntegralLimits(i);
        if(vUp_!=vMax_)prob_for_all_comp_is_zero=false;
    }
    if(prob_for_all_comp_is_zero)return 0;

    cerr<<"Error (in BremsStochastic/e): sum was not initialized correctly" << endl;
    cerr<<"ecut: " << cut_settings_->GetEcut() << "\t vcut: " <<  cut_settings_->GetVcut() << "\t energy: " << particle_->GetEnergy() << "\t type: " << particle_->GetName() << endl;
    return 0;
}

//----------------------------------------------------------------------------//

void Epairproduction::EnableDNdxInterpolation(){
    if(do_dndx_Interpolation_)return;

    double energy = particle_->GetEnergy();
    dndx_interpolant_1d_.resize(medium_->GetNumCompontents());
    dndx_interpolant_2d_.resize(medium_->GetNumCompontents());
    for(int i=0; i<(medium_->GetNumCompontents()); i++)
    {
        component_ = i;
        dndx_interpolant_2d_.at(i) =    new Interpolant(NUM1, particle_->GetLow(), BIGENERGY,  NUM1, 0, 1, boost::bind(&Epairproduction::FunctionToBuildDNdxInterpolant2D, this, _1 , _2), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false, order_of_interpolation_, true, false, false);
        dndx_interpolant_1d_.at(i) =    new Interpolant(NUM1, particle_->GetLow(), BIGENERGY,  boost::bind(&Epairproduction::FunctionToBuildDNdxInterpolant1D, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, true, false, false);
    }
    particle_->SetEnergy(energy);

    do_dndx_Interpolation_=true;
}
//----------------------------------------------------------------------------//

void Epairproduction::EnableDEdxInterpolation(){
    if(do_dedx_Interpolation_)return;

    double energy = particle_->GetEnergy();

    dedx_interpolant_ = new Interpolant(NUM1, particle_->GetLow(), BIGENERGY, boost::bind(&Epairproduction::FunctionToBuildDEdxInterpolant, this, _1), order_of_interpolation_, true, false, true, order_of_interpolation_, false, false, false);

    particle_->SetEnergy(energy);

    do_dedx_Interpolation_=true;
}

//----------------------------------------------------------------------------//


void Epairproduction::DisableDNdxInterpolation(){
    do_dndx_Interpolation_  =   false;
}

//----------------------------------------------------------------------------//

void Epairproduction::DisableDEdxInterpolation(){
    do_dedx_Interpolation_  =   false;
}

//----------------------------------------------------------------------------//

double Epairproduction::FunctionToDEdxIntegral(double variable){

    if(reverse_)
    {
        variable   =   1-variable;
    }

    return variable*EPair(variable, component_);
}

//----------------------------------------------------------------------------//

double Epairproduction::FunctionToDNdxIntegral(double variable){
    return  multiplier_ * EPair(variable, component_);
}
//----------------------------------------------------------------------------//

double Epairproduction::FunctionToIntegral(double r){

    double Fe, Fm, Le, Lm, Ye, Ym, s, b, k, g1, g2;
    double aux, aux1, aux2, r2, Z3;

    r       =   1-r; // only for integral optimization - do not forget to swap integration limits!
    r2      =   r*r;
    Z3      =   pow(medium_->GetNucCharge().at(component_) , -1./3);
    aux     =   (particle_->GetMass()*v_)/(2*ME);
    aux     *=  aux;
    s       =   aux*(1 - r2)/(1 - v_);
    b       =   (v_*v_)/(2*(1 - v_));
    Ye      =   (5 - r2 + 4*b*(1 + r2))/(2*(1 + 3*b)*log(3 + 1/s) - r2 - 2*b*(2 - r2));
    Ym      =   (4 + r2 + 3*b*(1 + r2))/((1 + r2)*(1.5 + 2*b)*log(3 + s) + 1 - 1.5*r2);
    aux     =   (1.5*ME)/(particle_->GetMass()*Z3);
    aux     *=  aux;
    aux1    =   (1 + s)*(1 + Ye);
    aux2    =   (2*ME*SQRTE*medium_->GetLogConstant().at(component_)*Z3) / (particle_->GetEnergy()*v_*(1 - r2));
    Le      =   log((medium_->GetLogConstant().at(component_)*Z3*sqrt(aux1)) / (1 + aux2*aux1)) - 0.5*log(1 + aux*aux1);
    Lm      =   log(((particle_->GetMass()/(1.5*ME))*medium_->GetLogConstant().at(component_)*Z3*Z3)/(1 + aux2*(1 + s)*(1 + Ym)));

    if ( Le > 0 )
    {
        if (1/s < HALF_PRECISION)
        {
            Fe = (1.5 - r2/2 + b*(1 + r2))/s*Le;
        }
        else
        {
            Fe = (((2 + r2)*(1 + b) + s*(3 + r2))*log(1 + 1/s) + (1 - r2 - b)/(1 + s) - (3 + r2))*Le;
        }
    }
    else
    {
        Fe = 0;
    }

    if ( Le > 0)
    {
        Fm = (((1 + r2)*(1 + 1.5*b) - (1 + 2*b)*(1 - r2)/s)*log(1 + s) + s*(1 - r2 - b)/(1 + s) + (1 + 2*b)*(1 - r2))*Lm;
    }

    else
    {
        Fm = 0;
    }

    if(medium_->GetNucCharge().at(component_)==1)
    {
        g1  =   4.4e-5;
        g2  =   4.8e-5;
    }
    else
    {
        g1  =   1.95e-5;
        g2  =   5.3e-5;
    }

    aux     =   particle_->GetEnergy()/particle_->GetMass();
    k       =   0.058*log(aux/(1 + g2*aux/Z3)) - 0.14;

    if(k<=0)
    {
        k   =   0;
    }

    else
    {
        k   =   (0.073*log(aux/(1 + g1*aux/(Z3*Z3))) - 0.26)/k;
    }

    if(k<0)
    {
        k   =   0;
    }

    aux     =   ALPHA*RE;
    aux     *=  aux/(1.5*PI);
    aux1    =   ME/particle_->GetMass();
    aux1    *=  aux1;

    if(lpm_effect_enabled_)
    {
        aux     *=  2*medium_->GetNucCharge().at(component_)
                *  (medium_->GetNucCharge().at(component_) + k)*((1 - v_)/v_)*lpm(r2, b, s)*(Fe + aux1*Fm);
    }
    else
    {
        aux     *=  2*medium_->GetNucCharge().at(component_)
                *  (medium_->GetNucCharge().at(component_) + k)*((1 - v_)/v_)*(Fe + aux1*Fm);
    }
    if(aux<0)
    {
        aux =   0;
    }

    return aux;
}



//----------------------------------------------------------------------------//

double Epairproduction::lpm(double r2, double b, double x)
{


    if(init_lpm_effect_)
    {
        init_lpm_effect_        =   false;
        double sum  =   0;

        for(int i=0; i<medium_->GetNumCompontents(); i++)
        {
            sum +=  medium_->GetNucCharge().at(i)*medium_->GetNucCharge().at(i)
                    *log(3.25*medium_->GetLogConstant().at(i)*pow(medium_->GetNucCharge().at(i), -1./3));
        }

        eLpm_    =   particle_->GetMass()/(ME*RE);
        eLpm_    *=  (eLpm_*eLpm_)*ALPHA*particle_->GetMass()
                /(2*PI*medium_->GetMolDensity()*particle_->GetCharge()*particle_->GetCharge()*sum);
    }

    double A, B, C, D, E, s;
    double s2, s36, s6, d1, d2, atan_, log1, log2;

    s       =   sqrt(eLpm_/(particle_->GetEnergy()*v_*(1 - r2)))/4;
    s6      =   6*s;
    atan_   =   s6*(x + 1);

    if(atan_>1/COMPUTER_PRECISION)
    {
        return 1;
    }

    s2      =   s*s;
    s36     =   36*s2;
    d1      =   s6/(s6 + 1);
    d2      =   s36/(s36 + 1);
    atan_   =   atan(atan_) - PI/2;
    log1    =   log((s36*(1 + x)*(1 + x) + 1)/(s36*x*x));
    log2    =   log((s6*(1 + x) + 1)/(s6*x));
    A       =   0.5*d2*(1 + 2*d2*x)*log1 - d2 + 6*d2*s*(1 + ((s36 - 1)/(s36 + 1))*x)*atan_;
    B       =   d1*(1 + d1*x)*log2 - d1;
    C       =   -d2*d2*x*log1 + d2 - (d2*d2*(s36 - 1)/(6*s))*x*atan_;
    D       =   d1-d1*d1*x*log2;
    E       =   -s6*atan_;

    return ((1 + b)*(A + (1 + r2)*B) + b*(C + (1 + r2)*D) + (1 - r2)*E)
            /(((2 + r2)*(1 +b ) + x*(3 + r2))*log(1 +1 /x) + (1 - r2 - b)/(1 + x) - (3 + r2));


}

//----------------------------------------------------------------------------//

double Epairproduction::EPair(double v, int component)
{

    if(do_dndx_Interpolation_)
    {
        SetIntegralLimits(component);

        if(v>=vUp_)
        {
            return max(dndx_interpolant_2d_.at(component)->interpolate(particle_->GetEnergy(), log(v/vUp_)/log(vMax_/vUp_)), 0.0);
        }
    }

    double rMax, aux, aux2;

    component_  =   component;
    v_          =   v;
    aux         =   1 - (4*ME)/(particle_->GetEnergy()*v_);
    aux2        =   1 - (6*particle_->GetMass()*particle_->GetMass())/(particle_->GetEnergy()*particle_->GetEnergy()*(1 - v_));

    if(aux>0 && aux2>0)
    {
        rMax    =   sqrt(aux)*aux2;
    }
    else
    {
        rMax    =   0;
    }

    aux =   max(1 - rMax , COMPUTER_PRECISION);

    return medium_->GetMolDensity()*medium_->GetAtomInMolecule().at(component_)
           *particle_->GetCharge()*particle_->GetCharge()
           *(integral_->Integrate(1 - rMax, aux, boost::bind(&Epairproduction::FunctionToIntegral, this, _1),2)
                + integral_->Integrate(aux, 1,  boost::bind(&Epairproduction::FunctionToIntegral, this, _1),4));

}

double Epairproduction::FunctionToBuildDNdxInterpolant1D(double energy){
    return dndx_interpolant_2d_.at(component_)->interpolate(energy,1.);
}

double Epairproduction::FunctionToBuildDNdxInterpolant2D(double energy, double v){
    particle_->SetEnergy(energy);
    SetIntegralLimits(component_);

    if(vUp_==vMax_)
    {
    return 0.;
    }

    v   =   vUp_*exp(v*log(vMax_/vUp_));

    return dndx_integral_.at(component_)->Integrate(vUp_,v, boost::bind(&Epairproduction::FunctionToDNdxIntegral, this, _1),4);
}

double Epairproduction::FunctionToBuildDEdxInterpolant(double energy){
    particle_->SetEnergy(energy);
    return CalculatedEdx();
}

Epairproduction::~Epairproduction()
{
    delete integral_for_dEdx_;
    for(unsigned int i = 0 ; i < dndx_integral_.size() ; i++ ){
        delete dndx_integral_[i];
    }

    dndx_integral_.clear();
}
//----------------------------------------------------------------------------//

