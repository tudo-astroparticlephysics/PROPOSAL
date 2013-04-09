#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Constants.h"
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

Bremsstrahlung::Bremsstrahlung()
    :dndx_integral_()
    ,dndx_interpolant_1d_()
    ,dndx_interpolant_2d_()
    ,prob_for_component_()
{
    dedx_integral_   =  new Integral(IROMB, IMAXS, IPREC);
    dedx_interpolant_=  NULL;

    dndx_integral_.resize(medium_->GetNumCompontents());

    for(int i =0; i<(medium_->GetNumCompontents()); i++)
    {
            dndx_integral_.at(i) =   new Integral(IROMB, IMAXS, IPREC);
    }
    do_dedx_Interpolation_  = false;
    do_dndx_Interpolation_  = false;
}

//----------------------------------------------------------------------------//

Bremsstrahlung::Bremsstrahlung(const Bremsstrahlung &brems)
{
    *this = brems;
}

//----------------------------------------------------------------------------//

Bremsstrahlung& Bremsstrahlung::operator=(const Bremsstrahlung &brems){
    return *this;
}

//----------------------------------------------------------------------------//

Bremsstrahlung::Bremsstrahlung(Particle* particle,
                             Medium* medium,
                             EnergyCutSettings* cut_settings)
    :lorenz_(false)
    ,lorenz_cut_(1.e6)
    ,dndx_integral_()
    ,dndx_interpolant_1d_()
    ,dndx_interpolant_2d_()
    ,eLpm_(0)
    ,prob_for_component_()

{
    particle_                   = particle;
    medium_                     = medium;
    cut_settings_               = cut_settings;
    vMax_                       = 0;
    vUp_                        = 0;
    vMin_                       = 0;
    ebig_                       = BIGENERGY;
    do_dedx_Interpolation_  = false;
    do_dndx_Interpolation_  = false;
    multiplier_                 = 1.;
    parametrization_            = 1;
    lpm_effect_enabled_         = false;
    init_lpm_effect_            = true;
    component_                  = 0;

    dedx_integral_   = new Integral(IROMB, IMAXS, IPREC);
    dedx_interpolant_=  NULL;
    dndx_integral_.resize(medium_->GetNumCompontents());

    for(int i =0; i<(medium_->GetNumCompontents()); i++)
    {
            dndx_integral_.at(i) =   new Integral(IROMB, IMAXS, IPREC);
    }

    prob_for_component_.resize(medium_->GetNumCompontents());

}

//----------------------------------------------------------------------------//


void Bremsstrahlung::SetIntegralLimits(int component){

    component_ = component;

    vMax_   =   1 - (3./4)*SQRTE*(particle_->GetMass()/particle_->GetEnergy())
                *pow((medium_->GetNucCharge().at(component_)) , 1./3);

    if(vMax_<0)
    {
        vMax_   =   0;
    }

    if(lorenz_)
    {
        vMax_   =   min(vMax_, lorenz_cut_/(particle_->GetEnergy()));
    }

    vMax_   =   min(vMax_, (1-(particle_->GetMass()/particle_->GetEnergy())));
    vUp_    =   min(vMax_, cut_settings_->GetCut( particle_->GetEnergy()));

}

//----------------------------------------------------------------------------//

double Bremsstrahlung::CalculatedEdx(){

    double sum  =   0;

    if(multiplier_<=0)
    {
        return 0;
    }

    if(do_dedx_Interpolation_)
    {
        return max(dedx_interpolant_->interpolate(particle_->GetEnergy()), 0.0);
    }

    for(int i=0; i<(medium_->GetNumCompontents()); i++)
    {
        SetIntegralLimits(i);
        sum +=  dedx_integral_->IntegrateOpened(0, vUp_, boost::bind(&Bremsstrahlung::FunctionToDEdxIntegral, this, _1));
    }

    return multiplier_*particle_->GetEnergy()*sum;
}
//----------------------------------------------------------------------------//

double Bremsstrahlung::CalculatedNdx(){

    double sum;
    if(multiplier_<=0)
    {
        return 0;
    }

    sum    =   0;

    for(int i=0; i<(medium_->GetNumCompontents()); i++)
    {

        if(do_dndx_Interpolation_)
        {
            sum    +=  max( dndx_interpolant_1d_.at(i)->interpolate(particle_->GetEnergy()) ,  0.0);
        }
        else
        {
            SetIntegralLimits(i);
            sum    +=  dndx_integral_.at(i)->IntegrateWithLog(vUp_, vMax_, boost::bind(&Bremsstrahlung::FunctionToDNdxIntegral, this, _1));
        }

    }

    return sum;
}
//----------------------------------------------------------------------------//

double Bremsstrahlung::CalculatedNdx(double rnd){

    if(multiplier_<=0)
    {
        return 0;
    }

    if(do_dndx_Interpolation_)
    {
        //rnd_    =   rnd;
    }

    double sum    =   0; //SUM was local we will see if this is necessary

    for(int i=0; i<(medium_->GetNumCompontents()); i++)
    {
        if(do_dndx_Interpolation_)
        {
            prob_for_component_.at(i)   =   max(dndx_interpolant_1d_.at(i)->interpolate(particle_->GetEnergy()) , 0.0);
        }
        else
        {
            SetIntegralLimits(i);
            prob_for_component_.at(i)   =   dndx_integral_.at(i)->IntegrateWithLog(vUp_, vMax_, boost::bind(&Bremsstrahlung::FunctionToDNdxIntegral, this, _1), rnd);
        }
        sum    +=  prob_for_component_.at(i);
    }

    return sum;

}
//----------------------------------------------------------------------------//
double Bremsstrahlung::CalculateStochasticLoss(){
    return 0;
}

double Bremsstrahlung::CalculateStochasticLoss(double rnd1, double rnd2){
    double rand;
    double rsum;

    double rnd_ = rnd1;
    double sum_ = this->CalculatedNdx(rnd1);
    rand    =   rnd2*sum_;
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

                return (particle_->GetEnergy())*(vUp_*exp(dndx_interpolant_2d_.at(i)->findLimit((particle_->GetEnergy()), (rnd_)*prob_for_component_.at(i))*log(vMax_/vUp_)));
            }

            else
            {
                SetIntegralLimits(i);
                return (particle_->GetEnergy())*dndx_integral_.at(i)->GetUpperLimit();
            }

        }
    }

    cout<<"Error (in BremsStochastic/e): sum was not initialized correctly" << endl;
    cout<<"ecut: " << cut_settings_->GetEcut() << "\t vcut: " <<  cut_settings_->GetVcut() << "\t energy: " << particle_->GetEnergy() << "\t type: " << particle_->GetName() << endl;
    return 0;
}
//----------------------------------------------------------------------------//

void Bremsstrahlung::EnableDNdxInterpolation(){
    if(do_dndx_Interpolation_)return;

    double energy = particle_->GetEnergy();
    dndx_interpolant_1d_.resize( medium_->GetNumCompontents() );
    dndx_interpolant_2d_.resize( medium_->GetNumCompontents() );
    for(int i=0; i<(medium_->GetNumCompontents()); i++)
    {
        component_ = i;
        dndx_interpolant_2d_.at(i) = new Interpolant(NUM1, particle_->GetLow(), BIGENERGY,  NUM1, 0, 1, boost::bind(&Bremsstrahlung::FunctionToBuildDNdxInterpolant2D, this, _1 , _2), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false, order_of_interpolation_, true, false, false);
        dndx_interpolant_1d_.at(i) = new Interpolant(NUM1, particle_->GetLow(), BIGENERGY, boost::bind(&Bremsstrahlung::FunctionToBuildDNdxInterpolant, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, true, false, false);

    }
    particle_->SetEnergy(energy);

    do_dndx_Interpolation_=true;
}
//----------------------------------------------------------------------------//

void Bremsstrahlung::EnableDEdxInterpolation()
{
    double energy = particle_->GetEnergy();
    dedx_interpolant_ = new Interpolant(NUM1, particle_->GetLow(), BIGENERGY, boost::bind(&Bremsstrahlung::FunctionToBuildDEdxInterpolant, this, _1), order_of_interpolation_, true, false, true, order_of_interpolation_, false, false, true);
    do_dedx_Interpolation_=true;
    particle_->SetEnergy(energy);
}

//----------------------------------------------------------------------------//


void Bremsstrahlung::DisableDNdxInterpolation()
{

    for(unsigned int i = 0 ; i < dndx_interpolant_1d_.size() ; i++ ){
        delete dndx_interpolant_1d_[i];
    }

    for(unsigned int i = 0 ; i < dndx_interpolant_2d_.size() ; i++ ){
        delete dndx_interpolant_2d_[i];
    }

    dndx_interpolant_1d_.clear();
    dndx_interpolant_2d_.clear();

    do_dndx_Interpolation_  =   false;

}

//----------------------------------------------------------------------------//

void Bremsstrahlung::DisableDEdxInterpolation()
{
    delete dedx_interpolant_;
    do_dedx_Interpolation_  =   false;
}

//----------------------------------------------------------------------------//

double Bremsstrahlung::FunctionToDEdxIntegral(double variable){
    return variable * ElasticBremsstrahlungCrossSection(variable, component_);
}

//----------------------------------------------------------------------------//

double Bremsstrahlung::FunctionToDNdxIntegral(double variable){
    return multiplier_ * ElasticBremsstrahlungCrossSection(variable, component_);

}
//----------------------------------------------------------------------------//

double Bremsstrahlung::KelnerKakoulinPetrukhinParametrization(double v, int i)
{
    double Z3       =   0;
    double result   =   0;
    double Dn       =   0;
    double s1       =   0;

    Z3  =   pow((medium_->GetNucCharge()).at(i), -1./3);

    int step;
    double d, da, dn, Fa, maxV;

    d       =   particle_->GetMass()*particle_->GetMass()
                *v/(2*(particle_->GetEnergy())*(1-v));
    s1      =   (medium_->GetLogConstant()).at(i)*Z3;
    da      =   log(1 + ME/(d*SQRTE*s1));
    Dn      =   1.54*pow((medium_->GetAtomicNum()).at(i), 0.27);
    s1      =   ME*Dn/((particle_->GetMass())*s1);
    dn      =   log(Dn/(1 + d*(Dn*SQRTE - 2)/particle_->GetMass()));
    maxV    =   ME*(particle_->GetEnergy() - particle_->GetMass())
                /((particle_->GetEnergy())
                  *(particle_->GetEnergy() - particle_->GetMomentum() + ME));

    if(v<maxV)
    {
        Fa  =   log(((particle_->GetMass())/d)/(d*(particle_->GetMass())/(ME*ME) + SQRTE)) -
                log(1 + ME/(d*SQRTE*medium_->GetBPrime().at(i)*(pow(medium_->GetNucCharge().at(i) , -2./3))));
    }
    else
    {
        Fa  =   0;
    }

    if((medium_->GetNucCharge()).at(i)==1)
    {
        step    =   0;
    }

    else
    {
        step    =   1;
    }


    result = ((4./3)*(1-v) + v*v)
            *(log((particle_->GetMass())/d)
              - 0.5 -da - dn + (dn*step + Fa)/(medium_->GetNucCharge().at(i)));

    return result;

}

//----------------------------------------------------------------------------//



double Bremsstrahlung::AndreevBezrukovBugaevParametrization(double v, int i)
{

    double aux      =   0;
    double Z3       =   0;
    double result   =   0;

    Z3 = pow((medium_->GetNucCharge()).at(i), -1./3);

    double aux1, aux2, a1, a2,zeta, qc, qmin, x1, x2, d1,d2, psi1, psi2;

    a1      =   111.7*Z3/ME;
    a2      =   724.2*Z3*Z3/ME;
    qc      =   1.9*MMU*Z3;
    aux     =   2*(particle_->GetMass())/qc;
    aux     *=  aux;
    zeta    =   sqrt(1+aux);
    qmin    =   pow((particle_->GetMass()),2)
                *v/((particle_->GetEnergy())*(1-v));

    x1      =   a1*qmin;
    x2      =   a2*qmin;

    if((medium_->GetNucCharge()).at(i)==1)
    {
        d1  =   0;
        d2  =   0;
    }
    else
    {
        aux1    =   log((particle_->GetMass())/qc);
        aux2    =   (zeta/2)*log((zeta+1)/(zeta-1));
        d1      =   aux1 + aux2;
        d2      =   aux1 + ((3 - pow(zeta , 2))*aux2 + aux)/2;
    }

    aux     =   (particle_->GetMass())*a1;
    aux1    =   log(pow(aux , 2)/(1 + pow(x1 , 2)));
    aux     =   (particle_->GetMass())*a2;
    aux2    =   log(pow(aux , 2)/(1 + pow(x2 , 2)));
    psi1    =   (1+ aux1)/2 + (1 + aux2)/(2*(medium_->GetNucCharge()).at(i));
    psi2    =   (2./3 + aux1)/2 +
                (2./3 + aux2)/(2*(medium_->GetNucCharge()).at(i));

    aux1    =   x1*atan(1/x1);
    aux2    =   x2*atan(1/x2);
    psi1    -=  aux1 + aux2/(medium_->GetNucCharge().at(i));
    aux     =   pow(x1 , 2);
    psi2    +=  2*aux*(1 - aux1 + 3./4*log(aux/(1 + aux)));
    aux     =   pow(x2 , 2);
    psi2    +=  2*aux*(1 - aux2 + 3./4*log(aux/(1 + aux)))
                /(medium_->GetNucCharge().at(i));

    psi1    -=  d1;
    psi2    -=  d2;
    result  =   (2-2*v + pow(v , 2))*psi1 - (2./3)*(1-v)*psi2;

    if(result<0)
    {
        result  =   0;
    }

    return result;

}


//----------------------------------------------------------------------------//


double Bremsstrahlung::PetrukhinShestakovParametrization(double v, int i)
{

    double Z3       =   0;
    double result   =   0;
    double d, Fd;

    Z3  =   pow((medium_->GetNucCharge()).at(i), -1./3);

    d   =   pow((particle_->GetMass()) , 2)
            * v/(2*(particle_->GetEnergy())*(1-v));

    Fd  =   189*Z3/ME;
    Fd  =   (particle_->GetMass())*Fd/(1 + SQRTE*d*Fd);

    if((medium_->GetNucCharge()).at(i)>10)
    {
        Fd  *=  (2./3)*Z3;
    }

    result  =   ((4./3)*(1-v) + pow(v , 2))*log(Fd);

    return result;

}

//----------------------------------------------------------------------------//



double Bremsstrahlung::CompleteScreeningCase(double v, int i)
{

    double aux      =   0;
    double Z3       =   0;
    double result   =   0;
    double Lr, fZ, Lp;

    Z3  =   pow((medium_->GetNucCharge()).at(i) , -1./3);

    aux =   ALPHA*(medium_->GetNucCharge().at(i));
    aux *=  aux;
    fZ  =   aux*(1/(1 + aux) + 0.20206 + aux*(-0.0369 + aux*(0.0083 - 0.002*aux)));

    //check rounding
    switch((int)((medium_->GetNucCharge()).at(i) + 0.5))
    {

        case 1:
        {
            Lr  =   5.31;
            Lp  =   6.144;
        }break;

        case 2:
        {
            Lr  =   4.79;
            Lp  =   5.621;
        }break;

        case 3:
        {
            Lr  =   4.74;
            Lp  =   5.805;
        }break;

        case 4:
        {
            Lr  =   4.71;
            Lp  =   5.924;
        }break;

        default:
        {
            Lr  =   log(184.15*Z3);
            Lp  =   log (1194*pow(Z3 , 2));
        }break;

    }

    result = (((4./3)*(1-v) + pow(v , 2))*
              ((medium_->GetNucCharge()).at(i)*(Lr - fZ) + Lp)
             + (1./9)*(1-v)*((medium_->GetNucCharge()).at(i) + 1))
            /(medium_->GetNucCharge()).at(i);

    return result;

}

//----------------------------------------------------------------------------//

double Bremsstrahlung::ElasticBremsstrahlungCrossSection(double v, int i){

    double aux      =   0;
    double Z3       =   0;
    double result   =   0;
    double Dn       =   0;
    double s1       =   0;

    Z3  =   pow((medium_->GetNucCharge()).at(i), -1./3);

    switch(parametrization_)
    {
        case 1:
        {
            result  =   KelnerKakoulinPetrukhinParametrization(v, i);
        }break;

        case 2:
        {
            result  =   AndreevBezrukovBugaevParametrization(v, i);
        }break;

        case 3:
        {
            result  =   PetrukhinShestakovParametrization(v, i);
        }break;

        default:
        {
            result  =   CompleteScreeningCase(v, i);
        }break;

    }

    aux =   2*(medium_->GetNucCharge()).at(i)*(ME/particle_->GetMass())*RE;
    aux *=  (ALPHA/v)*aux*result;

    if(lpm_effect_enabled_)
    {
        if(parametrization_!=1)
        {
            s1  =   (medium_->GetLogConstant()).at(i)*Z3;
            Dn  =   1.54*pow((medium_->GetAtomicNum()).at(i) , 0.27);
            s1  =   ME*Dn/((particle_->GetMass())*s1);
        }
        aux *=  lpm(v,s1);
    }

    double c2   =   pow(particle_->GetCharge() , 2);

    return medium_->GetMolDensity()*medium_->GetAtomInMolecule().at(i)*pow(c2 , 2)*aux;
}

//----------------------------------------------------------------------------//

double Bremsstrahlung::lpm(double v, double s1)
{
    if(init_lpm_effect_)
    {
        Integral* integral_temp = new Integral(IROMB,IMAXS,IPREC);

        lpm_effect_enabled_ = false;
        double sum      =   0;
        double e        =   particle_->GetEnergy();
        init_lpm_effect_    =   false;
        particle_->SetEnergy(BIGENERGY);

        for(int i=0; i < medium_->GetNumCompontents(); i++)
        {

           SetIntegralLimits(i);

           sum +=  integral_temp->IntegrateOpened(0, vUp_, boost::bind(&Bremsstrahlung::FunctionToDEdxIntegral, this, _1));
           sum +=  integral_temp->IntegrateWithLog(vUp_, vMax_, boost::bind(&Bremsstrahlung::FunctionToDEdxIntegral, this, _1));
        }

        eLpm_        =   ALPHA*(particle_->GetMass());
        eLpm_        *=  eLpm_/(4*PI*ME*RE*sum);


        particle_->SetEnergy(e);
        SetIntegralLimits(0);
        lpm_effect_enabled_=true;
        delete integral_temp;

    }

    double G, fi, xi, sp, h, s, s2, s3, ps, Gamma;

    const double fi1    =   1.54954;
    const double G1     =   0.710390;
    const double G2     =   0.904912;
    s1                  *=  s1*SQRT2;
    sp                  =   sqrt(eLpm_*v/(8*(particle_->GetEnergy())*(1-v)));
    h                   =   log(sp)/log(s1);

    if(sp < s1)
    {
        xi  =   2;
    }
    else if(sp < 1)
    {
        xi  =   1 + h - 0.08*(1 - h)*(1 - (1-h)*(1 - h))/log(s1);
    }
    else
    {
        xi  =   1;
    }

    s       =   sp/sqrt(xi);
    Gamma   =   RE*ME/(ALPHA*(particle_->GetMass())*v);
    Gamma   =   1 +4*PI*(medium_->GetMolDensity())*(medium_->GetSumCharge())*RE*pow(Gamma,2);
    s       *=  Gamma;
    s2      =   pow(s,2);
    s3      =   pow(s,3);

    if(s < fi1)
    {
        fi  =   1-exp(-6*s*(1 + (3-PI)*s) + s3/(0.623 + 0.796*s + 0.658*s2));
    }
    else
    {
        fi  =   1 - 0.012/pow(s2 , 2);
    }

    if(s < G1)
    {
        ps  =   1 - exp(-4*s - 8*s2/(1 + 3.936*s + 4.97*s2 - 0.05*s3 + 7.50*pow(s2 , 2)));
        G   =   3*ps - 2*fi;
    }
    else if (s < G2)
    {
        G   =   36*s2/(36*s2 + 1);
    }
    else
    {
        G   =   1 - 0.022/pow(s2 , 2);
    }

    return ((xi/3)*((v*v)*G/(Gamma*Gamma) + 2*(1 + (1-v)*(1-v))*fi/Gamma))/((4./3)*(1-v) + v*v);
}

//----------------------------------------------------------------------------//

double Bremsstrahlung::FunctionToBuildDEdxInterpolant(double energy){
    particle_->SetEnergy(energy);
    return CalculatedEdx();
}


//----------------------------------------------------------------------------//


double Bremsstrahlung::FunctionToBuildDNdxInterpolant(double energy)
{
    return dndx_interpolant_2d_.at(component_)->interpolate(energy, 1.);
}


//----------------------------------------------------------------------------//


double Bremsstrahlung::FunctionToBuildDNdxInterpolant2D(double energy , double v)
{
    particle_->SetEnergy(energy);
    SetIntegralLimits(component_);

    if(vUp_==vMax_)
    {
        return 0;
    }

    v   =   vUp_*exp(v*log(vMax_/vUp_));

    return dndx_integral_.at(component_)->IntegrateWithLog(vUp_, v, boost::bind(&Bremsstrahlung::FunctionToDNdxIntegral, this, _1));
}

//----------------------------------------------------------------------------//

void Bremsstrahlung::SetLorenz(bool lorenz){
    lorenz_ = lorenz;
}
//----------------------------------------------------------------------------//
void Bremsstrahlung::SetLorenzCut(double lorenz_cut){
    lorenz_cut_ = lorenz_cut;
}

//----------------------------------------------------------------------------//


Bremsstrahlung::~Bremsstrahlung()
{
    DisableDNdxInterpolation();
    DisableDEdxInterpolation();

    delete dedx_integral_;
    for(unsigned int i = 0 ; i < dndx_integral_.size() ; i++ ){
        delete dndx_integral_[i];
    }

    dndx_integral_.clear();
}

