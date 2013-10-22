#include <cmath>

#include <boost/math/special_functions/erf.hpp>

#include "PROPOSAL/Coefficients.h"                  //coefficients for calculating the power series approximation of the moliere function
#include "PROPOSAL/Constants.h"

#include "PROPOSAL/ScatteringMoliere.h"


#define HBAR    6.58211928e-22                          //hbar in MeV*s
#define C       0.577215664901532860606512090082402431	//Euler-Mascheroni constant

#define erf(x)      boost::math::erf(x)
#define erfInv(x)   boost::math::erf_inv(x)

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ScatteringMoliere::Scatter(double dr, Particle* part, Medium* med)
{
    double rnd1, rnd2, sx, tx, sy, ty, sz, tz, ax, ay, az;
    double x, y, z;

    dx_          =   dr;
    medium_      =   med;

    numComp_     =   medium_->GetNumComponents();

    p_           =   part->GetMomentum();
    m_           =   part->GetMass();

    Zi_.resize(numComp_);
    ki_.resize(numComp_);
    Ai_.resize(numComp_);

    for(int i = 0; i < numComp_; i++)
    {
        Zi_.at(i)    =   medium_->GetNucCharge().at(i);
        ki_.at(i)    =   medium_->GetAtomInMolecule().at(i);
        Ai_.at(i)    =   medium_->GetAtomicNum().at(i);
    }

    CalcBetaSq();

    CalcWeight();

    CalcChi0();
    CalcChiASq();
    CalcChiCSq();
    CalcB();

    //----------------------------------------------------------------------------//

    rnd1    =   GetRandom();
    rnd2    =   GetRandom();
    sx      =   (rnd1/SQRT3+rnd2)/2;
    tx      =   rnd2;

    rnd1    =   GetRandom();
    rnd2    =   GetRandom();

    sy      =   (rnd1/SQRT3+rnd2)/2;
    ty      =   rnd2;

    sz      =   sqrt(max(1.-(sx*sx+sy*sy), 0.));
    tz      =   sqrt(max(1.-(tx*tx+ty*ty), 0.));

    double sinth, costh,sinph,cosph;
    double theta, phi;

    sinth = part->GetSinTheta();
    costh = part->GetCosTheta();
    sinph = part->GetSinPhi();
    cosph = part->GetCosPhi();

    x   = part->GetX();
    y   = part->GetY();
    z   = part->GetZ();


    ax      =   sinth*cosph*sz+costh*cosph*sx-sinph*sy;
    ay      =   sinth*sinph*sz+costh*sinph*sx+cosph*sy;
    az      =   costh*sz-sinth*sx;

    x       +=  ax*dr;
    y       +=  ay*dr;
    z       +=  az*dr;


    ax      =   sinth*cosph*tz+costh*cosph*tx-sinph*ty;
    ay      =   sinth*sinph*tz+costh*sinph*tx+cosph*ty;
    az      =   costh*tz-sinth*tx;



    costh   =   az;
    sinth   =   sqrt(max(1.-costh*costh, 0.));

    if(sinth!=0.)
    {
        sinph   =   ay/sinth;
        cosph   =   ax/sinth;
    }

    if(costh>1.)
    {
        theta   =   acos(1.)*180./PI;
    }
    else if(costh<-1.)
    {
        theta   =   acos(-1.)*180./PI;
    }
    else
    {
        theta   =   acos(costh)*180./PI;
    }

    if(cosph>1)
    {
        phi =   acos(1.)*180./PI;
    }
    else if(cosph<-1.)
    {
        phi =   acos(-1)*180./PI;
    }
    else
    {
        phi =   acos(cosph)*180./PI;
    }

    if(sinph<0.)
    {
        phi =   360.-phi;
    }

    if(phi>=360.)
    {
        phi -=  360.;
    }

    part->SetX(x);
    part->SetY(y);
    part->SetZ(z);
    part->SetPhi(phi);
    part->SetTheta(theta);
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

ScatteringMoliere::ScatteringMoliere()
{
    MathMachine_ =   new MathModel();
}

ScatteringMoliere::ScatteringMoliere(const ScatteringMoliere &scattering)
    :dx_(scattering.dx_)
    ,betaSq_(scattering.betaSq_)
    ,p_(scattering.p_)
    ,m_(scattering.m_)
    ,numComp_(scattering.numComp_)
    ,chiCSq_(scattering.chiCSq_)
{
    Zi_      = scattering.Zi_;
    ki_      = scattering.ki_;
    Ai_      = scattering.Ai_;
    weight_  = scattering.weight_;

    chi0_    = scattering.chi0_;
    chiASq_  = scattering.chiASq_;
    B_       = scattering.B_;

    if(scattering.medium_ != NULL)
    {
        medium_ = new Medium(*scattering.medium_) ;
    }
    else
    {
        medium_ = NULL;
    }

    if(scattering.MathMachine_ != NULL)
    {
//        MathMachine->swap(*scattering.MathMachine) ;
    }
    else
    {
        MathMachine_ = NULL;
    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


ScatteringMoliere& ScatteringMoliere::operator=(const ScatteringMoliere &scattering){
    if (this != &scattering)
    {
      ScatteringMoliere tmp(scattering);
      swap(tmp);
    }
    return *this;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool ScatteringMoliere::operator==(const ScatteringMoliere &scattering) const
{
    if(this->Ai_      != scattering.Ai_)        return false;
    if(this->B_       != scattering.B_)         return false;
    if(this->chi0_    != scattering.chi0_)      return false;
    if(this->chiASq_  != scattering.chiASq_)    return false;
    if(this->ki_      != scattering.ki_)        return false;
    if(this->weight_  != scattering.weight_)    return false;
    if(this->Zi_      != scattering.Zi_)         return false;

    if(this->dx_         != scattering.dx_ )          return false;
    if(this->betaSq_     != scattering.betaSq_ )      return false;
    if(this->p_          != scattering.p_ )           return false;
    if(this->m_          != scattering.m_ )           return false;
    if(this->numComp_    != scattering.numComp_ )     return false;
    if(this->chiCSq_     != scattering.chiCSq_ )      return false;

    if(*(this->medium_)          != *(scattering.medium_) )       return false;
//    if(*(this->MathMachine)     != *(scattering.MathMachine) )  return false;
    return true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool ScatteringMoliere::operator!=(const ScatteringMoliere &scattering) const {
  return !(*this == scattering);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ScatteringMoliere::swap(ScatteringMoliere &scattering)
{
    using std::swap;

    swap(this->Zi_ , scattering.Zi_);
    swap(this->ki_ , scattering.ki_);
    swap(this->Ai_ , scattering.Ai_);
    swap(this->weight_ , scattering.weight_);
    swap(this->chi0_ , scattering.chi0_);
    swap(this->chiASq_ , scattering.chiASq_);
    swap(this->B_ , scattering.B_);

    swap(this->dx_ , scattering.dx_);
    swap(this->betaSq_ , scattering.betaSq_);
    swap(this->p_ , scattering.p_);
    swap(this->m_ , scattering.m_);
    swap(this->numComp_ , scattering.numComp_);
    swap(this->chiCSq_ , scattering.chiCSq_);

    if(scattering.medium_ != NULL)
    {
        medium_->swap(*scattering.medium_) ;
    }
    else
    {
        medium_ = NULL;
    }

    if(scattering.MathMachine_ != NULL)
    {
//        MathMachine->swap(*scattering.MathMachine) ;
    }
    else
    {
        MathMachine_ = NULL;
    }
}



//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------private member functions---------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//---------------------------calculate parameters-----------------------------//
//----------------------------------------------------------------------------//

void ScatteringMoliere::CalcBetaSq()
{
    betaSq_ = 1./( 1.+m_*m_/(p_*p_) );
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

void ScatteringMoliere::CalcWeight()
{
    weight_.resize(numComp_);

    double A = 0.;

    for(int i = 0; i < numComp_; i++)
    {
        A += ki_.at(i)*Ai_.at(i);
    }

    for(int i = 0; i < numComp_; i++)
    {
        weight_.at(i) = ki_.at(i)*Ai_.at(i)/A;
    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

void ScatteringMoliere::CalcChi0()
{
    chi0_.resize(numComp_);

    for(int i = 0; i < numComp_; i++)
    {
        chi0_.at(i) = ( ME*ALPHA*pow(Zi_.at(i)*128./(9.*PI*PI), 1./3.) )/p_ ;
    }
}

//----------------------------------------------------------------------------//

void ScatteringMoliere::CalcChiASq()
{
    chiASq_.resize(numComp_);

    for(int i = 0; i < numComp_; i++)
    {
        chiASq_.at(i) = chi0_.at(i)*chi0_.at(i)*( 1.13+3.76*ALPHA*ALPHA*Zi_.at(i)*Zi_.at(i)/betaSq_ );
    }
}

//----------------------------------------------------------------------------//

void ScatteringMoliere::CalcChiCSq()
{
    double y1 = 0.;
    double y2 = 0.;

    for(int i = 0; i < numComp_; i++)
    {
        //if case of an electron, replace Z² by Z(Z+1) to into account scatterings
        //on atomic electrons in the medium
        if(m_ == ME) y1 += weight_.at(i)*Zi_.at(i)*(Zi_.at(i)+1.);
        else y1 += weight_.at(i)*Zi_.at(i)*Zi_.at(i);
        y2 += weight_.at(i)*Ai_.at(i);
    }

    chiCSq_ = ( (4.*PI*NA*ALPHA*ALPHA*HBAR*HBAR*SPEED*SPEED)*(medium_->GetMassDensity()*medium_->GetRho()*dx_) / (p_*p_*betaSq_) ) * ( y1/y2 );
}


void ScatteringMoliere::CalcB()
{
    B_.resize(numComp_);

    for(int i = 0; i < numComp_; i++)
    {
        //calculate B-ln(B) = ln(chi_c^2/chi_a^2)+1-2*C via Newton-Raphson method
        double xn = 15.;

        for(int n = 0; n < 6; n++)
        {
            xn = xn*( (1.-log(xn)-log(chiCSq_/chiASq_.at(i))-1.+2.*C)/(1.-xn) );
        }

        B_.at(i) = xn;
    }
}

//----------------------------------------------------------------------------//
//--------------------------calculate distribution----------------------------//
//----------------------------------------------------------------------------//

double ScatteringMoliere::f1M(double x)
{
    // approximation for large numbers to avoid numerical errors
    if(x > 12.) return 0.5*sqrt(PI)/( pow(x, 1.5)*pow(1.-4.5/x, 2./3.) );

    double sum = c1[69];

    // Horner's method
    for(int p = 68; p >= 0; p--) sum = sum*x+c1[p];

    return sum;
}

//----------------------------------------------------------------------------//

double f2Mlarge(double x)
{
    double a = 0.00013567765224589459194769192063035;
    double b = -0.0022635502525409950842771866774683;
    double c = 0.0098037758070269476889935233998585;

    // the junction of both parametrizations is smoothed by an interpolation parabola
    if((x >= 4.25*4.25) && (x <= 6.5*6.5)) return a*x+b*sqrt(x)+c;

    double sum = 0;

    for(int p = 2; p < 13; p++)
    {
        sum += c2large[p]*(0.5*log(x)+s2large[p])*pow(x, -(p+0.5));
    }

    return sum;
}

double ScatteringMoliere::f2M(double x)
{
    // approximation for larger x to avoid numerical errors
    if(x > 4.25*4.25) return f2Mlarge(x);

    double sum = c2[69];

    for(int p = 68; p >= 0; p--) sum = sum*x+c2[p];

    return sum;
}

//----------------------------------------------------------------------------//

double ScatteringMoliere::f(double theta)
{
    double y1 = 0;
    double y2 = 0;

    for(int i = 0; i < numComp_; i++)
    {
        double x = theta*theta/(chiCSq_*B_.at(i));

        y1 += weight_.at(i)*Zi_.at(i)*Zi_.at(i)/sqrt( chiCSq_*B_.at(i)*PI )*( exp(-x) + f1M(x)/B_.at(i) + f2M(x)/(B_.at(i)*B_.at(i)) );
        y2 += weight_.at(i)*Zi_.at(i)*Zi_.at(i);
    }

    return y1/y2;
}

//----------------------------------------------------------------------------//
//----------------------calculate indefinite integral-------------------------//
//----------------------------------------------------------------------------//

double F1Mlarge(double x)
{
    double sum = C1large[14];

    // Horner's method
    for(int p = 13; p >= 0; p--) sum = C1large[p] + sum/x;

    return sum;
}

double ScatteringMoliere::F1M(double x)
{
    if(x > 12.) return F1Mlarge(x);

    double sum = c1[69]/(2.*69+1.);

    // Horner's method
    for(int p = 68; p >= 0; p--) sum = sum*x+c1[p]/(2.*p+1.);

    return sum*sqrt(x);
}

//----------------------------------------------------------------------------//

double F2Mlarge(double x)
{
    double a = -0.00026360133958801203364619158975302;
    double b = 0.0039965027465457608410459577896745;
    double c = -0.016305842044996649714549974419242;

    // the junction of both parametrizations is smoothed by an interpolation parabola
    if((x >= 4.25*4.25) && (x <= 6.5*6.5)) return a*x+b*sqrt(x)+c;

    double sum = 0;

    for(int p = 2; p < 13; p++)
    {
        sum += -0.5*c2large[p]/p*(0.5/p+0.5*log(x)+s2large[p])*pow(x, -p);
    }

    return sum;
}

double ScatteringMoliere::F2M(double x)
{
    if(x > 4.25*4.25) return F2Mlarge(x);

    double sum = c2[69]/(2.*69+1.);

    for(int p = 68; p >= 0; p--) sum = sum*x+c2[p]/(2.*p+1.);

    return sum*sqrt(x);
}

//----------------------------------------------------------------------------//

double ScatteringMoliere::F(double theta)
{
    double y1 = 0;
    double y2 = 0;

    for(int i = 0; i < numComp_; i++)
    {
        double x = theta*theta/(chiCSq_*B_.at(i));

        y1 += weight_.at(i)*Zi_.at(i)*Zi_.at(i)*( 0.5*erf(sqrt(x)) + sqrt(1./PI)*( F1M(x)/B_.at(i) + F2M(x)/(B_.at(i)*B_.at(i)) ) );
        y2 += weight_.at(i)*Zi_.at(i)*Zi_.at(i);
    }

    return (theta < 0.) ? (-1.)*y1/y2 : y1/y2;
}

//----------------------------------------------------------------------------//
//-------------------------generate random angle------------------------------//
//----------------------------------------------------------------------------//


double ScatteringMoliere::GetRandom()
{
    //  Generate random angles following Moliere's distribution by comparing a
    //  uniformly distributed random number with the integral of the distribution.
    //  Therefore, determine the angle where the integral is equal to the random number.


    //rndm element of ]-0.5,0.5]
    double rndm = (MathMachine_->RandomDouble()-0.5);


    // Newton-Raphson method:
    double theta_n;

        // guessing an initial value by assuming a gaussian distribution
        // only regarding the component j with maximum weight

    int j = 0;
    for(int i = 0; i+1 < numComp_; i++)
    {
        if( weight_.at(i+1) > weight_.at(i) ) j = i+1;
    }
    double theta_np1 = sqrt( chiCSq_*B_.at(j) )*erfInv(2.*rndm);

        // iterating until the number of correct digits is greater than 4

    do
    {
        theta_n = theta_np1;
        theta_np1 = theta_n - (F(theta_n)-rndm)/f(theta_n);

    } while( abs((theta_n-theta_np1)/theta_np1) > 1e-4 );

    return theta_np1;
}
