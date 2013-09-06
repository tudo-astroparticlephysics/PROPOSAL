#include <cmath>

#include "boost/bind.hpp"

#include "PROPOSAL/Coefficients.h"          //coefficients for calculating the power series approximation of the moliere function
#include "PROPOSAL/ScatteringMoliere.h"

#include "PROPOSAL/Constants.h"


#define E0CGS 4.803204197*1e-10                     //charge of an electron in (cm^(3/2)*g^(1/2))/s (Gauß-cgs)
#define C 0.577215664901532860606512090082402431	//Euler-Mascheroni constant


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

    Z_.resize(numComp_);
    ki_.resize(numComp_);
    Ai_.resize(numComp_);

    for(int i = 0; i < numComp_; i++)
    {
        Z_.at(i)     =   medium_->GetNucCharge().at(i);
        ki_.at(i)    =   medium_->GetAtomInMolecule().at(i);
        Ai_.at(i)    =   medium_->GetAtomicNum().at(i);
    }

    CalcBetaSq();

    CalcWeight();

    CalcChi0();
    CalcChiASq();
    CalcChiCSq();
    CalcB();


    thetaMax_    =  15.*sqrt(0.5*chiCSq_*20.);
    // vaguely 15 sigma of the gaussian aprroximation (for B set to 20)

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

//new constructur as expected in PROPOSAL
ScatteringMoliere::ScatteringMoliere()
{
    MathMachine_ =   new MathModel();
    IntMachine_  =   new Integral();
}

ScatteringMoliere::ScatteringMoliere(const ScatteringMoliere &scattering)
    :dx_(scattering.dx_)
    ,betaSq_(scattering.betaSq_)
    ,p_(scattering.p_)
    ,m_(scattering.m_)
    ,numComp_(scattering.numComp_)
    ,chiCSq_(scattering.chiCSq_)
    ,thetaMax_(scattering.thetaMax_)
{
    Z_       = scattering.Z_;
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

    if(scattering.IntMachine_ != NULL)
    {
        IntMachine_ = new Integral(*scattering.IntMachine_) ;
    }
    else
    {
        IntMachine_ = NULL;
    }
}

/* This is the old constructur. The current version of PROPOSAL does require no transfer parameter for the constructor.
ScatteringMoliere::ScatteringMoliere(double dr, Particle* part, Medium* med)
{
    MathMachine =   new MathModel();
    IntMachine  =   new Integral();

    dx          =   dr;
    medium      =   med;

    numComp     =   medium->GetNumComponents();

    p           =   part->GetMomentum();
    m           =   part->GetMass();

    Z.resize(numComp);
    ki.resize(numComp);
    Ai.resize(numComp);

    for(int i = 0; i < numComp; i++)
    {
        Z.at(i)     =   medium->GetNucCharge().at(i);
        ki.at(i)    =   medium->GetAtomInMolecule().at(i);
        Ai.at(i)    =   medium->GetAtomicNum().at(i);
    }

    CalcBetaSq();

    CalcWeight();

    CalcChi0();
    CalcChiASq();
    CalcChiCSq();
    CalcB();


    phiMax      =   15.*sqrt(0.5*chiCSq*20.);
    // vaguely 15 sigma of the gaussian aprroximation (for B set to 20)
}
*/

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
    if(this->Z_       != scattering.Z_)         return false;

    if(this->dx_         != scattering.dx_ )          return false;
    if(this->betaSq_     != scattering.betaSq_ )      return false;
    if(this->p_          != scattering.p_ )           return false;
    if(this->m_          != scattering.m_ )           return false;
    if(this->numComp_    != scattering.numComp_ )     return false;
    if(this->chiCSq_     != scattering.chiCSq_ )      return false;
    if(this->thetaMax_   != scattering.thetaMax_ )    return false;

    if(*(this->medium_)          != *(scattering.medium_) )       return false;
//    if(*(this->MathMachine)     != *(scattering.MathMachine) )  return false;
    if(*(this->IntMachine_)      != *(scattering.IntMachine_) )   return false;
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

    swap(this->Z_ , scattering.Z_);
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
    swap(this->thetaMax_ , scattering.thetaMax_);

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

    if(scattering.IntMachine_ != NULL)
    {
        IntMachine_->swap(*scattering.IntMachine_) ;
    }
    else
    {
        IntMachine_ = NULL;
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
        chi0_.at(i) = ( ME*ALPHA*pow(Z_.at(i)*128./(9.*PI*PI), 1./3.) )/p_ ;
    }
}

//----------------------------------------------------------------------------//

void ScatteringMoliere::CalcChiASq()
{
    chiASq_.resize(numComp_);

    for(int i = 0; i < numComp_; i++)
    {
        chiASq_.at(i) = chi0_.at(i)*chi0_.at(i)*( 1.13+3.76*Z_.at(i)*Z_.at(i)*ALPHA*ALPHA/betaSq_ );
    }
}

//----------------------------------------------------------------------------//

void ScatteringMoliere::CalcChiCSq()
{
    double y1 = 0.;
    double y2 = 0.;

    for(int i = 0; i < numComp_; i++)
    {
        y1 += weight_.at(i)*Z_.at(i)*(Z_.at(i)+1.);
        y2 += weight_.at(i)*Ai_.at(i);
    }

    chiCSq_ = ( (4.*PI*NA*pow(E0CGS, 4.)/(1e26*E0*E0))*(medium_->GetMassDensity()*medium_->GetRho()*dx_) / (p_*p_*betaSq_) ) * ( y1/y2 );
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
    if(x < -12.) return 0.5*sqrt(PI)/( pow(-x, 1.5)*pow(1.+4.5/x, 2./3.) );

    double sum = coeff1[69];

    // Horner's method
    for(int n = 68; n >= 0; n--) sum = sum*x+coeff1[n];

    return sum;
}

//----------------------------------------------------------------------------//

double ScatteringMoliere::f2M(double x)
{
    // approximation for large numbers to avoid numerical errors
    if(x < -19.) return 0.;

    double sum = coeff2[69];

    for(int n = 68; n >= 0; n--) sum = sum*x+coeff2[n];

    return sum;
}

//----------------------------------------------------------------------------//

double ScatteringMoliere::f(double theta)
{
    double y1 = 0;
    double y2 = 0;

    for(int i = 0; i < numComp_; i++)
    {
        double x = -theta*theta/(chiCSq_*B_.at(i));

        y1 += weight_.at(i)*Z_.at(i)*Z_.at(i)/sqrt( chiCSq_*B_.at(i)*PI )*( exp(x) + f1M(x)/B_.at(i) + f2M(x)/(B_.at(i)*B_.at(i)) );
        y2 += weight_.at(i)*Z_.at(i)*Z_.at(i);
    }

    return y1/y2;
}

//----------------------------------------------------------------------------//
//-------------------------generate random angle------------------------------//
//----------------------------------------------------------------------------//

int ScatteringMoliere::BinarySearch(int n, const double *array, double value)
{
    // Binary search in an array of n values to locate value.
    //
    // Array is supposed  to be sorted prior to this call.
    // If match is found, function returns position of element.
    // If no match found, function gives nearest element smaller than value.

    const double* pind;
    pind = lower_bound(array, array + n, value);

    if ( (pind != array + n) && (*pind == value) )
    {
        return (pind - array);
    }
    else
    {
        return ( pind - array - 1);
    }
}

//----------------------------------------------------------------------------//

double ScatteringMoliere::GetRandom()
{
   // Return a random number following this function shape in [-thetaMax,thetaMax]
   //
   //   The distribution function must be normalized to 1.
   //   For each bin the integral is approximated by a parabola.
   //   The parabola coefficients are stored as non persistent data members
   //   Getting one random number implies:
   //     - Generating a random number between 0 and 1
   //     - Looking in which bin in the normalized integral corresponds to
   //     - Evaluate the parabolic curve in the selected bin to find
   //       the corresponding x value.


    const int NumBin = 100;
    double dtheta = 2.*thetaMax_/NumBin;

    double* integral = new double[NumBin+1];
    double* alpha    = new double[NumBin];
    double* beta     = new double[NumBin];
    double* gamma    = new double[NumBin];

    integral[0] =   0;
    double integ;
    double x0, r1, r2, r3;

    int i;
    for (i = 0; i < NumBin; i++)
    {
        x0      =   -thetaMax_+i*dtheta;

        integ   =   IntMachine_->Integrate(x0, x0+dtheta, boost::bind(&ScatteringMoliere::f, this, _1), 2, 0);

        integral[i+1] = integral[i] + integ;

    //the integral r for each bin is approximated by a parabola
    //  x = alpha + beta*r +gamma*r²
    // compute the coefficients alpha, beta, gamma for each bin

        x0      =   -thetaMax_+i*dtheta;

        r2      =   integ;
        r1      =   IntMachine_->Integrate(x0,x0+0.5*dtheta, boost::bind(&ScatteringMoliere::f, this, _1), 2, 0);
        r3      =   2.*r2 - 4.*r1;

        if (abs(r3) > 1e-8)
        {
            gamma[i] = r3/(dtheta*dtheta);
        }
        else
        {
            gamma[i] = 0;
        }

        beta[i] =   r2/dtheta - gamma[i]*dtheta;
        alpha[i] =  x0;
        gamma[i] *= 2.;
    }


    // return random number
    int nbinmin =   0;
    int nbinmax =   (int)(2.*thetaMax_/dtheta)+2;
    if(nbinmax > NumBin) nbinmax=NumBin;

    double pmin =   integral[nbinmin];
    double pmax =   integral[nbinmax];

    double r, thetaRndm, xx, rr;

    do
    {
        r       =   pmin + (pmax-pmin)*MathMachine_->RandomDouble();

        int bin =   BinarySearch(NumBin, integral, r);

        rr      =   r - integral[bin];

        if(gamma[bin] != 0)
        {
            xx  =   (-beta[bin] + sqrt(beta[bin]*beta[bin]+2.*gamma[bin]*rr))/gamma[bin];
        }
        else
        {
            xx  =   rr/beta[bin];
        }
        thetaRndm =   alpha[bin] + xx;

    } while( (thetaRndm<-thetaMax_) || (thetaRndm > thetaMax_) );


    return thetaRndm;

}
