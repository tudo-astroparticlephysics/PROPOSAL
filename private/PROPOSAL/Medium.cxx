/*! \file   Medium.cxx
*   \brief  Source file for the medium routines.
*
*   For more details see the class documentation.
*
*   \date   2013.03.14
*   \author Jan-Hendrik Koehne
*/

#include <cmath>
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Integral.h"
#include "boost/bind.hpp"
#include <iomanip>
#include "PROPOSAL/Output.h"


using namespace std;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------Temp. Func for radiation length calc.--------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double fZ(unsigned int Z);
double Lrad(unsigned int Z);
double Lrad_dash(unsigned int Z);
double X0_inv(unsigned int Z, double M);


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Medium::Medium()
    :numCompontents_    (0)
    ,nucCharge_         ()
    ,atomicNum_         ()
    ,atomInMolecule_    ()
    ,sumCharge_         (0)
    ,ZA_                (0)
    ,I_                 (0)
    ,C1_                (0)
    ,C_                 (0)
    ,a_                 (0)
    ,m_                 (0)
    ,X0_                (0)
    ,X1_                (0)
    ,d0_                (0)
    ,r_                 (0)
    ,logConstant_       ()
    ,bPrime_            ()
    ,rho_               (1.)
    ,massDensity_       (0)
    ,molDensity_        (0)
    ,radiationLength_   (0)
    ,M_                 ()
    ,elementName_       ()
    ,name_              ("water")
    ,ecut_              (0)
    ,vcut_              (0)
    ,vCut_              (0)
    ,mN_                ()
    ,MM_                (0)
    ,sumNucleons_       (0)
    ,r0_                (0)
{
    log_warn("Standard constructor called: defaulting to water");
    InitWater();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Medium::Medium(const Medium &medium)
    :numCompontents_    (medium.numCompontents_)
    ,nucCharge_         (medium.nucCharge_ )
    ,atomicNum_         (medium.atomicNum_ )
    ,atomInMolecule_    (medium.atomInMolecule_ )
    ,sumCharge_         (medium.sumCharge_)
    ,ZA_                (medium.ZA_)
    ,I_                 (medium.I_)
    ,C1_                (medium.C1_)
    ,C_                 (medium.C_)
    ,a_                 (medium.a_)
    ,m_                 (medium.m_)
    ,X0_                (medium.X0_)
    ,X1_                (medium.X1_)
    ,d0_                (medium.d0_)
    ,r_                 (medium.r_)
    ,logConstant_       (medium.logConstant_ )
    ,bPrime_            (medium.bPrime_ )
    ,rho_               (medium.rho_)
    ,massDensity_       (medium.massDensity_)
    ,molDensity_        (medium.molDensity_)
    ,radiationLength_    (medium.radiationLength_)
    ,M_                 (medium.M_ )
    ,elementName_       (medium.elementName_ )
    ,name_              (medium.name_ )
    ,ecut_              (medium.ecut_)
    ,vcut_              (medium.vcut_)
    ,vCut_              (medium.vCut_)
    ,mN_                (medium.mN_ )
    ,MM_                (medium.MM_)
    ,sumNucleons_       (medium.sumNucleons_)
    ,r0_                (medium.r0_)
{

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Medium::Medium(string w, double rho)
    :numCompontents_    (0)
    ,nucCharge_         ()
    ,atomicNum_         ()
    ,atomInMolecule_    ()
    ,sumCharge_         (0)
    ,ZA_                (0)
    ,I_                 (0)
    ,C1_                (0)
    ,C_                 (0)
    ,a_                 (0)
    ,m_                 (0)
    ,X0_                (0)
    ,X1_                (0)
    ,d0_                (0)
    ,r_                 (0)
    ,logConstant_       ()
    ,bPrime_            ()
    ,massDensity_       (0)
    ,molDensity_        (0)
    ,radiationLength_   (0)
    ,M_                 ()
    ,elementName_       ()
    ,ecut_              (0)
    ,vcut_              (0)
    ,vCut_              (0)
    ,mN_                ()
    ,MM_                (0)
    ,sumNucleons_       (0)
    ,r0_                (0)
{
    name_   =   w;

    if (rho > 0 )
    {
        rho_  =   rho;
    }

    else
    {
        rho_  =   1;
    }

    if(EqualsIgnoreCase(w,"water"))
    {
        InitWater();
    }
    else if(EqualsIgnoreCase(w,"ice"))
    {
        InitIce();
    }
    else if(EqualsIgnoreCase(w,"salt"))
    {
        InitSalt();
    }
    else if(EqualsIgnoreCase(w,"standard_rock"))
    {
        InitStandardrock();
    }
    else if(EqualsIgnoreCase(w,"frejus_rock"))
    {
        InitFrejusrock();
    }
    else if(EqualsIgnoreCase(w,"iron"))
    {
        InitIron();
    }
    else if(EqualsIgnoreCase(w,"hydrogen"))
    {
        InitHydrogen();
    }
    else if(EqualsIgnoreCase(w,"lead"))
    {
        InitLead();
    }
    else if(EqualsIgnoreCase(w,"copper"))
    {
        InitCopper();
    }
    else if(EqualsIgnoreCase(w,"uranium"))
    {
        InitUranium();
    }
    else if(EqualsIgnoreCase(w,"air"))
    {
        InitAir();
    }
    else if(EqualsIgnoreCase(w,"mineral_oil"))
    {
        InitParaffin();
    }
    else if(EqualsIgnoreCase(w,"antares_water"))
    {
        InitAntaresWater();
    }
    else
    {
        printf("Warning (in Medium/Medium): unknown medium: defaulting to water");
        name_   =   "water";
        InitWater();
    }

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Medium& Medium::operator=(const Medium &medium)
{
    if (this != &medium)
    {
      Medium tmp(medium);
      swap(tmp);
    }
    return *this;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Medium::operator==(const Medium &medium) const
{
    if( numCompontents_   != medium.numCompontents_)return false;
    if( sumCharge_        != medium.sumCharge_)     return false;
    if( ZA_               != medium.ZA_)            return false;
    if( I_                != medium.I_)             return false;
    if( C1_               != medium.C1_)            return false;
    if( C_                != medium.C_)             return false;
    if( a_                != medium.a_)             return false;
    if( m_                != medium.m_)             return false;
    if( X0_               != medium.X0_)            return false;
    if( X1_               != medium.X1_)            return false;
    if( d0_               != medium.d0_)            return false;
    if( r_                != medium.r_)             return false;
    if( rho_              != medium.rho_)           return false;
    if( massDensity_      != medium.massDensity_)   return false;
    if( molDensity_       != medium.molDensity_)    return false;
    if( radiationLength_  != medium.radiationLength_)return false;
    if( ecut_             != medium.ecut_)          return false;
    if( vcut_             != medium.vcut_)          return false;
    if( vCut_             != medium.vCut_)          return false;
    if( MM_               != medium.MM_)            return false;
    if( sumNucleons_      != medium.sumNucleons_)   return false;
    if( r0_               != medium.r0_)            return false;

    if( M_.size()               != medium.M_.size() )               return false;
    if( logConstant_.size()     != medium.logConstant_.size() )     return false;
    if( bPrime_.size()          != medium.bPrime_.size() )          return false;
    if( mN_.size()              != medium.mN_.size() )              return false;
    if( nucCharge_.size()       != medium.nucCharge_.size() )       return false;
    if( atomicNum_.size()       != medium.atomicNum_.size() )       return false;
    if( atomInMolecule_.size()  != medium.atomInMolecule_.size() )  return false;
    if( elementName_.size()     != medium.elementName_.size() )     return false;

    if( name_.compare( medium.name_) != 0 ) return false;
    for(unsigned int i =    0; i < M_.size(); i++)
    {
        if(M_.at(i)             !=  medium.M_.at(i))            return false;
    }

    for(unsigned int i =    0; i < logConstant_.size(); i++)
    {
        if(logConstant_.at(i)   !=  medium.logConstant_.at(i))  return false;
    }

    for(unsigned int i =    0; i < bPrime_.size(); i++)
    {
        if(bPrime_.at(i)        !=  medium.bPrime_.at(i))       return false;
    }

    for(unsigned int i =    0; i < mN_.size(); i++)
    {
        if(mN_.at(i)            !=  medium.mN_.at(i))           return false;
    }

    for(unsigned int i =    0; i < nucCharge_.size(); i++)
    {
        if(nucCharge_.at(i)     !=  medium.nucCharge_.at(i))    return false;
    }

    for(unsigned int i =    0; i < atomicNum_.size(); i++)
    {
        if(atomicNum_.at(i)     !=  medium.atomicNum_.at(i))    return false;
    }

    for(unsigned int i =    0; i < atomInMolecule_.size(); i++)
    {
        if(atomInMolecule_.at(i)!=  medium.atomInMolecule_.at(i))   return false;
    }

    for(unsigned int i =    0; i < elementName_.size(); i++)
    {
        if(elementName_.at(i).compare( medium.elementName_.at(i))!=0)return false;
    }
    //else
    return true;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Medium::operator!=(const Medium &medium) const
{
  return !(*this == medium);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


ostream& operator<<(ostream& os, Medium const& medium)
{

    os<<"--------------------Medium( "<<&medium<<" )--------------------"<<endl;
    os<<medium.name_<<endl;
    os<<"\tnumber of components:\t\t\t\t"<<medium.numCompontents_<<endl;
    os<<"\tmass density [g/cm3]:\t\t\t\t"<<medium.massDensity_<<endl;
    os<<"\tmolecule density [number/cm3]:\t\t\t"<<medium.molDensity_<<endl;
    os<<"\t<Z/A>:\t\t\t\t\t\t"<<medium.ZA_<<endl;
    os<<"\tsum of nucleons of all nuclei:\t\t\t"<<medium.sumNucleons_<<endl;
    os<<"\tionization potential [eV]:\t\t\t"<<medium.I_<<endl;
    os<<"\trefraction index:\t\t\t\t"<<medium.r_<<endl;
    os<<"\taverage all-component nucleon weight [MeV]:\t"<<medium.MM_<<endl;
    os<<"\tmultiplicative density correction factor:\t"<<medium.rho_<<endl;
    os<<"\tsum of charges of all nuclei:\t\t\t"<<medium.sumCharge_<<endl;
    os<<"\tradiation Length:\t\t"<<medium.radiationLength_<<endl;

    os<<"\n# | Name | atomic number | number of atoms in a molecule | nucleus charge | average nucleon weight in a nucleus [MeV]\n"<<endl;

    for(int i = 0; i<medium.numCompontents_; i++)
    {
        os<<fixed<<std::setprecision(6)<<"\t"<<i <<"\t"<<medium.elementName_.at(i)<<"\t"<<medium.atomicNum_.at(i)<<"\t"
          <<medium.atomInMolecule_.at(i)<<"\t"<<medium.nucCharge_.at(i)<<"\t"
          <<medium.M_.at(i)<<endl;
    }

    os<<"------------------------------------------------------------";
    return os;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Medium::swap(Medium &medium)
{
    using std::swap;

    swap( numCompontents_   , medium.numCompontents_);
    swap( sumCharge_        , medium.sumCharge_);
    swap( ZA_               , medium.ZA_);
    swap( I_                , medium.I_);
    swap( C1_               , medium.C1_);
    swap( C_                , medium.C_);
    swap( a_                , medium.a_);
    swap( m_                , medium.m_);
    swap( X0_               , medium.X0_);
    swap( X1_               , medium.X1_);
    swap( d0_               , medium.d0_);
    swap( r_                , medium.r_);
    swap( rho_              , medium.rho_);
    swap( massDensity_      , medium.massDensity_);
    swap( molDensity_       , medium.molDensity_);
    swap( radiationLength_  , medium.radiationLength_);
    swap( ecut_             , medium.ecut_);
    swap( vcut_             , medium.vcut_);
    swap( vCut_             , medium.vCut_);
    swap( MM_               , medium.MM_);
    swap( sumNucleons_      , medium.sumNucleons_);
    swap( r0_               , medium.r0_);

    M_.swap( medium.M_ );
    elementName_.swap( medium.elementName_ );
    logConstant_.swap( medium.logConstant_ );
    bPrime_.swap( medium.bPrime_ );
    mN_.swap( medium.mN_ );
    nucCharge_.swap( medium.nucCharge_ );
    atomicNum_.swap( medium.atomicNum_ );
    atomInMolecule_.swap( medium.atomInMolecule_ );
    name_.swap( medium.name_ );

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------private member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Medium::Inita(int i)
{
    numCompontents_ =   i;

    nucCharge_.resize(numCompontents_);
    atomicNum_.resize(numCompontents_);
    atomInMolecule_.resize(numCompontents_);
    logConstant_.resize(numCompontents_);
    bPrime_.resize(numCompontents_);
    M_.resize(numCompontents_);
    elementName_.resize(numCompontents_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Medium::Initr()
{
    int i;
    bool flag       =   false;
    double aux1     =   0;
    double aux2     =   0;
    double aux3     =   0;
    massDensity_    *=  this->rho_;

    for(i = 0 ; i < numCompontents_ ; i++)
    {
        aux1                +=  atomInMolecule_.at(i)*nucCharge_.at(i);
        aux2                +=  atomInMolecule_.at(i)*atomicNum_.at(i);

        SetLogConstant(i);
        SetBPrime(i);

        M_.at(i)            =   (nucCharge_.at(i)*MP + (atomicNum_.at(i) - nucCharge_.at(i))*MN)
                                /atomicNum_.at(i);
        aux3                +=  atomInMolecule_.at(i)*atomicNum_.at(i)*M_.at(i);

        if(nucCharge_.at(i)!=1)
        {
            flag=true;
        }
    }

    sumCharge_      =   aux1;
    sumNucleons_    =   aux2;
    ZA_             =   aux1/aux2;
    molDensity_     =   massDensity_*NA/aux2;
    MM_             =   aux3/aux2;
    C1_             =   2*LOG10;
    r_              =   1.31;                   // only for ice - change if needed (sea water: 1.35)

    if(flag)
    {
        mN_.resize(numCompontents_);
        Integral *integral = new Integral(IROMB, IMAXS, IPREC);
        for(i=0; i<numCompontents_; i++)
        {
            if(nucCharge_.at(i)!=1)
            {
                r0_         =   pow(atomicNum_.at(i), 1./3);
                r0_         =   1.12*r0_ - 0.86/r0_;
                mN_.at(i)   =   1 - 4*PI*0.17*integral->Integrate(r0_ , -1.0, boost::bind(&Medium::FunctionToIntegral, this, _1),3 , 2.0)
                                /atomicNum_.at(i);
            }
        }

        delete integral;

    }


    //Calculation of the radiation length
    aux1 = 0;
    aux2 = 0;
    for(i = 0 ; i < numCompontents_ ; i++)
    {
        aux1 += X0_inv(nucCharge_.at(i),atomicNum_.at(i))*(atomInMolecule_.at(i)*atomicNum_.at(i));
        aux2 += atomInMolecule_.at(i)*atomicNum_.at(i);
    }
    radiationLength_ = aux2/aux1;
    radiationLength_ /= massDensity_;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Medium::SetLogConstant(int i)
{
    int z   =   RoundValue(nucCharge_.at(i));
    switch(z)
    {
        case 1: logConstant_.at(i)  =  202.4; break;
        case 2: logConstant_.at(i)  =  151.9; break;
        case 3: logConstant_.at(i)  =  159.9; break;
        case 4: logConstant_.at(i)  =  172.3; break;
        case 5: logConstant_.at(i)  =  177.9; break;
        case 6: logConstant_.at(i)  =  178.3; break;
        case 7: logConstant_.at(i)  =  176.6; break;
        case 8: logConstant_.at(i)  =  173.4; break;
        case 9: logConstant_.at(i)  =  170.0; break;
        case 10: logConstant_.at(i) =  165.8; break;
        case 11: logConstant_.at(i) =  165.8; break;
        case 12: logConstant_.at(i) =  167.1; break;
        case 13: logConstant_.at(i) =  169.1; break;
        case 14: logConstant_.at(i) =  170.8; break;
        case 15: logConstant_.at(i) =  172.2; break;
        case 16: logConstant_.at(i) =  173.4; break;
        case 17: logConstant_.at(i) =  174.3; break;
        case 18: logConstant_.at(i) =  174.8; break;
        case 19: logConstant_.at(i) =  175.1; break;
        case 20: logConstant_.at(i) =  175.6; break;
        case 21: logConstant_.at(i) =  176.2; break;
        case 22: logConstant_.at(i) =  176.8; break;
        case 26: logConstant_.at(i) =  175.8; break;
        case 29: logConstant_.at(i) =  173.1; break;
        case 32: logConstant_.at(i) =  173.0; break;
        case 35: logConstant_.at(i) =  173.5; break;
        case 42: logConstant_.at(i) =  175.9; break;
        case 50: logConstant_.at(i) =  177.4; break;
        case 53: logConstant_.at(i) =  178.6; break;
        case 74: logConstant_.at(i) =  177.6; break;
        case 82: logConstant_.at(i) =  178.0; break;
        case 92: logConstant_.at(i) =  179.8; break;
        default: logConstant_.at(i) =  182.7;
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Medium::SetBPrime(int i)
{
    int z   =   RoundValue(nucCharge_.at(i));
    switch(z)
    {
        case 1: bPrime_.at(i)  =  446; break;
        default: bPrime_.at(i) =  1429;
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------init Media---------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


/*
* initialize water
*/

void Medium::InitWater()
{
    Inita(2);

    elementName_.at(0)                =   "H";
    elementName_.at(1)                =   "O";
    nucCharge_.at(0)        =   1; // H
    nucCharge_.at(1)        =   8; // O
    atomicNum_.at(0)        =   1.00794;
    atomicNum_.at(1)        =   15.9994;
    atomInMolecule_.at(0)   =   2;
    atomInMolecule_.at(1)   =   1;
    I_                      =   75.0;
    C_                      =   -3.5017;
    a_                      =   0.09116;
    m_                      =   3.4773;
    X0_                     =   0.2400;
    X1_                     =   2.8004;
    d0_                     =   0;
    massDensity_            =   1.000;

    Initr();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


/*
* initialize ice
*/

void Medium::InitIce()
{
    Inita(2);

    elementName_.at(0)      =   "H";
    elementName_.at(1)      =   "O";
    nucCharge_.at(0)        =   1; // H
    nucCharge_.at(1)        =   8; // O
    atomicNum_.at(0)        =   1.00794;
    atomicNum_.at(1)        =   15.9994;
    atomInMolecule_.at(0)   =   2;
    atomInMolecule_.at(1)   =   1;
    I_                      =   75.0;
    C_                      =   -3.5017;
    a_                      =   0.09116;
    m_                      =   3.4773;
    X0_                     =   0.2400;
    X1_                     =   2.8004;
    d0_                     =   0;
    massDensity_            =   0.917;

    Initr();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


/*
* initialize salt (added by Ped)
*/

void Medium::InitSalt()
{
    Inita(2);

    elementName_.at(0)      =   "Na";
    elementName_.at(1)      =   "Cl";
    nucCharge_.at(0)        =   11;
    nucCharge_.at(1)        =   17;
    atomicNum_.at(0)        =   22.98977;
    atomicNum_.at(1)        =   35.4527;
    atomInMolecule_.at(0)   =   1;
    atomInMolecule_.at(1)   =   1;

    // Calculated by ESTAR detabase
    // (it could be 185 eV by the method of reference below)

    I_                      =   175.3;

    // C through X1 are based on
    // Atomic Data and Nuclear Data Tables 78, 183-356 (2001), Appendix A

    C_                      =   -4.5041;
    a_                      =   0.1632;
    m_                      =   3;
    X0_                     =   0.2;
    X1_                     =   3.0;
    d0_                     =   0;
    massDensity_            =   2.323; // Solid halite density

    Initr();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


/*
* initialize standard rock
*/

void Medium::InitStandardrock()
{
    Inita(1);
    // Ionization potential and density corrections
    // are close to those of calcium carbonate

    elementName_.at(0)      =   "Standard_Rock";
    nucCharge_.at(0)        =   11;
    atomicNum_.at(0)        =   22;
    atomInMolecule_.at(0)   =   1;
    I_                      =   136.4;
    C_                      =   -3.7738;
    a_                      =   0.08301;
    m_                      =   3.4120;
    X0_                     =   0.0492;
    X1_                     =   3.0549;
    d0_                     =   0;
    massDensity_            =   2.650;

    Initr();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


/*
* initialize Frejus rock
*/

void Medium::InitFrejusrock()
{
    Inita(1);

    elementName_.at(0)      =   "Frejus_Rock";
    nucCharge_.at(0)        =   10.12;
    atomicNum_.at(0)        =   20.34;
    atomInMolecule_.at(0)   =   1;
    I_                      =   149.0;
    C_                      =   -5.053;
    a_                      =   0.078;
    m_                      =   3.645;
    X0_                     =   0.288;
    X1_                     =   3.196;
    d0_                     =   0;
    massDensity_            =   2.740;

    Initr();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


/*
* initialize iron
*/

void Medium::InitIron()
{
    Inita(1);

    elementName_.at(0)      =   "Fe";
    nucCharge_.at(0)        =   26;
    atomicNum_.at(0)        =   55.845;
    atomInMolecule_.at(0)   =   1;
    I_                      =   286.0;
    C_                      =   -4.2911;
    a_                      =   0.14680;
    m_                      =   2.9632;
    X0_                     =   -0.0012;
    X1_                     =   3.1531;
    d0_                     =   0.12;
    massDensity_            =   7.874;

    Initr();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


/*
* initialize hydrogen
*/

void Medium::InitHydrogen()
{
    Inita(1);

    elementName_.at(0)      =   "H";
    nucCharge_.at(0)        =   1;
    atomicNum_.at(0)        =   1.00794;
    atomInMolecule_.at(0)   =   1;
    I_                      =   21.8;
    C_                      =   -3.0977;
    a_                      =   0.13483;
    m_                      =   5.6249;
    X0_                     =   0.4400;
    X1_                     =   1.8856;
    d0_                     =   0;
    massDensity_            =   0.07080;

    Initr();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//


/*
* initialize lead
*/

void Medium::InitLead()
{
    Inita(1);

    elementName_.at(0)      =   "Pb";
    nucCharge_.at(0)        =   82;
    atomicNum_.at(0)        =   207.2;
    atomInMolecule_.at(0)   =   1;
    I_                      =   823.0;
    C_                      =   -6.2018;
    a_                      =   0.09359;
    m_                      =   3.1608;
    X0_                     =   0.3776;
    X1_                     =   3.8073;
    d0_                     =   0.14;
    massDensity_            =   11.350;

    Initr();
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

/*
* initialize copper
*/

void Medium::InitCopper()
{
    Inita(1);

    elementName_.at(0)      =   "Cu";
    nucCharge_.at(0)        =   29;
    atomicNum_.at(0)        =   63.546;
    atomInMolecule_.at(0)   =   1;
    I_                      =   322.0;
    C_                      =   -4.4190;
    a_                      =   0.14339;
    m_                      =   2.9044;
    X0_                     =   -0.0254;
    X1_                     =   3.2792;
    d0_                     =   0.08;
    massDensity_            =   8.960;

    Initr();
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


/*
* initialize uranium
*/

void Medium::InitUranium()
{
    Inita(1);

    elementName_.at(0)      =   "U";
    nucCharge_.at(0)        =   92;
    atomicNum_.at(0)        =   238.0289;
    atomInMolecule_.at(0)   =   1;
    I_                      =   890.0;
    C_                      =   -5.8694;
    a_                      =   0.19677;
    m_                      =   2.8171;
    X0_                     =   0.2260;
    X1_                     =   3.3721;
    d0_                     =   0.14;
    massDensity_            =   18.950;

    Initr();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


/*
* initialize air
*/

void Medium::InitAir()
{
    const double fr1        =   2*78.1;
    const double fr2        =   2*21.0;
    const double fr3        =   0.9;
    const double fra        =   fr1+fr2+fr3;

    Inita(3);

    elementName_.at(0)      =   "N";
    elementName_.at(1)      =   "O";
    elementName_.at(2)      =   "Ar";
    nucCharge_.at(0)        =   7; // N
    nucCharge_.at(1)        =   8; // O
    nucCharge_.at(2)        =   18; // Ar
    atomicNum_.at(0)        =   14.0067;
    atomicNum_.at(1)        =   15.9994;
    atomicNum_.at(2)        =   39.948;
    atomInMolecule_.at(0)   =   fr1/fra;
    atomInMolecule_.at(1)   =   fr2/fra;
    atomInMolecule_.at(2)   =   fr3/fra;
    I_                      =   85.7;
    C_                      =   -10.5961;
    a_                      =   0.10914;
    m_                      =   3.3994;
    X0_                     =   1.7418;
    X1_                     =   4.2759;
    d0_                     =   0;
    massDensity_            =   1.205e-3; // dry, 1 atm

    Initr();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


/*
* initialize mineral oil or paraffin CH3(CH2)~23CH3 (added by Ped)
*/

void Medium::InitParaffin()
{
    Inita(2);

    elementName_.at(0)      =   "C";
    elementName_.at(1)      =   "H";
    nucCharge_.at(0)        =   6; // C
    nucCharge_.at(1)        =   1; // H
    atomicNum_.at(0)        =   12.0011;
    atomicNum_.at(1)        =   1.0079;
    atomInMolecule_.at(0)   =   25;
    atomInMolecule_.at(1)   =   52;
    I_                      =   55.9;
    C_                      =   -2.9551;
    a_                      =   0.1209;
    m_                      =   3.4288;
    X0_                     =   0.1289;
    X1_                     =   2.5084;
    d0_                     =   0;
    massDensity_            =   0.93;

    Initr();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


/*
* initialize ANTARES water
* Sea water (Mediterranean Sea, ANTARES place)
* ==========================================================================
* WATER DENSITY CHANGES WITH THE DEPTH FROM 1.0291 g/cm^3 AT SURFACE
* UP TO 1.0404 g/cm^3 AT THE SEA BED
* (ANTARES-Site/2000-001 and references therein)
*
* The error which is caused by this simplified approach (average value for
* density) does not exceed 0.5% (much less, in fact) that is comparable with
*  an error which comes from uncertainties with the muon cross-sections.
*==========================================================================
*/

void Medium::InitAntaresWater()
{
    // added by Claudine Colnard,
    // Institute Nikhef, The Netherlands,
    // ANTARES collaboration.

    Inita(8);

    elementName_.at(0)      =   "H";
    elementName_.at(1)      =   "O";
    elementName_.at(2)      =   "Na";
    elementName_.at(3)      =   "K";
    elementName_.at(4)      =   "Mg";
    elementName_.at(5)      =   "Ca";
    elementName_.at(6)      =   "Cl";
    elementName_.at(7)      =   "S";
    nucCharge_.at(0)        =   1;      // H
    nucCharge_.at(1)        =   8;      // O
    nucCharge_.at(2)        =   11;     // Na
    nucCharge_.at(3)        =   19;     // K
    nucCharge_.at(4)        =   12;     // Mg
    nucCharge_.at(5)        =   20;     // Ca
    nucCharge_.at(6)        =   17;     // Cl
    nucCharge_.at(7)        =   16;     // S
    atomicNum_.at(0)        =   1.008;  //  Chemical composition of the seawater
    atomicNum_.at(1)        =   15.999; //  according to
    atomicNum_.at(2)        =   22.99;  //  A.Okada, Astropart. Phys. 2 (1994) 393
    atomicNum_.at(3)        =   39.10;  //  and references therein
    atomicNum_.at(4)        =   24.31;  //  corrected for Mediterranean Sea, ANTARES place
    atomicNum_.at(5)        =   40.08;  //  according to salinity  38.44+-0.02 g/kg,
    atomicNum_.at(6)        =   35.45;  //  as cited in J.Brunner, ANTARES-Site/2000-001
    atomicNum_.at(7)        =   32.07;  //  instead of 35.0 g/kg as cited in A.Okada, ...
    atomInMolecule_.at(0)   =   2;      //  (so, n[2-7] have been just multiplied by 1.098)
    atomInMolecule_.at(1)   =   1.00884;
    atomInMolecule_.at(2)   =   0.00943;
    atomInMolecule_.at(3)   =   0.000209;
    atomInMolecule_.at(4)   =   0.001087;
    atomInMolecule_.at(5)   =   0.000209;
    atomInMolecule_.at(6)   =   0.01106;
    atomInMolecule_.at(7)   =   0.00582;

     // All the same as for pure water

    I_                      =   75.0;
    C_                      =   -3.5017;
    a_                      =   0.09116;
    m_                      =   3.4773;
    X0_                     =   0.2400;
    X1_                     =   2.8004;
    d0_                     =   0;
    massDensity_            =   1.03975;

    // J.Brunner, ANTARES-Site/2000-001, the mean value
    // for sea water density at the ANTARES place between
    // sea bed D = 2400 m (1.0404 g/cm^3) and middle of
    // detector D = 2126 m (1.0391 g/cm^3).
    // Ro = 1.0341
    // for sea water density at the ANTARES place between
    // surface D = 0 m (1.0291 g/cm^3) and middle of
    // detector D = 2126 m (1.0391 g/cm^3)

    Initr();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------Functions to integrate----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Medium::FunctionToIntegral(double r)
{
    const double a  =   0.54;

    return r*r/(1+exp((r-r0_)/a));
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Medium::SetNumComponents(int numCompontents){
    numCompontents_ = numCompontents;
}

void Medium::SetNucCharge(std::vector<double> nucCharge){
    nucCharge_ = nucCharge;
}

void Medium::SetAtomicNum(std::vector<double> atomicNum){
    atomicNum_ = atomicNum;
}

void Medium::SetAtomInMolecule(std::vector<double> atomInMolecule){
    atomInMolecule_ = atomInMolecule;
}

void Medium::SetSumCharge(double sumCharge){
    sumCharge_ = sumCharge;
}

void Medium::SetZA(double ZA){
    ZA_ = ZA;
}

void Medium::SetI(double I){
    I_ = I;
}

void Medium::SetC1(double C1){
    C1_ = C1;
}

void Medium::SetC(double C){
    C_ = C;
}

void Medium::SetA(double a){
    a_ = a;
}

void Medium::SetM(double m){
    m_ = m;
}

void Medium::SetX0(double X0){
    X0_ = X0;
}

void Medium::SetX1(double X1){
    X1_ = X1;
}

void Medium::SetD0(double d0){
    d0_ = d0;
}

void Medium::SetR(double r){
    r_ = r;
}

void Medium::SetRho(double rho){
    rho_ = rho;
}

void Medium::SetMassDensity(double massDensity){
    massDensity_ = massDensity;
}

void Medium::SetMolDensity(double molDensity){
    molDensity_ = molDensity;
}

void Medium::SetAverageNucleonWeight(std::vector<double> M){
    M_ = M;
}

void Medium::SetElementName(std::vector<std::string> E){
    elementName_ = E;
}

void Medium::SetName(std::string name){
    name_ = name;
}

void Medium::SetMN(std::vector<double> mN){
    mN_ = mN;
}

void Medium::SetMM(double MM){
    MM_ = MM;
}

void Medium::SetSumNucleons(double sumNucleons){
    sumNucleons_ = sumNucleons;
}

void Medium::SetR0(double r0){
    r0_ = r0;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Destructor---------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Medium::~Medium()
{
    nucCharge_.clear();
    atomicNum_.clear();
    atomInMolecule_.clear();
    logConstant_.clear();
    bPrime_.clear();
    M_.clear();
    elementName_.clear();
    mN_.clear();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------Temporary Functions-------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double fZ(unsigned int Z)
{
    double a_sq = ALPHA*ALPHA*Z*Z;

    return a_sq*( 1./(1.+a_sq)+0.20206-0.0369*a_sq+0.0083*a_sq*a_sq-0.002*a_sq*a_sq*a_sq );
}

double Lrad(unsigned int Z)
{
    if(Z > 4) return log(184.15*pow(Z, -1./3.)); //Elements Z>4

    if(Z == 1) return 5.31;     //Hydrogen
    if(Z == 2) return 4.79;     //Helium
    if(Z == 3) return 4.74;     //Lithium
    if(Z == 4) return 4.71;     //Beryllium
    return 0;
}

double Lrad_dash(unsigned int Z)
{
    if(Z > 4) return log(1194.*pow(Z, -2./3.)); //Elements Z>4

    if(Z == 1) return 6.144;     //Hydrogen
    if(Z == 2) return 5.621;     //Helium
    if(Z == 3) return 5.805;     //Lithium
    if(Z == 4) return 5.924;     //Beryllium
    return 0;
}

double X0_inv(unsigned int Z, double M)
{
    return 4.*ALPHA*RE*RE*NA/M*( Z*Z*(Lrad(Z)-fZ(Z))+Z*Lrad_dash(Z) );
}
