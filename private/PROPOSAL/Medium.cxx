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


using namespace std;

//----------------------------------------------------------------------------//

Medium::Medium(){}

//----------------------------------------------------------------------------//

Medium::Medium(const Medium &medium)
{
    *this = medium;
}
//----------------------------------------------------------------------------//

Medium& Medium::operator=(const Medium &medium){
    return *this;
}

//----------------------------------------------------------------------------//

Medium::~Medium(){}

//----------------------------------------------------------------------------//


Medium::Medium(string w, double rho)
{
    name_   =   w;

    if (rho > 0 )
    {
        this->rho_  =   rho;
    }

    else
    {
        this->rho_  =   1;
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
    else if(EqualsIgnoreCase(w,"standard rock"))
    {
        InitStandardrock();
    }
    else if(EqualsIgnoreCase(w,"frejus rock"))
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
    else if(EqualsIgnoreCase(w,"uranium"))
    {
        InitUranium();
    }
    else if(EqualsIgnoreCase(w,"air"))
    {
        InitAir();
    }
    else if(EqualsIgnoreCase(w,"mineral oil"))
    {
        InitParaffin();
    }
    else if(EqualsIgnoreCase(w,"antares water"))
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
                mN_.at(i)   =   1 - 4*PI*0.17*integral->IntegrateWithSubstitution(r0_ , -1.0, FunctionToIntegral , 2.0)
                                /atomicNum_.at(i);
            }
        }

        delete integral;

    }
}

//----------------------------------------------------------------------------//

void Medium::SetLogConstant(int i)
{
    int z   =   RoundValue(nucCharge_.at(i));
    switch(z)
    {
        case 1: logConstant_.at(i)  =  202.4;
        case 2: logConstant_.at(i)  =  151.9;
        case 3: logConstant_.at(i)  =  159.9;
        case 4: logConstant_.at(i)  =  172.3;
        case 5: logConstant_.at(i)  =  177.9;
        case 6: logConstant_.at(i)  =  178.3;
        case 7: logConstant_.at(i)  =  176.6;
        case 8: logConstant_.at(i)  =  173.4;
        case 9: logConstant_.at(i)  =  170.0;
        case 10: logConstant_.at(i) =  165.8;
        case 11: logConstant_.at(i) =  165.8;
        case 12: logConstant_.at(i) =  167.1;
        case 13: logConstant_.at(i) =  169.1;
        case 14: logConstant_.at(i) =  170.8;
        case 15: logConstant_.at(i) =  172.2;
        case 16: logConstant_.at(i) =  173.4;
        case 17: logConstant_.at(i) =  174.3;
        case 18: logConstant_.at(i) =  174.8;
        case 19: logConstant_.at(i) =  175.1;
        case 20: logConstant_.at(i) =  175.6;
        case 21: logConstant_.at(i) =  176.2;
        case 22: logConstant_.at(i) =  176.8;
        case 26: logConstant_.at(i) =  175.8;
        case 29: logConstant_.at(i) =  173.1;
        case 32: logConstant_.at(i) =  173.0;
        case 35: logConstant_.at(i) =  173.5;
        case 42: logConstant_.at(i) =  175.9;
        case 50: logConstant_.at(i) =  177.4;
        case 53: logConstant_.at(i) =  178.6;
        case 74: logConstant_.at(i) =  177.6;
        case 82: logConstant_.at(i) =  178.0;
        case 92: logConstant_.at(i) =  179.8;
        default: logConstant_.at(i) =  182.7;
    }
}

//----------------------------------------------------------------------------//

void Medium::SetBPrime(int i)
{
    int z   =   RoundValue(nucCharge_.at(i));
    switch(z)
    {
        case 1: bPrime_.at(i)  =  446;
        default: bPrime_.at(i) =  1429;
    }
}

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

/*
* initialize standard rock
*/

void Medium::InitStandardrock()
{
    Inita(1);
    // Ionization potential and density corrections
    // are close to those of calcium carbonate

    elementName_.at(0)      =   "Standard Rock";
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

/*
* initialize Frejus rock
*/

void Medium::InitFrejusrock()
{
    Inita(1);

    elementName_.at(0)      =   "Frejus Rock";
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

double Medium::FunctionToIntegral(double r)
{
    const double a  =   0.54;

    return r*r/(1+exp((r-r0_)/a));
}

//----------------------------------------------------------------------------//
    // Setter

    void Medium::SetNumCompontents(int numCompontents){
        numCompontents_ = numCompontents;
    }
//----------------------------------------------------------------------------//
    void Medium::SetNucCharge(std::vector<double> nucCharge){
        nucCharge_ = nucCharge;
    }
//----------------------------------------------------------------------------//
    void Medium::SetAtomicNum(std::vector<double> atomicNum){
        atomicNum_ = atomicNum;
    }
//----------------------------------------------------------------------------//
    void Medium::SetAtomInMolecule(std::vector<double> atomInMolecule){
        atomInMolecule_ = atomInMolecule;
    }
//----------------------------------------------------------------------------//
    void Medium::SetSumCharge(double sumCharge){
        sumCharge_ = sumCharge;
    }
//----------------------------------------------------------------------------//
    void Medium::SetZA(double ZA){
        ZA_ = ZA;
    }
//----------------------------------------------------------------------------//
    void Medium::SetI(double I){
        I_ = I;
    }
//----------------------------------------------------------------------------//
    void Medium::SetC1(double C1){
        C1_ = C1;
    }
//----------------------------------------------------------------------------//
    void Medium::SetC(double C){
        C_ = C;
    }
//----------------------------------------------------------------------------//
    void Medium::SetA(double a){
        a_ = a;
    }
//----------------------------------------------------------------------------//
    void Medium::SetM(double m){
        m_ = m;
    }
//----------------------------------------------------------------------------//
    void Medium::SetX0(double X0){
        X0_ = X0;
    }
//----------------------------------------------------------------------------//
    void Medium::SetX1(double X1){
        X1_ = X1;
    }
//----------------------------------------------------------------------------//
    void Medium::SetD0(double d0){
        d0_ = d0;
    }
//----------------------------------------------------------------------------//
    void Medium::SetR(double r){
        r_ = r;
    }
//----------------------------------------------------------------------------//
    void Medium::SetRho(double rho){
        rho_ = rho;
    }
//----------------------------------------------------------------------------//
    void Medium::SetMassDensity(double massDensity){
        massDensity_ = massDensity;
    }
//----------------------------------------------------------------------------//
    void Medium::SetMolDensity(double molDensity){
        molDensity_ = molDensity;
    }
//----------------------------------------------------------------------------//
    void Medium::SetAverageNucleonWeight(std::vector<double> M){
        M_ = M;
    }
//----------------------------------------------------------------------------//
    void Medium::SetElementName(std::vector<std::string> E){
        elementName_ = E;
    }
//----------------------------------------------------------------------------//
    void Medium::SetName(std::string name){
        name_ = name;
    }
//----------------------------------------------------------------------------//
    void Medium::SetMN(std::vector<double> mN){
        mN_ = mN;
    }
//----------------------------------------------------------------------------//
    void Medium::SetMM(double MM){
        MM_ = MM;
    }
//----------------------------------------------------------------------------//
    void Medium::SetSumNucleons(double sumNucleons){
        sumNucleons_ = sumNucleons;
    }
//----------------------------------------------------------------------------//
    void Medium::SetR0(double r0){
        r0_ = r0;
    }

//----------------------------------------------------------------------------//

