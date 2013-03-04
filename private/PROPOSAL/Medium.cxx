/*! \file   Medium.cxx
*   \brief  Source file for the medium routines.
*
*   For more details see the class documentation.
*
*   \date   21.06.2010
*   \author Jan-Hendrik Koehne
*/


#include "PROPOSAL/Medium.h"
#include <cmath>
#include "PROPOSAL/methods.h"

using namespace std;

//----------------------------------------------------------------------------//

Medium::Medium(){}

//----------------------------------------------------------------------------//

Medium::~Medium(){}

//----------------------------------------------------------------------------//


Medium::Medium(string w, double ecut, double vcut, double rho)
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

    if(equalsIgnoreCase(w,"water"))
    {
        initWater();
    }
    else if(equalsIgnoreCase(w,"ice"))
    {
        initIce();
    }
    else if(equalsIgnoreCase(w,"salt"))
    {
        initSalt();
    }
    else if(equalsIgnoreCase(w,"standard rock"))
    {
        initStandardrock();
    }
    else if(equalsIgnoreCase(w,"frejus rock"))
    {
        initFrejusrock();
    }
    else if(equalsIgnoreCase(w,"iron"))
    {
        initIron();
    }
    else if(equalsIgnoreCase(w,"hydrogen"))
    {
        initHydrogen();
    }
    else if(equalsIgnoreCase(w,"lead"))
    {
        initLead();
    }
    else if(equalsIgnoreCase(w,"uranium"))
    {
        initUranium();
    }
    else if(equalsIgnoreCase(w,"air"))
    {
        initAir();
    }
    else if(equalsIgnoreCase(w,"mineral oil"))
    {
        initParaffin();
    }
    else if(equalsIgnoreCase(w,"antares water"))
    {
        initAntaresWater();
    }
    else
    {
        printf("Warning (in Medium/Medium): unknown medium: defaulting to water");
        name_   =   "water";
        initWater();
    }

    this->ecut_ =   ecut;
    this->vcut_ =   vcut;
    vCut_       =   1.;
}


//----------------------------------------------------------------------------//

void Medium::inita(int i)
{
    numCompontents_ =   i;

    nucCharge_.resize(numCompontents_);
    atomicNum_.resize(numCompontents_);
    atomInMolecule_.resize(numCompontents_);
    logConstant_.resize(numCompontents_);
    bPrime_.resize(numCompontents_);
    M_.resize(numCompontents_);
    E_.resize(numCompontents_);
}

//----------------------------------------------------------------------------//


void Medium::initr()
{
    int i;
    bool flag       =   false;
    double aux1     =   0;
    double aux2     =   0;
    double aux3     =   0;
    massDensity_    *=  this->rho_;

    for(i=0; i<numCompontents_; i++)
    {
        aux1                +=  atomInMolecule_.at(i)*nucCharge_.at(i);
        aux2                +=  atomInMolecule_.at(i)*atomicNum_.at(i);
        logConstant_.at(i)  =   elB(i);
        bPrime_.at(i)       =   elP(i);
        M_.at(i)            =   (nucCharge_.at(i)*MP + (atomicNum_.at(i) - nucCharge_.at(i))*MN)/atomicNum_.at(i);
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
    ecut_           =   ME/sqrt(1-1/(r_*r_));   // in order to emit Cerenkov photons

    if(flag)
    {
        mN_.resize(numCompontents_);
        Integral *integral_ = new Integral(IROMB, IMAXS, IPREC);

        for(i=0; i<numCompontents_; i++)
        {
            if(nucCharge_.at(i)!=1)
            {
                r0_         =   pow(atomicNum_.at(i), 1./3);
                r0_         =   1.12*r0_ - 0.86/r0_;
                mN_.at(i)   =   1 - 4*PI*0.17*integral_->integrateWithSubstitution(r0_ , -1.0, this , 2.0)/atomicNum_.at(i);
            }
        }
    }
}

//----------------------------------------------------------------------------//

double Medium::elB(int i)
{
    int z   =   roundValue(nucCharge_.at(i));

    switch(z)
    {
        case 1: return 202.4;
        case 2: return 151.9;
        case 3: return 159.9;
        case 4: return 172.3;
        case 5: return 177.9;
        case 6: return 178.3;
        case 7: return 176.6;
        case 8: return 173.4;
        case 9: return 170.0;
        case 10: return 165.8;
        case 11: return 165.8;
        case 12: return 167.1;
        case 13: return 169.1;
        case 14: return 170.8;
        case 15: return 172.2;
        case 16: return 173.4;
        case 17: return 174.3;
        case 18: return 174.8;
        case 19: return 175.1;
        case 20: return 175.6;
        case 21: return 176.2;
        case 22: return 176.8;
        case 26: return 175.8;
        case 29: return 173.1;
        case 32: return 173.0;
        case 35: return 173.5;
        case 42: return 175.9;
        case 50: return 177.4;
        case 53: return 178.6;
        case 74: return 177.6;
        case 82: return 178.0;
        case 92: return 179.8;
        default: return 182.7;
    }
}

//----------------------------------------------------------------------------//

double Medium::elP(int i)
{
    int z   =   roundValue(nucCharge_.at(i));

    switch(z)
    {
        case 1: return 446;
        default: return 1429;
    }
}

//----------------------------------------------------------------------------//

/*
* initialize water
*/

void Medium::initWater()
{
    inita(2);

    E_.at(0)                =   "H";
    E_.at(1)                =   "O";
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

    initr();
}

//----------------------------------------------------------------------------//

/*
* initialize ice
*/

void Medium::initIce()
{
    inita(2);

    E_.at(0)                =   "H";
    E_.at(1)                =   "O";
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

    initr();
}

//----------------------------------------------------------------------------//

/*
* initialize salt (added by Ped)
*/

void Medium::initSalt()
{
    inita(2);

    E_.at(0)                =   "Na";
    E_.at(1)                =   "Cl";
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

    initr();
}

//----------------------------------------------------------------------------//

/*
* initialize standard rock
*/

void Medium::initStandardrock()
{
    inita(1);
    // Ionization potential and density corrections
    // are close to those of calcium carbonate

    E_.at(0)                =   "Standard Rock";
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

    initr();
}

//----------------------------------------------------------------------------//

/*
* initialize Frejus rock
*/

void Medium::initFrejusrock()
{
    inita(1);

    E_.at(0)                =   "Frejus Rock";
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

    initr();
}


//----------------------------------------------------------------------------//

/*
* initialize iron
*/

void Medium::initIron()
{
    inita(1);

    E_.at(0)                =   "Fe";
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

    initr();
}

//----------------------------------------------------------------------------//

/*
* initialize hydrogen
*/

void Medium::initHydrogen()
{
    inita(1);

    E_.at(0)                =   "H";
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

    initr();
}

//----------------------------------------------------------------------------------------------------//

/*
* initialize lead
*/

void Medium::initLead()
{
    inita(1);

    E_.at(0)                =   "Pb";
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

    initr();
}

//----------------------------------------------------------------------------//

/*
* initialize uranium
*/

void Medium::initUranium()
{
    inita(1);

    E_.at(0)                =   "U";
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

    initr();
}

//----------------------------------------------------------------------------//

/*
* initialize air
*/

void Medium::initAir()
{
    const double fr1        =   2*78.1;
    const double fr2        =   2*21.0;
    const double fr3        =   0.9;
    const double fra        =   fr1+fr2+fr3;

    inita(3);

    E_.at(0)                =   "N";
    E_.at(1)                =   "O";
    E_.at(2)                =   "Ar";
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

    initr();
}

//----------------------------------------------------------------------------//

/*
* initialize mineral oil or paraffin CH3(CH2)~23CH3 (added by Ped)
*/

void Medium::initParaffin()
{
    inita(2);

    E_.at(0)                =   "C";
    E_.at(1)                =   "H";
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

    initr();
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

void Medium::initAntaresWater()
{
    // added by Claudine Colnard,
    // Institute Nikhef, The Netherlands,
    // ANTARES collaboration.

    inita(8);

    E_.at(0)                =   "H";
    E_.at(1)                =   "O";
    E_.at(2)                =   "Na";
    E_.at(3)                =   "K";
    E_.at(4)                =   "Mg";
    E_.at(5)                =   "Ca";
    E_.at(6)                =   "Cl";
    E_.at(7)                =   "S";
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

    initr();
}

//----------------------------------------------------------------------------//


double Medium::vCut(double E)
{
    double aux;

    aux =   ecut_/E;

    if(ecut_>0)
    {
        if(vcut_>0 && vcut_<=1)
        {
            if(aux<vcut_)
            {
                vCut_   =   aux;
            }
            else
            {
                vCut_   =   vcut_;
            }
        }
        else
        {
            vCut_   =   aux;
        }
    }
    else
    {
        if(vcut_>0 && vcut_<=1)
        {
            vCut_   =   vcut_;
        }
        else
        {
            vCut_   =   1.;
        }
    }

    return vCut_;
}

//----------------------------------------------------------------------------//

double Medium::function(double r)
{
    const double a  =   0.54;

    return r*r/(1+exp((r-r0_)/a));
}



