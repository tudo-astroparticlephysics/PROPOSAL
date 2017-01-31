#ifndef CONSTANTS_H
#define CONSTANTS_H

#ifndef PI
    #define PI 3.141592653589793
#endif

#ifndef IROMB // romb # for integration
    #define IROMB 5
#endif

#ifndef IMAXS // max number of int. steps
    #define IMAXS 40
#endif

#ifndef IPREC // integration precision
    #define IPREC 1.e-6//100
#endif

#ifndef IPREC2 // integration precision
    #define IPREC2 1.e-6*10//100
#endif

#ifndef NUM1 // number of interpolation
    #define NUM1 100
#endif

#ifndef NUM2 // points specified mainly
    #define NUM2 200
#endif

#ifndef NUM3 // in Propagate.cxx
    #define NUM3 1000
#endif

#ifndef COMPUTER_PRECISION
    #define COMPUTER_PRECISION 1.e-10
#endif

#ifndef HALF_PRECISION
    #define HALF_PRECISION 1.e-5 //std::sqrt(computerPrecision);
#endif

#ifndef GEOMETRY_PRECISION
    #define GEOMETRY_PRECISION 1.e-9 //std::sqrt(computerPrecision);
#endif

#ifndef ALPHA // fine structure constant
    #define ALPHA 0.007297352533285885
#endif

#ifndef ME // electron mass (MeV)
    #define ME 0.510998902
#endif

#ifndef RY // Rydberg energy (eV)
    #define RY 13.60569172
#endif

#ifndef IONK // in the ionization formula (MeV*cm2/g)
    #define IONK 0.307075
#endif

#ifndef SPEED // speed of light (cm/s)
    #define SPEED 2.99792458e10
#endif

#ifndef RE // classical electron radius (cm)
    #define RE 2.817940285e-13
#endif

#ifndef E0
    #define E0 1.602176462e-19                         //charge of an electron in A*s (SI)
#endif

#ifndef RM // classical Muon radius(cm) :   RM = (1.602176487*pow(10,-19))/(4*PI*8.854187817*pow(10,-12)*Mmu_)*pow(10,-4);
    #define RM 1.362849110866631e-15
#endif

#ifndef NA // Avogadro's number (1/mol)
    #define NA 6.02214199e23
#endif

#ifndef MMU // muon mass (MeV)
    #define MMU 105.658389
#endif

#ifndef LMU // muon lifetime (sec)
    #define LMU 2.19703e-6
#endif

#ifndef MTAU // tau mass (MeV)
    #define MTAU 1777.03
#endif

#ifndef LTAU // tau lifetime (sec)
    #define LTAU 290.6e-15
#endif

#ifndef MPI // pion mass (MeV)
    #define MPI 139.57018
#endif

#ifndef MP // proton mass (MeV)
    #define MP 938.271998
#endif

#ifndef MN // neutron mass (MeV)
    #define MN 939.56533
#endif

#ifndef MRH // rho-770 mass (MeV)
    #define MRH 769.3
#endif

#ifndef MA1 // a1-1260 mass (MeV)
    #define MA1 1230
#endif

#ifndef MRS // rho-1450 mass (MeV)
    #define MRS 1465
#endif

#ifndef GF // Fermi coupling const staticant (MeV^-2)
    #define GF 1.16639e-11
#endif

#ifndef MW // W+- boson mass (MeV)
    #define MW 80419.
#endif

#ifndef MZ // Z0 boson mass (MeV)
    #define MZ 91188.2
#endif

#ifndef XW // sin^2(mixing angle at Mz)
    #define XW 0.23117
#endif

#ifndef GW // W+- boson width (MeV)
    #define GW 2120.
#endif

#ifndef GZ // Z0 boson width (MeV)
    #define GZ 2495.2
#endif

#ifndef DS2 // Sun neutrino mass difference (eV^2)
    #define DS2 6.9e-5
#endif

#ifndef TT2 // tan^2(th) of the mixing angle
    #define TT2 0.43
#endif

#ifndef DE2 // Earth neutrino mass difference (eV^2)
    #define DE2 2.6e-3
#endif

#ifndef ST2 // sin^2(2 th) of the mixing angle
    #define ST2 1.0
#endif

#ifndef MMON // monopole mass (MeV)
    #define MMON 1.e5
#endif

#ifndef CMON // monopole charge (in units of e)
    #define CMON 68.51799988 //CMON=1/(2*ALPHA);
#endif

#ifndef MSTAU // stau mass (MeV)
    #define MSTAU 1.e5
#endif

#ifndef LSTAU // stau lifetime (sec)
    #define LSTAU -1
#endif

#ifndef XRES // resolution of particle position (cm)
    #define XRES 1.e-3
#endif

#ifndef BIGENERGY // used for radiation length (MeV)
    #define BIGENERGY 1.e14
#endif

#ifndef ELOW // bounds of parameterizations
    #define ELOW 0
#endif

#ifndef NLOW // bounds of parameterizations
    #define NLOW  0.510998902
#endif

#ifndef LOG10 // log(10)
    #define LOG10 2.302585092994046
#endif

#ifndef SQRT2 // sqrt(2)
    #define SQRT2 1.414213562373095
#endif

#ifndef SQRT3 // sqrt(3)
    #define SQRT3 1.732050807568877
#endif

#ifndef SQRTE // sqrt(e)
    #define SQRTE 1.648721270700128
#endif


#endif // CONSTANTS_H
