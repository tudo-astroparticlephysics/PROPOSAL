#pragma once

#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace PROPOSAL
{
    // numbers
    extern const double PI;
    extern const double LOG10; // log(10)
    extern const double SQRT2; // sqrt(2)
    extern const double SQRT3; // sqrt(3)
    extern const double SQRTE; // sqrt(e)
    extern const double EULER_MASCHERONI; // Euler-Mascheroni constant

    // integration parameters
    extern const int IROMB; // romb # for integration
    extern const int IMAXS; // max number of int. steps
    extern const double IPREC; // integration precision
    extern const double IPREC2; // integration precision

    // interpolation parameters
    extern const int NUM1; // number of interpolation in cross section
    extern const int NUM2; // number of interpolation in continuous randomization
    extern const int NUM3; // number of interpolation in propagate

    extern const double BIGENERGY; // upper energy bound for Interpolation (MeV)

    // precision parameters
    extern const double COMPUTER_PRECISION;
    extern const double HALF_PRECISION; //std::sqrt(computerPrecision);
    extern const double GEOMETRY_PRECISION;
    extern const double PARTICLE_POSITION_RESOLUTION; // resolution of particle position (cm)

    // physical constants
    extern const double ALPHA; // fine structure constant
    extern const double RY; // Rydberg energy (eV)
    extern const double NA; // Avogadro's number (1/mol)
    extern const double SPEED; // speed of light (cm/s)
    extern const double IONK; // Ionization Constant = 4*PI*NA*RE*RE*ME (MeV*cm2/g)
    extern const double HBAR; // hbar in MeV*s

    // particle constants

    extern const double ME; // electron mass (MeV)
    extern const double RE; // classical electron radius (cm)

    extern const double MMU; // muon mass (MeV)
    extern const double LMU; // muon lifetime (sec)

    extern const double MTAU; // tau mass (MeV)
    extern const double LTAU; // tau lifetime (sec)

    extern const double MPI; // charged pion mass (MeV)
    extern const double LPI; // charged pion lifetime (sec)

    extern const double MPI0; // pion 0 mass (MeV)
    extern const double LPI0; // pion 0 lifetime (sec)

    extern const double MKAON; // charged kaon mass (MeV)
    extern const double LKAON; // charged kaon lifetime (sec)

    extern const double MP; // proton mass (MeV)
    extern const double MN; // neutron mass (MeV)

    extern const double MRH; // rho-770 mass (MeV)
    extern const double MA1; // a1-1260 mass (MeV)
    extern const double MRS; // rho-1450 mass (MeV)

    extern const double MMON; // monopole mass (MeV)
    extern const double CMON; // monopole charge (in units of e) = 1/(2*ALPHA)

    extern const double MSMP; // Stable massive particle mass (MeV)
    extern const double MSTAU; // stau mass (MeV)

    extern const double STABLE_PARTICLE; // lifetime of stable particle, -1 because of history

}

#endif // CONSTANTS_H
