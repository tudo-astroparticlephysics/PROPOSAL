
#include "PROPOSAL/Constants.h"

// numbers
const double PROPOSAL::PI               = 3.141592653589793;
const double PROPOSAL::LOG10            = 2.302585092994046;                      // log(10)
const double PROPOSAL::SQRT2            = 1.414213562373095;                      // sqrt(2)
const double PROPOSAL::SQRT3            = 1.732050807568877;                      // sqrt(3)
const double PROPOSAL::SQRTE            = 1.648721270700128;                      // sqrt(e)
const double PROPOSAL::EULER_MASCHERONI = 0.577215664901532860606512090082402431; // Euler-Mascheroni constant

// integration parameters
const int PROPOSAL::IROMB     = 5;          // romb # for integration
const int PROPOSAL::IMAXS     = 40;         // max number of int. steps
const double PROPOSAL::IPREC  = 1.e-6;      // integration precision
const double PROPOSAL::IPREC2 = 1.e-6 * 10; // integration precision

// precision parameters
const double PROPOSAL::COMPUTER_PRECISION           = 1.e-10;
const double PROPOSAL::HALF_PRECISION               = 1.e-5; // std::sqrt(computerPrecision);
const double PROPOSAL::GEOMETRY_PRECISION           = 1.e-9;
const double PROPOSAL::PARTICLE_POSITION_RESOLUTION = 1.e-3; // resolution of particle position (cm)

// physical constants
const double PROPOSAL::ALPHA = 0.0072973525664; // fine structure constant
const double PROPOSAL::RY    = 13.605693009;    // Rydberg energy (eV)
const double PROPOSAL::NA    = 6.022140857e23;  // Avogadro's number (1/mol)
const double PROPOSAL::SPEED = 2.99792458e10;   // speed of light (cm/s)
const double PROPOSAL::IONK  = 0.307075;        // Ionization Constant = 4*PI*NA*RE*RE*ME (MeV*cm2/g)
const double PROPOSAL::HBAR  = 6.58211928e-22;  // hbar in MeV*s

// particle constants

const double PROPOSAL::MP = 938.2720813; // proton mass (MeV)
const double PROPOSAL::MN = 939.565413;  // neutron mass (MeV)

const double PROPOSAL::ME = 0.5109989461;     // electron mass (MeV)
const double PROPOSAL::RE = 2.8179403227e-13; // classical electron radius (cm)

const double PROPOSAL::MMU = 105.6583745;  // muon mass (MeV)
const double PROPOSAL::LMU = 2.1969811e-6; // muon lifetime (sec)

const double PROPOSAL::MTAU = 1776.86;   // tau mass (MeV)
const double PROPOSAL::LTAU = 290.3e-15; // tau lifetime (sec)

const double PROPOSAL::MPI = 139.57018; // charged pion mass (MeV)
const double PROPOSAL::LPI = 2.6033e-8; // charged pion lifetime (sec)

const double PROPOSAL::MPI0 = 134.9766; // pion 0 mass (MeV)
const double PROPOSAL::LPI0 = 8.52e-17; // pion 0 lifetime (sec)

const double PROPOSAL::MKAON0 = 497.614;   // uncharged kaon mass (MeV)
const double PROPOSAL::MKAON  = 493.677;   // charged kaon mass (MeV)
const double PROPOSAL::LKAON  = 1.2380e-8; // charged kaon lifetime (sec)

const double PROPOSAL::MRH = 775.26; // rho-770 mass (MeV)
const double PROPOSAL::MA1 = 1230;   // a1-1260 mass (MeV)
const double PROPOSAL::MRS = 1465;   // rho-1450 mass (MeV)

const double PROPOSAL::MMON = 1.e5;        // monopole mass (MeV)
const double PROPOSAL::CMON = 68.51799988; // monopole charge (in units of e) = 1/(2*ALPHA)

const double PROPOSAL::MSMP  = 1.e5; // Stable massive particle mass (MeV)
const double PROPOSAL::MSTAU = 1.e5; // stau mass (MeV)

const double PROPOSAL::STABLE_PARTICLE = -1.; // lifetime of stable particle, -1 because of history
