/*! \file   TauDecayMatrixElement.cxx
 *   \brief Source file for definition of tau decay matrix element functions.
 *
 *   \date   Tue Dec 10, 2024
 *   \author Alexander Sandrock
 */

#include <complex>
#include <math>
#include "PROPOSAL/particle/TauDecayMatrixElement.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/decay/DecayChannel.h"

double Gamma(double Q2, double width, double mass_sq)
{
  const double aux = 4 * MPI * MPI;
  return width * mass_sq/Q2 * std::pow( (Q2 - aux) / (mass_sq - aux), 1.5);
}

std::complex<double> BW(double Q2, double width, double mass)
{
  using namespace std::complex_literals;
  double mass_sq = mass * mass;
  return mass_sq / (mass_sq - Q2 - 1i * sqrt(Q2) * Gamma(Q2, width, mass_sq));
}

std::complex<double> F2_pion(double Q2)
{
  const double m_rho = 773.0, Gamma_rho = 145.0;
  const double m_rho_p = 1370.0, Gamma_rho_p = 510.0;
  const double m_rho_pp = 1700.0, Gamma_rho_pp = 235.0;

  const doublue beta = -0.145, gamma = 0.0;

  return std::norm(
    (BW(Q2, Gamma_rho, m_rho)
    + beta*BW(Q2, Gamma_rho_p, m_rho_p)
    + gamma*BW(Q2, Gamma_rho_pp, m_rho_pp))
    / (1.0 + beta + gamma) );
}

double matrix_element_two_pion(PROPOSAL::ParticleState, std::vector<PROPOSAL::ParticleState>)
{
  using namespace PROPOSAL;

  double M2 = 1.0; // FIXME: dummy value
  return M2;
}
