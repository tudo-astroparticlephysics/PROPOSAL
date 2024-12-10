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

double matrix_element_two_pion(PROPOSAL::ParticleState parent, std::vector<PROPOSAL::ParticleState> daughters)
{
  using namespace PROPOSAL;

  // cf. Yung-Su Tsai, Phys. Rev. D 4, 2821 (1971)
  // p1 = tau, p2 = tau neutrino, q1 = pion, q2 = pion, Q = q1 - q2
  // NB: Q here is *not* the momentum for the formfactor, which corresponds to (q1 + q2) in this notation
  double p1_Q, p2_Q, p1_p2, Q2;
  p1_Q = parent.energy * (daughters[0].energy - daughters[1].energy)
    - parent.direction * (daughters[0].direction - daughters[1].direction);
  p2_Q = daughters[2].energy * (daughters[0].energy - daughters[1].energy)
    - daughters[2].direction * (daughters[0].direction - daughters[1].direction);
  p1_p2 = parent.energy * daughters[2].energy - parent.direction * daughters[2].direction;
  Q2 = 2 * MPI * MPI - 2 * (daughters[0].energy * daughters[1].energy
    - daughters[0].direction * daughters[1].direction);

  return (2 * p1_Q * p2_Q - p1_p2 * Q2) * F2_pion(2 * MPI * MPI + 2 * p1_p2);
}
