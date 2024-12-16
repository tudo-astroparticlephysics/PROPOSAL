/*! \file   TauDecayMatrixElement.cxx
 *   \brief Source file for definition of tau decay matrix element functions.
 *
 *   \date   Tue Dec 10, 2024
 *   \author Alexander Sandrock
 */

#include <complex>
#include <cmath>
#include "PROPOSAL/particle/TauDecayMatrixElement.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/decay/DecayChannel.h"

using namespace PROPOSAL;

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

double F2_pion(double Q2)
{
  const double m_rho = 773.0, Gamma_rho = 145.0;
  const double m_rho_p = 1370.0, Gamma_rho_p = 510.0;
  const double m_rho_pp = 1700.0, Gamma_rho_pp = 235.0;

  const double beta = -0.145, gamma = 0.0;

  return std::norm(
    (BW(Q2, Gamma_rho, m_rho)
    + beta*BW(Q2, Gamma_rho_p, m_rho_p)
    + gamma*BW(Q2, Gamma_rho_pp, m_rho_pp))
    / (1.0 + beta + gamma) );
}

std::complex<double> B_rho(double s)
{
  const double m_rho = 773.0, Gamma_rho = 145.0;
  const double m_rho_p = 1370.0, Gamma_rho_p = 510.0;
  const double beta = -0.145;

  return (BW(s, Gamma_rho, m_rho) + beta * BW(s, Gamma_rho_p, m_rho_p))/(1.0 + beta);
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

double matrix_element_three_pion(PROPOSAL::ParticleState parent, std::vector<PROPOSAL::ParticleState> daughters)
{
  using namespace PROPOSAL;

  // cf. Kühn & Santamaria, Z. Phys. C 48, 445 (1990)
  // p1 = Tau == parent, p2 = NuTau == daughters[3],
  // q1 = PiMinus, q2 = Pi0, q3 = Pi0 or q1 = PiMinus, q2 = PiMinus, q3 = PiPlus
  double Q2, s1, s2, s3;

  s1 = 2* MPI * MPI + 2 * (daughters[1].energy * daughters[2].energy
     - daughters[1].direction * daughters[2].direction);
  s2 = 2* MPI * MPI + 2 * (daughters[0].energy * daughters[2].energy
     - daughters[0].direction * daughters[2].direction);
  s3 = 2* MPI * MPI + 2 * (daughters[0].energy * daughters[1].energy
     - daughters[0].direction * daughters[1].direction);
  Q2 = s1 + s2 + s3 - 3* MPI * MPI;

  double B2_s1, B2_s2, Bs1_Bs2;
  B2_s1 = std::norm(B_rho(s1));
  B2_s2 = std::norm(B_rho(s2));
  Bs1_Bs2 = std::real(B_rho(s1) * std::conj(B_rho(s2)));

  double p1_q1, p1_q2, p1_q3;
  p1_q1 = parent.energy*daughters[0].energy - parent.direction * daughters[0].direction;
  p1_q2 = parent.energy*daughters[1].energy - parent.direction * daughters[1].direction;
  p1_q3 = parent.energy*daughters[2].energy - parent.direction * daughters[2].direction;

  return B2_s2 * ((MTAU*MTAU - Q2)/Q2 * (s2*s2 + 4*s1*s2 - 6*MPI*MPI*s2 + 2*Q2*s2 + 4*s1*s1
      - 12*MPI*MPI*s1 - 4*Q2*s1 + 9*std::pow(MPI, 4) - 10*Q2*MPI*MPI + Q2*Q2)
    + 8*p1_q1*p1_q2/(Q2*Q2)*(s2 + 2*s1 - 3*MPI*MPI - Q2)*(s2 + 2*s1 - 3*MPI*MPI + Q2)
    + 8*p1_q1*p1_q3/(Q2*Q2)*(s2 + 2*s1 - 3*MPI*MPI - 3*Q2)*(s2 + 2*s1 - 3*MPI*MPI + Q2)
    + 8*p1_q2*p1_q3/(Q2*Q2)*(s2 + 2*s1 - 3*MPI*MPI - 3*Q2)*(s2 + 2*s1 - 3*MPI*MPI - Q2)
    + 4*p1_q1*p1_q1/(Q2*Q2)*std::pow(s2 + 2*s1 - 3*MPI*MPI + Q2, 2)
    + 4*p1_q2*p1_q2/(Q2*Q2)*std::pow(s2 + 2*s1 - 3*MPI*MPI - Q2, 2)
    + 4*p1_q3*p1_q3/(Q2*Q2)*std::pow(s2 + 2*s1 - 3*MPI*MPI - 3*Q2, 2) )
  + B2_s1 * ((MTAU*MTAU - Q2)/Q2 * (4*s2*s2 + 4*s1*s2 - 12*MPI*MPI*s2 - 4*Q2*s2 + s1*s1
      - 6*MPI*MPI*s1 + 2*Q2*s1 + 9*std::pow(MPI, 4) - 10*Q2*MPI*MPI + Q2*Q2)
    + 8*p1_q1*p1_q2*(2*s2 + s1 - 3*MPI*MPI - Q2)*(2*s2 + s1 - 3*MPI*MPI + Q2)
    + 8*p1_q1*p1_q3*(2*s2 + s1 - 3*MPI*MPI - 3*Q2)*(2*s2 + s1 - 3*MPI*MPI - Q2)
    + 8*p1_q2*p1_q3*(2*s2 + s1 - 3*MPI*MPI - 3*Q2)*(2*s2 + s1 - 3*MPI*MPI + Q2)
    + 4*p1_q1*p1_q1/(Q2*Q2)*std::pow(2*s2 + s1 - 3*MPI*MPI - Q2 ,2)
    + 4*p1_q2*p1_q2/(Q2*Q2)*std::pow(2*s2 + s1 - 3*MPI*MPI + Q2 ,2)
    + 4*p1_q3*p1_q3/(Q2*Q2)*std::pow(2*s2 + s1 - 3*MPI*MPI - 3*Q2 ,2) )
  + Bs1_Bs2 * (2*(MTAU*MTAU - Q2)/Q2*(2*s2*s2 + 5*s1*s2 - 9*MPI*MPI*s2 + Q2*s2 + 2*s1*s1
      - 9*MPI*MPI*s1 + Q2*s1 + 9*std::pow(MPI, 4) - 8*Q2*MPI*MPI - Q2*Q2)
    + 16*p1_q1*p1_q2*(2*s2*s2 + 5*s1*s2 - 9*MPI*MPI*s2 + 2*s1*s1 - 9*MPI*MPI*s1
      + 9*std::pow(MPI, 4) + Q2*Q2)
    + 16*p1_q1*p1_q3*(2*s2*s2 + 5*s1*s2 - 9*MPI*MPI*s2 - 4*Q2*s2 + 2*s1*s1 - 9*MPI*MPI*s1
      - 5*Q2*s1 + 9*std::pow(MPI, 4) + 9*Q2*MPI*MPI)
    + 16*p1_q2*p1_q3*(2*s2*s2 + 5*s1*s2 - 9*MPI*MPI*s2 - 5*Q2*s2 + 2*s1*s1 - 9*MPI*MPI*s1
      - 4*Q2*s1 + 9*std::pow(MPI, 4) + 9*Q2*MPI*MPI)
    + 8*p1_q1*p1_q1/(Q2*Q2)*(s2 + 2*s1 - 3*MPI*MPI + Q2)*(2*s2 + s1 - 3*MPI*MPI - Q2)
    + 8*p1_q2*p1_q2/(Q2*Q2)*(s2 + 2*s1 - 3*MPI*MPI - Q2)*(2*s2 + s1 - 3*MPI*MPI + Q2)
    + 8*p1_q3*p1_q3/(Q2*Q2)*(s2 + 2*s1 - 3*MPI*MPI - 3*Q2)*(2*s2 + s1 - 3*MPI*MPI - 3*Q2) )
  ;
}
