#pragma once

#ifndef Ionization_H
#define Ionization_H

#include "PROPOSAL/crossection/CrossSections.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/Interpolant.h"


namespace PROPOSAL
{
    class Ionization;
}

std::ostream& operator<<(std::ostream& os, PROPOSAL::Ionization const &ioniz);

namespace PROPOSAL{

class Ionization: public CrossSections
{
protected:

    // double beta_;
    // double gamma_;

    Integral*   integral_;
    Interpolant* dedx_interpolant_;

    Interpolant* dndx_interpolant_1d_;
    Interpolant* dndx_interpolant_2d_;

//----------------------------------------------------------------------------//

    /*!
     * sets the parameters \n
     * \f$v_{min}=\frac{10^{-6}I}{E_p}\f$,
     * \f$v_{Max}=min\Big(v_{Max}, 1-\frac{m_p}{E_p}\Big)\f$
     * with \f$v_{Max}=\frac{2m_e(\gamma^2-1)}{1+2\gamma
     * \frac{m_e}{m_p}+\Big(\frac{m_e}{m_p}\Big)^2E}\f$ and
     * \f$v_{up}=min(v_{Max}, v_{Cut})\f$
     */

    CrossSections::IntegralLimits SetIntegralLimits(double energy, int component);

//----------------------------------------------------------------------------//
    /*!
     * \brief inelastic electron bremsstrahlung correction to
     * dEdx - Interface to Integral:
     *
     * \f[f(r)= v \frac{d^2 N}{dv dx} \cdot \Delta_{inel}\f]
     *
     * \param   v   relative energy loss
     * \return  inelastic electron bremsstrahlung correction to dEdx
     */


    double FunctionToDEdxIntegral(double energy, double variable);

//----------------------------------------------------------------------------//

    /*!
     * \brief density correction term calculation;
     *
     * if \f$X=\frac{\ln\frac{\gamma}{\beta}}{\ln10} < X_0\f$  :
     * \f$\delta= \delta_0 10^{2(X-X_0)}\f$
     *
     * if \f$ X<X_1\f$ :  \f$\delta= c_1X+c+a(X1-X)^m\f$
     * either :  \f$\delta = c_1x+c  \f$
     *
     * \return density correction term
     */

    double Delta(double beta, double gamma);

//----------------------------------------------------------------------------//

    /*!
     * this is what d2Ndvdx is equal to in the first approximation:
     * \f[\frac{d^2N}{dxdv}=\frac{Kz^2Z\rho}{2A\beta^2Ev^2}
     * \Big[1-\beta^2 \frac{v}{v_{max}}+\frac{1}{2}\Big
     * (\frac{v}{1+\frac{1}{\gamma}}\Big)^2 \Big]\f]
     *
     * \param   v   relative energy loss
     * \return  d2N/dvdx
     */

    double D2Ndvdx(double energy, double v);

//----------------------------------------------------------------------------//
    /*!
     * this is the inelastic electron bremsstrahlung correction to
     * the first approximation of the \f$ \frac{d2N}{dvdx}\f$:
     * \f[\Delta_{inelastic} =\frac{\alpha}{2\pi} \Big[
     * \ln\big(1+\frac{2vE}{m_e}\big)\Big[ 2\ln\Big(\frac{1-\frac{v}{v_{max}}}
     * {1-v}\Big)+\ln\Big(\frac{2\gamma m_e(1-v)}{mv}\Big)\Big]-\ln^2
     * \Big(\frac{1-\frac{v}{v_{max}}}{1-v}\Big) \Big]\f]
     *
     * \param   v   relative energy loss
     * \return  \f[\Delta_{inelastic}\f]
     */

    double InelCorrection(double energy, double v);

//----------------------------------------------------------------------------//

    double FunctionToBuildDEdxInterpolant( double energy);

//----------------------------------------------------------------------------//

    double FunctionToBuildDNdxInterpolant(double energy);

//----------------------------------------------------------------------------//

    double FunctionToBuildDNdxInterpolant2D( double energy , double v);

//----------------------------------------------------------------------------//

    double CalculateStochasticLoss(double energy, double rnd1);

//----------------------------------------------------------------------------//

public:

//----------------------------------------------------------------------------//

    Ionization();
    Ionization(const Ionization&);
    Ionization(PROPOSALParticle&, Medium* medium, EnergyCutSettings* cut_settings);
    Ionization& operator=(const Ionization&);
    bool operator==(const Ionization &ioniz) const;
    bool operator!=(const Ionization &ioniz) const;
    friend std::ostream& operator<<(std::ostream& os, Ionization const &ioniz);

//----------------------------------------------------------------------------//
    void swap(Ionization &ioniz);

//----------------------------------------------------------------------------//

    /*!
     * \brief contribution of ionization to -dE/dx.
     *
     * \f[\frac{dE'}{dx}=c_i\Big[\rho \cdot \frac{dE}{dx} +
     * E \int_{v_{min}}^{v_{up}} v \frac{d^2N}{dvdx}\Delta_{inelastic}dv\Big]\f]
     * with
     * \f[ \frac{dE}{dx} =Kz^2\frac{Z}{A}\frac{1}{2\beta^2}\Big[ \ln(2v_{up}m_e E)
     * +2\ln\Big(\frac{\beta \gamma}{10^{-6}I}\Big)+\Big(
     * \frac{v_{up}}{2\Big(1+\frac{1}{\gamma}\Big)}
     * \Big)^2-\beta^2\Big(1+\frac{v_{up}}{v_{Max}}\Big)-\delta \Big] \f]
     *
     * \return dE/dx [E/cm]
     */

    double CalculatedEdx(double energy);

//----------------------------------------------------------------------------//

    double CalculatedNdx(double energy);


//----------------------------------------------------------------------------//

    double CalculatedNdx(double energy, double rnd);

//----------------------------------------------------------------------------//

    double CalculateStochasticLoss(double energy, double rnd1, double rnd2);

//----------------------------------------------------------------------------//

    void EnableDNdxInterpolation( std::string path ="", bool raw=false);

//----------------------------------------------------------------------------//

    void EnableDEdxInterpolation( std::string path ="", bool raw=false);

//----------------------------------------------------------------------------//

    void DisableDNdxInterpolation();

//----------------------------------------------------------------------------//

    void DisableDEdxInterpolation();

//----------------------------------------------------------------------------//

    double FunctionToDNdxIntegral(double energy, double variable);

//----------------------------------------------------------------------------//

    void ValidateOptions();

//----------------------------------------------------------------------------//
//-----------------------------Getter and Setter------------------------------//
//----------------------------------------------------------------------------//
    //Getter

	const Interpolant* GetDedxInterpolant() const {
		return dedx_interpolant_;
	}

	const Interpolant* GetDndxInterpolant1d() const {
		return dndx_interpolant_1d_;
	}

	const Interpolant* GetDndxInterpolant2d() const {
		return dndx_interpolant_2d_;
	}

	const Integral* GetIntegral() const {
		return integral_;
	}

//----------------------------------------------------------------------------//
    //Setter

    // void SetParametrization(ParametrizationType::Enum parametrization = ParametrizationType::IonizBetheBloch);


//----------------------------------------------------------------------------//
    //Destructor
    ~Ionization(){}

};

}

#endif //Ionization_H
