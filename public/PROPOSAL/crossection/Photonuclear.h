#pragma once

// #include <vector>

#include "PROPOSAL/crossection/CrossSections.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/Interpolant.h"

namespace PROPOSAL
{
    class Photonuclear;
}

std::ostream& operator<<(std::ostream& os, PROPOSAL::Photonuclear const &photo);

namespace PROPOSAL{

namespace ShadowingType
{
    enum Enum
    {
        Dutta               = 141,
        ButkevichMikhailov  = 142
    };
}

class Photonuclear: public CrossSections
{
protected:

    int         component_;     /// nucleon in the medium on which the bremsstahlung occur
    bool        init_measured_;
    bool        init_hardbb_;
    int         hmax_;
    double      v_;
    bool        do_photo_interpolation_;
    bool        hard_component_;
    ShadowingType::Enum        shadow_;

    Integral*   integral_;
    Integral*   integral_for_dEdx_;
    std::vector<Integral*>    dndx_integral_;

    Interpolant* dedx_interpolant_; //formerly interpolateJ_

    Interpolant* interpolant_measured_;
    std::vector<Interpolant*> interpolant_hardBB_;
    std::vector<Interpolant*> dndx_interpolant_1d_; //Stochastic dNdx()
    std::vector<Interpolant*> dndx_interpolant_2d_; //Stochastic dNdx()
    std::vector<Interpolant*> photo_interpolant_;  //!< Interpolates function used by dNdx and dEdx

    std::vector<double> prob_for_component_; //!< probability for each medium component to interact with the particle (formerly H_)

    std::vector<double> woodSaxonPotential_; //!< Wood-Saxon potential factor (needed for Butkevich Shadowing)


//----------------------------------------------------------------------------//

    /*!
     * parameterized inelastic nuclear scattering
     * with two different approaches:
     *
     * 1. to assume that there is an interaction between a real photon with a nucleus
     * Parameterizations using this strategy:
     * Kokoulin and Petrukhin in Proc. Vulcano Workshop 1996
     * Rhode in Nuclear Physics B (Proc. Suppl.) 35 1994
     * Bezrukov and Bugaev in Sov.J.Nucl.Phys. 33 1981
     * Derrick et al (ZEUS Collaboration) in Z.Phys. C 63 1994
     *
     * 2. to integrate over the momentum squared of the exchange photon Q2
     * Parameterizations using this strategy:
     * Abramowicz et al in Phys.Let. B 269 1991
     * Abramowicz et al in DESY Reports 251 1997
     * Butkevich and Mikheyev in JETP 95 2002
     * Reno and Sarcevic and Su in Astrop. Phys. 24 2005
     *
     * \param   v   relative energy loss
     * \param   i   crossection component
     * \return  photonuclear cross section
     */


    double ParametrizationOfRealPhotonAssumption(double energy, double v, int i);
    double ParametrizationOfQ2Integration(double energy, double v, int i);

//----------------------------------------------------------------------------//

    /*!
     * shadowing effect parametrized by Dutta and Butkevich/Mikhailov
     * for Parametrization with Q2 Integration
     * and Bezrukov/Bugaev parametrization
     * for Parametrization with real photon assumption
     *
     * \return  shadowing factor
     */

    double ShadowEffect(double x , double nu);
    double ShadowBezrukovBugaev(double sgn, double atomic_number);

//----------------------------------------------------------------------------//

    /*!
     * parametrized cross section of the interaction
     * between a nucleus and a real photon
     *
     * a parametrization is used only for high energy transfer
     * Caldwell in Phys.Rev.Let. 42 1979
     *
     * \param   nu  Energy of thr "real" photon
     * \return  function value
     */

    double PhotoNucleusCrossSectionCaldwell(double nu);

    double PhotoNucleusCrossSectionKokoulin(double nu);
    double PhotoNucleusCrossSectionRhode(double nu);
    double PhotoNucleusCrossSectionBezrukovBugaev(double nu);
    double PhotoNucleusCrossSectionZeus(double nu, double medium_average_nucleon_weight);

//----------------------------------------------------------------------------//

    /*!
     * parametrized photonuclear cross section - interface to Integral
     *
     * \param   Q2  square of the 4-momentum
     * \return  function value
     */

    double FunctionToIntegralALLM91(double energy, double Q2);
    double FunctionToIntegralALLM97(double energy, double Q2);
    double FunctionToIntegralButMik(double energy, double Q2);
    double FunctionToIntegralRSS(double energy, double Q2);


//----------------------------------------------------------------------------//

    /*!
     * Set the component of the primary
     * sets the parameters \f$v_{min}=\frac{1}{E_{p}}\Big(m_{\pi}+
     * \frac{m_{\pi}^2}{2M_N}\Big)\f$, \f$v_{Max}=min\Big(v_{Max},
     * 1-\frac{m_p}{E_p}\Big)\f$
     * with \f$v_{Max}=1-\frac{M_N\Big[1+\Big(\frac{m_p}{M_N}\Big)^2\Big]}
     * {2E_p}\f$ if \f$m_p<m_{\pi}\f$ otherwise
     * \f$ v_{Max}=1\f$;
     * \f$v_{up}=min(v_{Max}, v_{Cut}(E))\f$
     *
     * \param   i   crossection component
     */
    CrossSections::IntegralLimits SetIntegralLimits(double energy, int component);

//----------------------------------------------------------------------------//


    double FunctionToDEdxIntegral(double energy, double variable);

//----------------------------------------------------------------------------//

    /*!
     * this is what the photonuclear interaction cross section is equal to,
     * there are different cases in which the cross section
     * \f$ \sigma_{\gamma,N}\f$ is calculated,
     * intermediate result is:
     * \f[f(v) = \frac{\alpha}{2\pi}A \sigma_{\gamma, N}v\Big[\frac{3}{4}G
     * \Big[k\ln\Big(1+\frac{m_1}{t}\Big)+ \frac{4\mu^2}{m_1}\ln
     * \Big(1+\frac{m_1}{t}\Big)-\frac{km_1}{m_1+t}-\frac{2\mu^2}{t} \Big] +
     * \frac{1}{4}\Big[\Big(k+\frac{2\mu^2}{m_2}\Big)\ln\Big(1+\frac{m_2}{t}\Big)
     * -\frac{2\mu^2}{t}\Big]
     * + \frac{1}{2}\frac{\mu^2}{t}\Big[\frac{3}{4}G\frac{m_1-4t}{m_1+t}+
     * \frac{1}{4}\frac{m_2}{t}\ln\Big(1+\frac{t}{m_2}\Big)\Big] \Big]\f]
     * the final result is:
     * \f[\rho_{mol}N(z^2)^2 \int_{min}^{max} \frac{d\sigma}{dv} dv\f]
     * with
     * \f[\frac{d\sigma}{dv} = \rho_{mol}N(z^2)^2\cdot f(v)\f]
     *
     * \param   v   relative energy loss
     * \param   i   crossection component
     * \return  \f$\rho_{mol}N(z^2)^2 \int_{min}^{max} \frac{d\sigma}{dv} dv\f$
     */

    double PhotoN(double energy, double v, int i);


//----------------------------------------------------------------------------//

    /*!
     * call this to replace photon-nucleon cross section
     * approximation formula with
     * a more exact at low energies parametrization of the experimental data
     *
     */

    void SetMeasured();

//----------------------------------------------------------------------------//
    /*!
     * \brief returns measured value of the photon nucleon cross section
     *
     * \param   e   energy
     * \return  value of the photon nucleon cross section
     */

    double MeasuredSgN(double e);

//----------------------------------------------------------------------------//
    /*!
     * \brief initializes hard part of BB (by Bugaev and Shlepin)
     * parametrization calculation
     */

    void EnableHardBB();

//----------------------------------------------------------------------------//
    /*!
     * \brief interpolated value of the hard part of bb cross section
     *
     * \param   e   energy
     * \param   v   relative energy loss
     * \return  hard part of bb cross section
     */

    double HardBB(double energy, double v);

//----------------------------------------------------------------------------//
    double FunctionToBuildDEdxInterpolant( double energy);

//----------------------------------------------------------------------------//

    double FunctionToBuildPhotoInterpolant( double energy , double v);

//----------------------------------------------------------------------------//
    double FunctionToBuildDNdxInterpolant1D(double energy);

//----------------------------------------------------------------------------//
    double FunctionToBuildDNdxInterpolant2D( double energy, double v);

//----------------------------------------------------------------------------//
    double CalculateStochasticLoss(double energy, double rnd1);

//----------------------------------------------------------------------------//

    void CalculateWoodSaxonPotential();

//----------------------------------------------------------------------------//

    double FunctionToWoodSaxonPotentialIntegral(const double r0, double r);

//----------------------------------------------------------------------------//

public:

//----------------------------------------------------------------------------//

    Photonuclear();
    Photonuclear(const Photonuclear&);
    Photonuclear(PROPOSALParticle&, Medium* medium, EnergyCutSettings* cut_settings);
    Photonuclear& operator=(const Photonuclear&);
    bool operator==(const Photonuclear &photo) const;
    bool operator!=(const Photonuclear &photo) const;
    friend std::ostream& operator<<(std::ostream& os, Photonuclear const &photo);

//----------------------------------------------------------------------------//
    void swap(Photonuclear &photo);

//----------------------------------------------------------------------------//
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

    void EnablePhotoInterpolation( std::string path ="", bool raw=false);

//----------------------------------------------------------------------------//

    void DisableDNdxInterpolation();

//----------------------------------------------------------------------------//

    void DisableDEdxInterpolation();

//----------------------------------------------------------------------------//
    void DisablePhotoInterpolation();

//----------------------------------------------------------------------------//

    double FunctionToDNdxIntegral(double energy, double variable);

//----------------------------------------------------------------------------//

    // void ValidateOptions();

//----------------------------------------------------------------------------//
//-----------------------------Getter and Setter------------------------------//
//----------------------------------------------------------------------------//
    //Getter

    int GetComponent() const { return component_; }

    int GetHmax() const { return hmax_; }

    bool GetInitHardbb() const { return init_hardbb_; }

    bool GetInitMeasured() const { return init_measured_; }

    double GetV() const {
        return v_;
    }

    ShadowingType::Enum GetShadow() const {
        return shadow_;
    }

    bool GetHardComponent() const {
        return hard_component_;
    }

	Interpolant* GetDedxInterpolant() const {
		return dedx_interpolant_;
	}

	std::vector<Integral*> GetDndxIntegral() const {
		return dndx_integral_;
	}

	std::vector<Interpolant*> GetDndxInterpolant1d() const {
		return dndx_interpolant_1d_;
	}

	std::vector<Interpolant*> GetDndxInterpolant2d() const {
		return dndx_interpolant_2d_;
	}

	Integral* GetIntegral() const {
		return integral_;
	}

	Integral* GetIntegralForDEdx() const {
		return integral_for_dEdx_;
	}

	std::vector<Interpolant*> GetInterpolantHardBb() const {
		return interpolant_hardBB_;
	}

	Interpolant* GetInterpolantMeasured() const {
		return interpolant_measured_;
	}

	std::vector<double> GetProbForComponent() const {
		return prob_for_component_;
	}

    std::vector<double> GetWoodSaxonPotential() const {
        return woodSaxonPotential_;
    }

//----------------------------------------------------------------------------//
    //Setter

    void SetComponent(int component){ component_ = component; }

	void SetHmax(int hmax){ hmax_ = hmax; }

	void SetInitHardbb(bool initHardbb){ init_hardbb_ = initHardbb; }

	void SetInitMeasured(bool initMeasured){ init_measured_ = initMeasured; }

    // void SetParametrization(ParametrizationType::Enum parametrization = ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich);

//----------------------------------------------------------------------------//
    //Destructor
    ~Photonuclear(){}


};

}
