#ifndef Photonuclear_H
#define Photonuclear_H


#include "CrossSections.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Interpolant.h"
#include <vector>

class Photonuclear: public CrossSections
{
protected:

    int         component_;     /// nucleon in the medium on which the bremsstahlung occur
    bool        init_measured_;
    bool        init_hardbb_;
    int         hmax_;
    double      v_;
    bool        do_photo_interpolation_;   
    int         shadow_;
    bool        hard_component_;
    int         parametrization_family_;

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


//----------------------------------------------------------------------------//
    // Now: parametrization_ = 1,  Former: form=1 and bb=1 Kokoulin
    // Now: parametrization_ = 2,  Former: form=2 and bb=1 Kokoulin + hard component
    // Now: parametrization_ = 3,  Former: form=1 and bb=2 Rhode
    // Now: parametrization_ = 4,  Former: form=2 and bb=2 Rhode + hard component
    // Now: parametrization_ = 5,  Former: form=1 and bb=3 Bezrukov/Bugaev
    // Now: parametrization_ = 6,  Former: form=2 and bb=3 Bezrukov/Bugaev + hard component
    // Now: parametrization_ = 7,  Former: form=1 and bb=4 Zeus
    // Now: parametrization_ = 8,  Former: form=2 and bb=4 Zeus + hard component
    // Now: parametrization_ = 9,  Former: form=3 and bb=1 shadow=1 ALLM 91
    // Now: parametrization_ = 10, Former: form=3 and bb=1 shadow=2 ALLM 91
    // Now: parametrization_ = 11, Former: form=3 and bb=2 shadow=1 ALLM 97
    // Now: parametrization_ = 12, Former: form=3 and bb=2 shadow=2 ALLM 97
    // Now: parametrization_ = 13, Former: form=4 and bb=1 shadow=1 Butkevich/Mikhailov
    // Now: parametrization_ = 14, Former: form=4 and bb=1 shadow=2 Butkevich/Mikhailov

    double KokoulinParametrization(double v, int i);
    double RhodeParametrization(double v, int i);
    double BezrukovBugaevParametrization(double v, int i);
    double ZeusParametrization(double v, int i);
    double ALLM91Parametrization(double v, int i);
    double ALLM97Parametrization(double v, int i);
    double ButkevichMikhailovParametrization(double v, int i);

    double ShadowEffect(double x , double nu);


//----------------------------------------------------------------------------//

    /*!
     * parametrized photonuclear cross section - interface to Integral
     *
     * \param   Q2  square of the 4-momentum
     * \return  function value
     */

    double FunctionToIntegralALLM91(double Q2);
    double FunctionToIntegralALLM97(double Q2);
    double FunctionToIntegralButMik(double Q2);


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
    void SetIntegralLimits(int component);

//----------------------------------------------------------------------------//


    double FunctionToDEdxIntegral(double variable);

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

    double PhotoN(double v, int i);


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

    double HardBB(double e, double v);

//----------------------------------------------------------------------------//
    double FunctionToBuildDEdxInterpolant(double energy);

//----------------------------------------------------------------------------//

    double FunctionToBuildPhotoInterpolant(double energy , double v);

//----------------------------------------------------------------------------//
    double FunctionToBuildDNdxInterpolant1D(double energy);

//----------------------------------------------------------------------------//
    double FunctionToBuildDNdxInterpolant2D(double energy, double v);

//----------------------------------------------------------------------------//
    double CalculateStochasticLoss(double rnd1);

//----------------------------------------------------------------------------//

public:

//----------------------------------------------------------------------------//

    Photonuclear();
    Photonuclear(const Photonuclear&);
    Photonuclear& operator=(const Photonuclear&);
    bool operator==(const Photonuclear &photo) const;
    bool operator!=(const Photonuclear &photo) const;
    Photonuclear(Particle* particle, Medium* medium, EnergyCutSettings* cut_settings);
    friend std::ostream& operator<<(std::ostream& os, Photonuclear const &photo);

//----------------------------------------------------------------------------//
    void swap(Photonuclear &photo);

//----------------------------------------------------------------------------//
    double CalculatedEdx();

//----------------------------------------------------------------------------//

    double CalculatedNdx();

//----------------------------------------------------------------------------//

    double CalculatedNdx(double rnd);

//----------------------------------------------------------------------------//

    double CalculateStochasticLoss(double rnd1, double rnd2);

//----------------------------------------------------------------------------//

    void EnableDNdxInterpolation(std::string path ="", bool raw=false);

//----------------------------------------------------------------------------//

    void EnableDEdxInterpolation(std::string path ="", bool raw=false);

//----------------------------------------------------------------------------//

    void EnablePhotoInterpolation(std::string path ="", bool raw=false);

//----------------------------------------------------------------------------//

    void DisableDNdxInterpolation();

//----------------------------------------------------------------------------//

    void DisableDEdxInterpolation();

//----------------------------------------------------------------------------//
    void DisablePhotoInterpolation();

//----------------------------------------------------------------------------//

    double FunctionToDNdxIntegral(double variable);

//----------------------------------------------------------------------------//

    boost::program_options::options_description CreateOptions();

//----------------------------------------------------------------------------//

    void ValidateOptions();

//----------------------------------------------------------------------------//
    //Getter

	int GetComponent() const {
		return component_;
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

	int GetHmax() const {
		return hmax_;
	}

	bool GetInitHardbb() const {
		return init_hardbb_;
	}

	bool GetInitMeasured() const {
		return init_measured_;
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

	double GetV() const {
		return v_;
	}

    int GetShadow() const {
        return shadow_;
    }

    bool GetHardComponent() const {
        return hard_component_;
    }
//----------------------------------------------------------------------------//
    //Setter
    void SetParametrization(int parametrization=1);
	void SetComponent(int component);
	void SetDedxInterpolant(Interpolant* dedxInterpolant);
	void SetDndxIntegral(std::vector<Integral*> dndxIntegral);
	void SetDndxInterpolant1d(std::vector<Interpolant*> dndxInterpolant1d);
	void SetDndxInterpolant2d(std::vector<Interpolant*> dndxInterpolant2d);
	void SetHmax(int hmax);
	void SetInitHardbb(bool initHardbb);
	void SetInitMeasured(bool initMeasured);
	void SetIntegral(Integral* integral);
	void SetIntegralForDEdx(Integral* integralForDEdx);
	void SetInterpolantHardBb(std::vector<Interpolant*> interpolantHardBb);
	void SetInterpolantMeasured(Interpolant* interpolantMeasured);
	void SetProbForComponent(std::vector<double> probForComponent);
    void SetV(double v);
    void SetHardComponent(bool hard);
    void SetShadow(int shadow);

//----------------------------------------------------------------------------//
    //Destructor
    ~Photonuclear(){}


};

#endif //Photonuclear_H
