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

    Integral*   integral_;
    Integral*   integral_for_dEdx_;
    Interpolant* dedx_interpolant_; //formerly interpolateJ_

    Interpolant* interpolant_measured_;
    std::vector<Interpolant*> interpolant_hardBB_;

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


    // Now: parametrization_ = 1,hard_component_= false / Former: form=1 and bb=1 Kokoulin
    // Now: parametrization_ = 1,hard_component_= true  / Former: form=2 and bb=1 Kokoulin + hard component
    // Now: parametrization_ = 2,hard_component_= false / Former: form=1 and bb=2 Rhode
    // Now: parametrization_ = 2,hard_component_= true  / Former: form=2 and bb=2 Rhode + hard component
    // Now: parametrization_ = 3,hard_component_= false / Former: form=1 and bb=3 Bezrukov/Bugaev
    // Now: parametrization_ = 3,hard_component_= true  / Former: form=2 and bb=3 Bezrukov/Bugaev + hard component
    // Now: parametrization_ = 4,hard_component_= false / Former: form=1 and bb=4 Zeus
    // Now: parametrization_ = 4,hard_component_= true  / Former: form=2 and bb=4 Zeus + hard component
    // Now: parametrization_ = 5,shadow_ = 1 or = 2     / Former: form=3 and bb=1 ALLM 91
    // Now: parametrization_ = 6,shadow_ = 1 or = 2     / Former: form=3 and bb=2 ALLM 97
    // Now: parametrization_ = 7,shadow_ = 1 or = 2     / Former: form=4 and bb=1 Butkevich/Mikhailov
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

    double FunctionToDNdxIntegral(double variable);

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

public:

//----------------------------------------------------------------------------//

    Photonuclear();
    Photonuclear(const Photonuclear&);
    Photonuclear& operator=(const Photonuclear&);
    Photonuclear(Particle* particle, Medium* medium, EnergyCutSettings* cut_settings);



    double CalculatedEdx();

//----------------------------------------------------------------------------//

    double CalculatedNdx();


//----------------------------------------------------------------------------//

    double CalculatedNdx(double rnd);

//----------------------------------------------------------------------------//

    double CalculateStochasticLoss();

//----------------------------------------------------------------------------//

    double CalculateStochasticLoss(double rnd1, double rnd2);

//----------------------------------------------------------------------------//

    void EnableDNdxInterpolation();

//----------------------------------------------------------------------------//

    void EnableDEdxInterpolation();

//----------------------------------------------------------------------------//

    double FunctionToBuildDEdxInterpolant(double energy);

//----------------------------------------------------------------------------//

    void DisableDNdxInterpolation();

//----------------------------------------------------------------------------//

    void DisableDEdxInterpolation();

//----------------------------------------------------------------------------//

    ~Photonuclear(){}

};
