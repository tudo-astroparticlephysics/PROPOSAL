#include "PROPOSAL/CrossSections.h"
#include "PROPOSAL/Integral.h"


class Epairproduction: public CrossSections
{
protected:

    int component_;
    double v_;
    bool reverse_;
    double eLpm_;


    Integral*   integral_;
    Integral*   integral_for_dEdx_;

//----------------------------------------------------------------------------//

    /*!
    this is the calculation of the d2Sigma/dvdRo - interface to Integral,
    the function which is defined here is:
    \f[ f(r) =return= \alpha^2r_e^2 \frac{2Z}{1,5\pi}(Z+k)
    \frac{1-v}{v}lpm(r^2,b,s)(F_e+\frac{m_e^2}{m_p^2}F_m)\f]
    */

    double FunctionToIntegral(double variable);

//----------------------------------------------------------------------------//

    /*!
    Landau Pomeranchuk Migdal effect evaluation,
    if the Landau Pomeranchuk Migdal effect is considered in
    the calculation, function is modified by a factor
    \f[lpm=return=\frac{(1+b)(A+(1+r^2)B)+b(C+(1+r^2)D)+(1-r^2)E}{[(2+r^2)(1+b)
    +x(3+r^2)]\ln\Big(1+\frac{1}{x}\Big)+\frac{1-r^2-b}{1+x}-(3+r^2)}\f]
    */

    double lpm(double r2, double b, double x);

//----------------------------------------------------------------------------//

    /*!
    this is the calculation of the dSigma/dv:
    \f[e_{Pair}=return = \rho N_Z z^2 \Big[ \int_{1-r_{max}}^{aux}
    f(r)dr + \int^1_{aux} f(r) dr \Big]\f]
    with \f$ aux=max(1-r_{Max}, ComputerPrecision)\f$ and \f$r_{max} =
    \sqrt{1-\frac{4m_e}{E_p v}}\Big(1-\frac{6m_p^2}{E_p^2(1-v)}\Big)\f$
    */

    double EPair(double v, int component);

//----------------------------------------------------------------------------//

    /*!
    pair production energy losses - interface to Integral;
    \f[ return= v\cdot e_{Pair}(v, components)\f]
    */

    double FunctionToContinuousIntegral(double variable);

//----------------------------------------------------------------------------//

    double FunctionToStochasticalIntegral(double variable);

//----------------------------------------------------------------------------//

    /*!
    Set the component of the primary;
    sets the parameters \f$v_{Min}=4 \frac{m_e}{E_p}\f$,
    \f$v_{Max}= max(-\frac{3}{4}\sqrt{e}\frac{m_p}{E_p}Z^{\frac{1}{3}},
    1-6\Big(\frac{m_p}{E_P}\Big)^2, 1-\frac{m_p}{E_p} )\f$
    and \f$v_{Up}=max( v_{max} , v_{Cut}(E_p))\f$
    */

    void SetIntegralLimits(int component);

//----------------------------------------------------------------------------//

public:

//----------------------------------------------------------------------------//

    Epairproduction();
    Epairproduction(const Epairproduction&);
    Epairproduction& operator=(const Epairproduction&);
    Epairproduction(Particle* particle, Medium* medium, EnergyCutSettings* cut_settings);

//----------------------------------------------------------------------------//

    /*!
    contribution of pair production to -dE/dx,
    \f[\frac{dE}{dx}=return= multiplier \cdot E_p \sum_{numC}
    \int_{v_{Min}}^{v_{Up}} v e_{Pair}(v, components) dv\f]
    \return [E/cm]
    */

    double CalculatedEdx();

//----------------------------------------------------------------------------//

    double CalculatedNdx();

//----------------------------------------------------------------------------//

    double CalculatedNdx(double rnd);

//----------------------------------------------------------------------------//

    double CalculateStochasticLoss();

//----------------------------------------------------------------------------//

    void EnableStochasticInerpolation();

//----------------------------------------------------------------------------//

    void EnableContinuousInerpolation();

//----------------------------------------------------------------------------//




    //Setter


    ~Epairproduction(){}

};
