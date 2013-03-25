#include "PROPOSAL/CrossSections.h"
#include "PROPOSAL/Integral.h"


class Epairproduction: public CrossSections
{
protected:

    int component_;
    double v_;
    bool reverse_;

    Integral*   integral_;
    Integral*   integral_for_dEdx_;



public:
//----------------------------------------------------------------------------//

    Epairproduction();
    Epairproduction(const Epairproduction&);
    Epairproduction& operator=(const Epairproduction&);
    Epairproduction(Particle* particle, Medium* medium, EnergyCutSettings* cut_settings);

//----------------------------------------------------------------------------//

    void SetIntegralLimits(int component);

//----------------------------------------------------------------------------//

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

    double FunctionToContinuousIntegral(double variable);

//----------------------------------------------------------------------------//

    double FunctionToStochasticalIntegral(double variable);

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

    //Setter


    ~Epairproduction(){}

};
