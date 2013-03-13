#include <vector>
#include "CrossSections.h"

class Medium;
class PROPOSALParticle;

class Bremsstrahlung: public CrossSections
{
protected:

    bool init_;
    double eLpm_;
    double xo_;

    bool lorenz_;
    double lorenzCut_;  	/// in [MeV] // - set to 1.e6 in Constructor

    PROPOSALParticle*   particle_;
    Medium*             medium_;


    std::vector<Integral*>   integrals_;

public:
//----------------------------------------------------------------------------//

    Bremsstrahlung();
    Bremsstrahlung(const Bremsstrahlung&);
    Bremsstrahlung& operator=(const Bremsstrahlung&);

//----------------------------------------------------------------------------//

    void SetIntegralLimits();

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

    ~Bremsstrahlung(){}

};
