#include "CrossSections.h"

class Medium;
class PROPOSALParticle;

class Bremsstrahlung: protected CrossSections
{
protected:

    //Integrallimits
    double vMax_;
    double vUp_;
    double vMin_;

    //Parametrization
    int form_;

    bool init_;
    double eLpm_;
    double xo_;

    bool lorenz_;
    double lorenzCut_;  	/// in [MeV] // - set to 1.e6 in Constructor

    PROPOSALParticle*   particle_;
    Medium*             medium_;

    // Interpolation flags
    bool doContinuousInterpolation_;
    bool doStochasticInterpolation_;


public:
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

};
