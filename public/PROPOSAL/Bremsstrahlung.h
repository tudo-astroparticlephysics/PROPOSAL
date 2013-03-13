#include "CrossSections.h"

class Medium;
class PROPOSALParticle;

class Bremsstrahlung: protected CrossSections
{
protected:

    double vMax_;			// Set to 0 in Constructor
    double vUp_;			// Set to 0 in Constructor
    double vMin_;			// Set to 0 in Constructor

    int form_;			// Set to 1 in Constrcutor - Switches the parametrization which is used.

    bool init_;		// set to true in Constructor
    double eLpm_;
    double xo_;

    bool lorenz_;		// Set to false in constructor
    double lorenzCut_;  	/// in [MeV] // - set to 1.e6 in Constructor

    PROPOSALParticle*   particle_;  //Tomasz
    Medium*             medium_;    //Tomasz
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
};
