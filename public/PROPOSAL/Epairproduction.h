#include "CrossSections.h"

class Epairproduction: public CrossSections
{
protected:


public:
//----------------------------------------------------------------------------//

    Epairproduction();
    Epairproduction(const Epairproduction&);
    Epairproduction& operator=(const Epairproduction&);

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

    ~Epairproduction(){}

};
