#include "CrossSections.h"


class Photonuclear: public CrossSections
{
protected:


public:
//----------------------------------------------------------------------------//

    Photonuclear();
    Photonuclear(const Photonuclear&);
    Photonuclear& operator=(const Photonuclear&);

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

    ~Photonuclear(){}

};
