#include "CrossSections.h"

class Decay: public CrossSections
{
protected:


public:
//----------------------------------------------------------------------------//

    Decay();
    Decay(const Decay&);
    Decay& operator=(const Decay& decay);
    Decay(Particle* particle, Medium* medium, EnergyCutSettings* cut_settings);
    bool operator==(const Decay &decay) const;
    bool operator!=(const Decay &decay) const;

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

    double CalculateStochasticLoss(double rnd1, double rnd2);
//----------------------------------------------------------------------------//

    void EnableDNdxInterpolation();

//----------------------------------------------------------------------------//

    void EnableDEdxInterpolation();

//----------------------------------------------------------------------------//

    void DisableDNdxInterpolation();

//----------------------------------------------------------------------------//

    void DisableDEdxInterpolation();

//----------------------------------------------------------------------------//

    double FunctionToDEdxIntegral(double variable);

//----------------------------------------------------------------------------//

    double FunctionToDNdxIntegral(double variable);

//----------------------------------------------------------------------------//

    void swap(Decay &decay);

//----------------------------------------------------------------------------//
    ~Decay(){}

    //Getter
//----------------------------------------------------------------------------//


};
