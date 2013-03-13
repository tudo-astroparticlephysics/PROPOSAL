#include "PROPOSAL/CrossSections.h"
#include "PROPOSAL/Constants.h"


CrossSections::CrossSections()
    :vMax_(0)
    ,vUp_ (0)
    ,vMin_(0)
    ,ebig_(BIGENERGY)
    ,doContinuousInterpolation_(false)
    ,doStochasticInterpolation_(false)
    ,multiplier_(1.)
    ,parametrization_(1)
{

}

//----------------------------------------------------------------------------//

void CrossSections::SetParametrizationLimit(double ebig){
    ebig_=ebig;
}

//----------------------------------------------------------------------------//

void CrossSections::SetMultiplier(double multiplier){
    multiplier_=multiplier;
}

//----------------------------------------------------------------------------//

void CrossSections::SetVMin(double vMin){
    vMin_=vMin;
}

//----------------------------------------------------------------------------//

void CrossSections::SetVMax(double vMax){
    vMax_=vMax;
}

//----------------------------------------------------------------------------//

void CrossSections::SetVUp(double vUp){
    vUp_=vUp;
}

//----------------------------------------------------------------------------//

void CrossSections::SetParametrization(int parametrization){
    parametrization_=parametrization;
}

//----------------------------------------------------------------------------//
