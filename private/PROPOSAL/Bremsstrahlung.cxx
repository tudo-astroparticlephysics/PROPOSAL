#include "PROPOSAL/Bremsstrahlung.h"

void Bremsstrahlung::SetIntegralLimits(){
}

//----------------------------------------------------------------------------//

double Bremsstrahlung::CalculatedEdx(){
    return 0;
}
//----------------------------------------------------------------------------//

double Bremsstrahlung::CalculatedNdx(){
    return 0;
}
//----------------------------------------------------------------------------//

double Bremsstrahlung::CalculatedNdx(double rnd){
    return 0;
}
//----------------------------------------------------------------------------//

double Bremsstrahlung::CalculateStochasticLoss(){
    return 0;
}
//----------------------------------------------------------------------------//

void Bremsstrahlung::EnableStochasticInerpolation(){
    doStochasticInterpolation_=true;
}
//----------------------------------------------------------------------------//

void Bremsstrahlung::EnableContinuousInerpolation(){
    doContinuousInterpolation_=true;
}
