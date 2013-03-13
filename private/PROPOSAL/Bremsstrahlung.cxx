#include "PROPOSAL/Bremsstrahlung.h"

Bremsstrahlung::Bremsstrahlung(){

}
//----------------------------------------------------------------------------//

Bremsstrahlung::Bremsstrahlung(const Bremsstrahlung &model)
{
    *this = model;
}
//----------------------------------------------------------------------------//

Bremsstrahlung& Bremsstrahlung::operator=(const Bremsstrahlung &model){
    return *this;
}

//----------------------------------------------------------------------------//

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

