#include "PROPOSAL/Photonuclear.h"

Photonuclear::Photonuclear(){

}
//----------------------------------------------------------------------------//

Photonuclear::Photonuclear(const Photonuclear &photo)
{
    *this = photo;
}
//----------------------------------------------------------------------------//

Photonuclear& Photonuclear::operator=(const Photonuclear &photo){
    return *this;
}

//----------------------------------------------------------------------------//

void Photonuclear::SetIntegralLimits(int component){
}

//----------------------------------------------------------------------------//

double Photonuclear::CalculatedEdx(){
    return 0;
}
//----------------------------------------------------------------------------//

double Photonuclear::CalculatedNdx(){
    return 0;
}
//----------------------------------------------------------------------------//

double Photonuclear::CalculatedNdx(double rnd){
    return 0;
}
//----------------------------------------------------------------------------//

double Photonuclear::CalculateStochasticLoss(){
    return 0;
}
//----------------------------------------------------------------------------//

void Photonuclear::EnableStochasticInerpolation(){
    doStochasticInterpolation_=true;
}
//----------------------------------------------------------------------------//

void Photonuclear::EnableContinuousInerpolation(){
    doContinuousInterpolation_=true;
}

//----------------------------------------------------------------------------//

double Photonuclear::FunctionToContinuousIntegral(double variable){
    return 0;
}

//----------------------------------------------------------------------------//

double Photonuclear::FunctionToStochasticalIntegral(double variable){
    return 0;
}
//----------------------------------------------------------------------------//
