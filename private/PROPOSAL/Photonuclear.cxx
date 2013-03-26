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


double Photonuclear::PhotoN(double v, int i){
    return 0;
}

//----------------------------------------------------------------------------//


void Photonuclear::SetMeasured(){

}

//----------------------------------------------------------------------------//

double Photonuclear::MeasuredSgN(double e){
    return 0;
}

//----------------------------------------------------------------------------//

void Photonuclear::EnableHardBB(){

}

//----------------------------------------------------------------------------//

double Photonuclear::HardBB(double e, double v){
    return 0;
}

//----------------------------------------------------------------------------//

double Photonuclear::FunctionToIntegral(double Q2){
    return 0;
}

//----------------------------------------------------------------------------//
