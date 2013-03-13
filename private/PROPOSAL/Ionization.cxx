#include "PROPOSAL/Ionization.h"

Ionization::Ionization(){

}
//----------------------------------------------------------------------------//

Ionization::Ionization(const Ionization &ioniz)
{
    *this = ioniz;
}
//----------------------------------------------------------------------------//

Ionization& Ionization::operator=(const Ionization &ioniz){
    return *this;
}

//----------------------------------------------------------------------------//

void Ionization::SetIntegralLimits(){
}

//----------------------------------------------------------------------------//

double Ionization::CalculatedEdx(){
    return 0;
}
//----------------------------------------------------------------------------//

double Ionization::CalculatedNdx(){
    return 0;
}
//----------------------------------------------------------------------------//

double Ionization::CalculatedNdx(double rnd){
    return 0;
}
//----------------------------------------------------------------------------//

double Ionization::CalculateStochasticLoss(){
    return 0;
}
//----------------------------------------------------------------------------//

void Ionization::EnableStochasticInerpolation(){
    doStochasticInterpolation_=true;
}
//----------------------------------------------------------------------------//

void Ionization::EnableContinuousInerpolation(){
    doContinuousInterpolation_=true;
}

//----------------------------------------------------------------------------//

double Ionization::FunctionToContinuousIntegral(double variable){
    return 0;
}

//----------------------------------------------------------------------------//

double Ionization::FunctionToStochasticalIntegral(double variable){
    return 0;
}
//----------------------------------------------------------------------------//
