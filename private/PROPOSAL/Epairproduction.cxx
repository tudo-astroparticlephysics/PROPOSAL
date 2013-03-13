#include "PROPOSAL/Epairproduction.h"

Epairproduction::Epairproduction(){

}
//----------------------------------------------------------------------------//

Epairproduction::Epairproduction(const Epairproduction &epair)
{
    *this = epair;
}
//----------------------------------------------------------------------------//

Epairproduction& Epairproduction::operator=(const Epairproduction &epair){
    return *this;
}

//----------------------------------------------------------------------------//

void Epairproduction::SetIntegralLimits(){
}

//----------------------------------------------------------------------------//

double Epairproduction::CalculatedEdx(){
    return 0;
}
//----------------------------------------------------------------------------//

double Epairproduction::CalculatedNdx(){
    return 0;
}
//----------------------------------------------------------------------------//

double Epairproduction::CalculatedNdx(double rnd){
    return 0;
}
//----------------------------------------------------------------------------//

double Epairproduction::CalculateStochasticLoss(){
    return 0;
}
//----------------------------------------------------------------------------//

void Epairproduction::EnableStochasticInerpolation(){
    doStochasticInterpolation_=true;
}
//----------------------------------------------------------------------------//

void Epairproduction::EnableContinuousInerpolation(){
    doContinuousInterpolation_=true;
}

//----------------------------------------------------------------------------//

double Epairproduction::FunctionToContinuousIntegral(double variable){
    return 0;
}

//----------------------------------------------------------------------------//

double Epairproduction::FunctionToStochasticalIntegral(double variable){
    return 0;
}
//----------------------------------------------------------------------------//
