#include "PROPOSAL/Decay.h"

Decay::Decay(){

}
//----------------------------------------------------------------------------//

Decay::Decay(const Decay &decay)
{
    *this = decay;
}
//----------------------------------------------------------------------------//

Decay& Decay::operator=(const Decay &decay){
    return *this;
}

//----------------------------------------------------------------------------//

void Decay::SetIntegralLimits(int component){
}

//----------------------------------------------------------------------------//

double Decay::CalculatedEdx(){
    return 0;
}
//----------------------------------------------------------------------------//

double Decay::CalculatedNdx(){
    return 0;
}
//----------------------------------------------------------------------------//

double Decay::CalculatedNdx(double rnd){
    return 0;
}
//----------------------------------------------------------------------------//

double Decay::CalculateStochasticLoss(){
    return 0;
}
//----------------------------------------------------------------------------//

void Decay::EnableDNdxInterpolation(){
    do_dndx_Interpolation_=true;
}
//----------------------------------------------------------------------------//

void Decay::EnableDEdxInterpolation(){
    do_dedx_Interpolation_=true;
}

//----------------------------------------------------------------------------//

double Decay::FunctionToDEdxIntegral(double variable){
    return 0;
}

//----------------------------------------------------------------------------//

double Decay::FunctionToDNdxIntegral(double variable){
    return 0;
}
//----------------------------------------------------------------------------//
