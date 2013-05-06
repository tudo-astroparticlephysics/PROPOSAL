#include "PROPOSAL/Decay.h"

Decay::Decay()
    :CrossSections ( )

{
    name_   = "Decay";

}
//----------------------------------------------------------------------------//

Decay::Decay(const Decay &decay)
    :CrossSections  ( decay )
{
}
//----------------------------------------------------------------------------//
Decay::Decay(Particle* particle, Medium* medium, EnergyCutSettings* cut_settings)
    :CrossSections          ( particle, medium, cut_settings )

{
    name_   = "Decay";

}
//----------------------------------------------------------------------------//

Decay& Decay::operator=(const Decay &decay){
    if (this != &decay)
    {
      Decay tmp(decay);
      swap(tmp);
    }
    return *this;
}


//----------------------------------------------------------------------------//

bool Decay::operator==(const Decay &decay) const
{
    if( this->CrossSections::operator !=(decay) )   return false;

    //else
    return true;
}

//----------------------------------------------------------------------------//

bool Decay::operator!=(const Decay &decay) const
{
    return !(*this == decay);
}
//----------------------------------------------------------------------------//

void Decay::swap(Decay &decay)
{
    using std::swap;

    this->CrossSections::swap(decay);

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

double Decay::CalculateStochasticLoss(double rnd1, double rnd2){
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

void Decay::DisableDNdxInterpolation(){
    do_dndx_Interpolation_=false;

}

//----------------------------------------------------------------------------//

void Decay::DisableDEdxInterpolation(){
    do_dedx_Interpolation_=false;
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
