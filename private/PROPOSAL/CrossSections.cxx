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
    particle_       = new Particle();
    medium_         = new Medium();
    cut_settings_   = new EnergyCutSettings();
}
//----------------------------------------------------------------------------//

CrossSections::CrossSections(Particle* particle,
                             Medium* medium,
                             EnergyCutSettings* cut_settings)
    :vMax_(0)
    ,vUp_ (0)
    ,vMin_(0)
    ,ebig_(BIGENERGY)
    ,doContinuousInterpolation_(false)
    ,doStochasticInterpolation_(false)
    ,multiplier_(1.)
    ,parametrization_(1)
{
    particle_       = particle;
    medium_         = medium;
    cut_settings_   = cut_settings;
}
//----------------------------------------------------------------------------//

void CrossSections::SetParametrizationLimit(double ebig){
    ebig_ = ebig;
}

//----------------------------------------------------------------------------//

void CrossSections::SetMultiplier(double multiplier){
    multiplier_ = multiplier;
}

//----------------------------------------------------------------------------//

void CrossSections::SetVMin(double vMin){
    vMin_ = vMin;
}

//----------------------------------------------------------------------------//

void CrossSections::SetVMax(double vMax){
    vMax_ = vMax;
}

//----------------------------------------------------------------------------//

void CrossSections::SetVUp(double vUp){
    vUp_ = vUp;
}

//----------------------------------------------------------------------------//

void CrossSections::SetParametrization(int parametrization){
    parametrization_ = parametrization;
}

//----------------------------------------------------------------------------//

void CrossSections::SetParticle(Particle *particle){
    particle_ = particle;
}

//----------------------------------------------------------------------------//

void CrossSections::SetMedium(Medium *medium){
    medium_ = medium;
}
//----------------------------------------------------------------------------//

void CrossSections::SetEnergyCutSettings(EnergyCutSettings *cuts){
    cut_settings_ = cuts;
}
