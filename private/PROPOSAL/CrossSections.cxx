#include "PROPOSAL/CrossSections.h"
#include "PROPOSAL/Constants.h"


CrossSections::CrossSections(){
    elow_=0;                //!< lower bound of parameterizations
    nlow_=ME;               //!< maximal number of parametrization points
    ebig_=BIGENERGY;        //!< upper bound of parameterizations


    // Interpolation flags
    doContinuousInterpolation_ = false;
    doStochasticInterpolation_ = false;

    //CrossSection multiplier
    multiplier_ = 1.;
}

void CrossSections::SetParameterizationLimits(double elow,
                                              double nlow,
                                              double ebig){

    elow_=elow;
    nlow_=nlow;
    ebig_=ebig;
}

void CrossSections::SetMultiplier(double multiplier){

    multiplier_=multiplier;

}

