#include "PROPOSAL/CrossSections.h"
#include "PROPOSAL/Constants.h"


CrossSections::CrossSections(){
    elow_=0;                //!< lower bound of parameterizations
    nlow_=ME;               //!< maximal number of parametrization points
    ebig_=BIGENERGY;        //!< upper bound of parameterizations
}

void CrossSections::SetParameterizationLimits(double elow,
                                              double nlow,
                                              double ebig){

    elow_=elow;
    nlow_=nlow;
    ebig_=ebig;
}
