/*
 * CrossSections.h
 *
 *  Created on: 2013.03.12
 *      Author: koehne
 */

#ifndef CrossSections_H
#define CrossSections_H

#include <vector>
#include "Integral.h"



/*! \class CrossSections CrossSections.h "CrossSections.h"
    \brief This is a pure virtual class
 */



class CrossSections
{


protected:


    // bounds of parameterizations


public:

    CrossSections();


    double elow_;
    double nlow_;
    double ebig_;

    // Memberfunctions

//----------------------------------------------------------------------------//

    virtual void SetIntegralLimits() = 0;

//----------------------------------------------------------------------------//

    virtual double CalculatedEdx() = 0;

//----------------------------------------------------------------------------//

    virtual double CalculatedNdx() = 0;


//----------------------------------------------------------------------------//

    virtual double CalculatedNdx(double rnd) = 0;

//----------------------------------------------------------------------------//

    virtual double CalculateStochasticLoss() = 0;

//----------------------------------------------------------------------------//

    virtual void EnableStochasticInerpolation() = 0;

//----------------------------------------------------------------------------//

    virtual void EnableContinuousInerpolation() = 0;

//----------------------------------------------------------------------------//

    virtual double FunctionToContinuousIntegral(double integrand) = 0;

//----------------------------------------------------------------------------//

    virtual double FunctionToStochasticalIntegral(double integrand) = 0;


//----------------------------------------------------------------------------//

    void SetParameterizationLimits(double elow=0.,
                                   double nlow=ME,
                                   double ebig=BIGENERGY);

//----------------------------------------------------------------------------//



    // Setter


    
//----------------------------------------------------------------------------//

    // destructors

    ///@brief Crush this CrossSections.
    virtual ~CrossSections(){}


};



#endif //CrossSections_H
