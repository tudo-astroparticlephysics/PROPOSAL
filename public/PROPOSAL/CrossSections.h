/*
 * CrossSections.h
 *
 *  Created on: 2013.03.12
 *      Author: koehne
 */

#ifndef CrossSections_H
#define CrossSections_H


#include "Integral.h"



/*! \class CrossSections CrossSections.h "CrossSections.h"
    \brief This is a pure virtual class
 */



class CrossSections
{


protected:

    // bounds of parameterizations
    double elow_;
    double nlow_;
    double ebig_;


    // Interpolation flags
    bool doContinuousInterpolation_;
    bool doStochasticInterpolation_;

    //CrossSection multiplier
    double multiplier_;

public:

    //Constructor
    CrossSections();

//----------------------------------------------------------------------------//

    // Memberfunctions

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

    virtual double FunctionToContinuousIntegral(double variable) = 0;

//----------------------------------------------------------------------------//

    virtual double FunctionToStochasticalIntegral(double variable) = 0;


//----------------------------------------------------------------------------//

    void SetParameterizationLimits(double elow=0.,
                                   double nlow=ME,
                                   double ebig=BIGENERGY);
    
//----------------------------------------------------------------------------//

    // Setter

    void SetMultiplier(double multiplier=1.);


//----------------------------------------------------------------------------//

    // Getter

    double GetElow() const
    {
        return elow_;
    }

//----------------------------------------------------------------------------//

    double GetNlow() const
    {
        return nlow_;
    }

//----------------------------------------------------------------------------//

    double GetEbig() const
    {
        return ebig_;
    }

//----------------------------------------------------------------------------//

    double GetMultiplier() const
    {
        return multiplier_;
    }

//----------------------------------------------------------------------------//

    // destructor

    virtual ~CrossSections(){}

};



#endif //CrossSections_H
