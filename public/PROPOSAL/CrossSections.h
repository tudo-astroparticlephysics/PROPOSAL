/*
 * CrossSections.h
 *
 *  Created on: 2013.03.12
 *      Author: Jan-Hendrik Köhne
 */

#ifndef CrossSections_H
#define CrossSections_H


#include "PROPOSAL/MathModel.h"
#include "PROPOSAL/Constants.h"




/*! \class CrossSections CrossSections.h "CrossSections.h"
    \brief This is a pure virtual class
 */



class CrossSections : public MathModel
{


protected:

    //bounds of integration
    double vMax_;   //!< upper bound of integration
    double vUp_;    //!< lower bound of integration
    double vMin_;   //!< lowest physical possible bound of integration

    double ebig_;   //!< upper bound of parameterizations

    // Interpolation flags
    bool doContinuousInterpolation_;
    bool doStochasticInterpolation_;

    //CrossSection multiplier
    double multiplier_;

    int parametrization_;

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

    void SetParametrizationLimit(double ebig=BIGENERGY);
    
//----------------------------------------------------------------------------//

    // Setter

    void SetMultiplier(double multiplier=1.);

//----------------------------------------------------------------------------//

    void SetVMin(double vMin=0);

//----------------------------------------------------------------------------//

    void SetVMax(double vMax=0);

//----------------------------------------------------------------------------//

    void SetVUp(double vUp=0);

//----------------------------------------------------------------------------//

    void SetParametrization(int parametrization=1);

//----------------------------------------------------------------------------//
    // Getter

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

    double GetVMax() const
    {
        return vMax_;
    }

//----------------------------------------------------------------------------//

    double GetVMin() const
    {
        return vMin_;
    }

//----------------------------------------------------------------------------//

    double GetVUp() const
    {
        return vUp_;
    }

//----------------------------------------------------------------------------//

    double GetParametrization() const
    {
        return parametrization_;
    }

//----------------------------------------------------------------------------//

    // destructor

    virtual ~CrossSections(){}

};



#endif //CrossSections_H
