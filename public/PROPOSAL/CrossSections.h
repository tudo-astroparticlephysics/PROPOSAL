/*
 * CrossSections.h
 *
 *  Created on: 2013.03.12
 *      Author: koehne
 */

#ifndef CrossSections_H
#define CrossSections_H

/*! \class CrossSections CrossSections.h "CrossSections.h"
    \brief This is a pure virtual class
 */


class CrossSections
{


protected:

    CrossSections();

public:

    // Memberfunctions

//----------------------------------------------------------------------------//

    virtual void SetIntegralLimits() const = 0;

//----------------------------------------------------------------------------//

    virtual double DEdx() const = 0;

//----------------------------------------------------------------------------//

    virtual double DNdx() const = 0;


//----------------------------------------------------------------------------//

    virtual double DNdx(double rnd) const = 0;

//----------------------------------------------------------------------------//

    virtual double CalculateStochasticLoss() const = 0;

//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//

    // Setter


    
//----------------------------------------------------------------------------//

    // destructors

    ///@brief Crush this CrossSections.
    virtual ~CrossSections();


};



#endif //CrossSections_H
