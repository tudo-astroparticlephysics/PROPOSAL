/*
 * Energy2Loss.h
 *
 *  Created on: 24.06.2010
 *      Author: koehne
 */

#ifndef ENERGY2LOSS_H_
#define ENERGY2LOSS_H_

#include "PROPOSAL/PhysicsModel.h"
#include "PROPOSAL/CrossSections.h"

/**
 * class contains functions for evaluation of
 * the spread of the continuous energy losses
 */

class Energy2LossX;
class Energy2LossE;
class PROPOSALParticle;
class Medium;




class Energy2Loss :public  PhysicsModel
{

public:

    Energy2LossX *e2lx;
    Energy2LossE *e2le;
    PROPOSALParticle *particle_;
    Medium *medium_;
    CrossSections *cros;
    Energy2Loss *e2loss;


//----------------------------------------------------------------------------//
    //Constructors

    /**
    * initializes an integral class separate from that in Propagate
    */

    /**
    * creates internal reference to superclass, to be called from subclasses
    */

    Energy2Loss(Energy2Loss *e2loss_);
//----------------------------------------------------------------------------//

    /**
    * initializes Energy2Loss classes, creates internal reference to CrossSections
    */

    Energy2Loss(CrossSections *cros);

//----------------------------------------------------------------------------//

    Energy2Loss(){}

 //----------------------------------------------------------------------------//

    // Destructor

    ~Energy2Loss(){}


};


#endif /* ENERGY2LOSS_H_ */
