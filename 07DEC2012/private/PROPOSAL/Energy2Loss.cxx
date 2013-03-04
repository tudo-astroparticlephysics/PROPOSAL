/*
 * Energy2Loss.cxx
 *
 *  Created on: 24.06.2010
 *      Author: koehne
 */

#include "PROPOSAL/Energy2Loss.h"
#include "PROPOSAL/Energy2LossX.h"
#include "PROPOSAL/Energy2LossE.h"


Energy2Loss::Energy2Loss(Energy2Loss *e2loss_)
{
    particle_       =   e2loss_->particle_;
    medium_         =   e2loss_->medium_;
    this->e2loss    =   e2loss_;
    this->cros      =   e2loss_->cros;
}

//----------------------------------------------------------------------------//

Energy2Loss::Energy2Loss(CrossSections *cros)
{
    this->particle_ =   cros->particle_;
    this->medium_   =   cros->medium_;
    this->cros      =   cros;
    e2lx            =   new Energy2LossX(this);
    e2le            =   new Energy2LossE(this);
}


