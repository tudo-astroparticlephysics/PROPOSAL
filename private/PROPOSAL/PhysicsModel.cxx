/*! \file   PhysicsModel.cxx
*   \brief  Source file for the physics model routines.
*
*   For more details see the class documentation.
*
*   \date   23.06.2010
*   \author Jan-Hendrik Koehne
*/


#include "PROPOSAL/PhysicsModel.h"

double PhysicsModel::ebig_;
double PhysicsModel::elow_;
double PhysicsModel::nlow_;
using namespace std;

// Some defaultvariables


PhysicsModel::PhysicsModel()
{

    elow_=0;                //!< lower bound of parameterizations
    nlow_=ME;               //!< maximal number of parametrization points
    ebig_=BIGENERGY;        //!< upper bound of parameterizations



}

//----------------------------------------------------------------------------//

PhysicsModel::~PhysicsModel() {}


