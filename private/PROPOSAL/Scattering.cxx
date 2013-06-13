/*! \file   Scattering.cxx
*   \brief  Source file for the scattering routines.
*
*   For more details see the class documentation.
*
*   \date   2013.06.13
*   \author Tomasz Fuchs
*/


#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include "PROPOSAL/Particle.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Scattering.h"

using namespace std;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Scattering::Scattering( )
{

}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Scattering& Scattering::operator=(const Scattering &scattering){
    if (this != &scattering)
    {
      Scattering tmp(scattering);
      swap(tmp);
    }
    return *this;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Scattering::operator==(const Scattering &scattering) const
{
    /*
    if(  propagated_distance_   != particle.propagated_distance_)   return false;
    if(  x_                     != particle.x_)                     return false;
    if(  y_                     != particle.y_)                     return false;
    if(  z_                     != particle.z_)                     return false;
    if(  t_                     != particle.t_)                     return false;
    if(  theta_                 != particle.theta_)                 return false;
    if(  phi_                   != particle.phi_)                   return false;
    if(  costh_                 != particle.costh_)                 return false;
    if(  sinth_                 != particle.sinth_)                 return false;
    if(  cosph_                 != particle.cosph_)                 return false;
    if(  sinph_                 != particle.sinph_)                 return false;
    if(  momentum_              != particle.momentum_)              return false;
    if(  square_momentum_       != particle.square_momentum_)       return false;
    if(  energy_                != particle.energy_)                return false;
    if(  mass_                  != particle.mass_)                  return false;
    if(  lifetime_              != particle.lifetime_)              return false;
    if(  charge_                != particle.charge_)                return false;
    if(  low_                   != particle.low_)                   return false;
    if(  type_                  != particle.type_)                  return false;
    if(  parent_particle_id_    != particle.parent_particle_id_)    return false;
    if(  particle_id_           != particle.particle_id_)           return false;
    if(  xi_                    != particle.xi_)                    return false;
    if(  yi_                    != particle.yi_)                    return false;
    if(  zi_                    != particle.zi_)                    return false;
    if(  ti_                    != particle.ti_)                    return false;
    if(  ei_                    != particle.ei_)                    return false;
    if(  xf_                    != particle.xf_)                    return false;
    if(  yf_                    != particle.yf_)                    return false;
    if(  zf_                    != particle.zf_)                    return false;
    if(  tf_                    != particle.tf_)                    return false;
    if(  ef_                    != particle.ef_)                    return false;
    if(  xc_                    != particle.xc_)                    return false;
    if(  yc_                    != particle.yc_)                    return false;
    if(  zc_                    != particle.zc_)                    return false;
    if(  tc_                    != particle.tc_)                    return false;
    if(  ec_                    != particle.ec_)                    return false;
    if(  elost_                 != particle.elost_)                 return false;

    if(  name_.compare(particle.name_) != 0)        return false;
*/
    //else
    return true;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Scattering::operator!=(const Scattering &scattering) const {
  return !(*this == scattering);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Scattering::swap(Scattering &scattering)
{
    using std::swap;
/*
    swap( propagated_distance_   , particle.propagated_distance_);
    swap( x_                     , particle.x_);
    swap( y_                     , particle.y_);
    swap( z_                     , particle.z_);
    swap( t_                     , particle.t_);
    swap( theta_                 , particle.theta_);
    swap( phi_                   , particle.phi_);
    swap( costh_                 , particle.costh_);
    swap( sinth_                 , particle.sinth_);
    swap( cosph_                 , particle.cosph_);
    swap( sinph_                 , particle.sinph_);
    swap( momentum_              , particle.momentum_);
    swap( square_momentum_       , particle.square_momentum_);
    swap( energy_                , particle.energy_);
    swap( mass_                  , particle.mass_);
    swap( lifetime_              , particle.lifetime_);
    swap( charge_                , particle.charge_);
    swap( low_                   , particle.low_);
    swap( type_                  , particle.type_);
    swap( parent_particle_id_    , particle.parent_particle_id_);
    swap( particle_id_           , particle.particle_id_);
    swap( xi_                    , particle.xi_);
    swap( yi_                    , particle.yi_);
    swap( zi_                    , particle.zi_);
    swap( ti_                    , particle.ti_);
    swap( ei_                    , particle.ei_);
    swap( xf_                    , particle.xf_);
    swap( yf_                    , particle.yf_);
    swap( zf_                    , particle.zf_);
    swap( tf_                    , particle.tf_);
    swap( ef_                    , particle.ef_);
    swap( xc_                    , particle.xc_);
    swap( yc_                    , particle.yc_);
    swap( zc_                    , particle.zc_);
    swap( tc_                    , particle.tc_);
    swap( ec_                    , particle.ec_);
    swap( elost_                 , particle.elost_);

    name_.swap(particle.name_);
*/

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------private member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//




//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

