#include "PROPOSAL/CrossSections.h"
#include "PROPOSAL/Constants.h"


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


CrossSections::CrossSections( )
    :vMax_                  ( 0 )
    ,vUp_                   ( 0 )
    ,vMin_                  ( 0 )
    ,ebig_                  ( BIGENERGY )
    ,rnd_                   ( 0 )
    ,do_dedx_Interpolation_ ( false )
    ,do_dndx_Interpolation_ ( false )
    ,multiplier_            ( 1. )
    ,parametrization_       ( 1 )
    ,lpm_effect_enabled_    ( false )
    ,init_lpm_effect_       ( true )
    ,order_of_interpolation_( 5 )
    ,sum_of_rates_          ( 0 )
{
    particle_       = new Particle();
    medium_         = new Medium();
    cut_settings_   = new EnergyCutSettings(-1,-1);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


CrossSections::CrossSections(Particle* particle,
                             Medium* medium,
                             EnergyCutSettings* cut_settings)
    :vMax_                  ( 0 )
    ,vUp_                   ( 0 )
    ,vMin_                  ( 0 )
    ,ebig_                  ( BIGENERGY )
    ,rnd_                   ( 0 )
    ,do_dedx_Interpolation_ ( false )
    ,do_dndx_Interpolation_ ( false )
    ,multiplier_            ( 1. )
    ,parametrization_       ( 1 )
    ,lpm_effect_enabled_    ( false )
    ,init_lpm_effect_       ( true )
    ,order_of_interpolation_( 5 )
    ,sum_of_rates_          ( 0 )
{
    particle_       = particle;
    medium_         = medium;
    cut_settings_   = cut_settings;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


CrossSections::CrossSections(const CrossSections& crossSections)
    :name_                     ( crossSections.name_ )
    ,vMax_                     ( crossSections.vMax_ )
    ,vUp_                      ( crossSections.vUp_ )
    ,vMin_                     ( crossSections.vMin_ )
    ,ebig_                     ( crossSections.ebig_ )
    ,rnd_                      ( crossSections.rnd_ )
    ,do_dedx_Interpolation_    ( crossSections.do_dedx_Interpolation_ )
    ,do_dndx_Interpolation_    ( crossSections.do_dndx_Interpolation_ )
    ,multiplier_               ( crossSections.multiplier_ )
    ,parametrization_          ( crossSections.parametrization_ )
    ,lpm_effect_enabled_       ( crossSections.lpm_effect_enabled_ )
    ,init_lpm_effect_          ( crossSections.init_lpm_effect_ )
    ,order_of_interpolation_   ( crossSections.order_of_interpolation_ )
    ,sum_of_rates_             ( crossSections.sum_of_rates_ )
{
    particle_                 = new Particle( *crossSections.particle_ );
    medium_                   = new Medium( *crossSections.medium_ );
    cut_settings_             = new EnergyCutSettings( *crossSections.cut_settings_ );
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool CrossSections::operator==(const CrossSections &crossSections) const
{
    if( vMax_                     != crossSections.vMax_ )                  return false;
    if( vUp_                      != crossSections.vUp_ )                   return false;
    if( vMin_                     != crossSections.vMin_ )                  return false;
    if( ebig_                     != crossSections.ebig_ )                  return false;
    if( rnd_                      != crossSections.rnd_ )                   return false;
    if( do_dedx_Interpolation_    != crossSections.do_dedx_Interpolation_ ) return false;
    if( do_dndx_Interpolation_    != crossSections.do_dndx_Interpolation_ ) return false;
    if( multiplier_               != crossSections.multiplier_ )            return false;
    if( parametrization_          != crossSections.parametrization_ )       return false;
    if( lpm_effect_enabled_       != crossSections.lpm_effect_enabled_ )    return false;
    if( init_lpm_effect_          != crossSections.init_lpm_effect_ )       return false;
    if( order_of_interpolation_   != crossSections.order_of_interpolation_ )return false;
    if( *cut_settings_            != *crossSections.cut_settings_ )         return false;
    if( *particle_                != *crossSections.particle_ )             return false;
    if( *medium_                  != *crossSections.medium_ )               return false;
    if( sum_of_rates_             != crossSections.sum_of_rates_ )          return false;
    if( name_.compare(crossSections.name_) != 0 )                           return false;



    //else
    return true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool CrossSections::operator!=(const CrossSections &crossSections) const
{
    return !(*this == crossSections);
}

std::ostream& operator<<(std::ostream& os, CrossSections const &crossSections)
{
    os<<"\tParametrization:\t" << crossSections.parametrization_ << std::endl;
    os<<"\tEbig:\t\t"<< crossSections.ebig_ << std::endl;
    os<<"\tMultiplier:\t\t"<<crossSections.multiplier_<<std::endl;
    os<<std::endl;
    os<<"\tParticle:\t\t"<<crossSections.particle_<<std::endl;
    if(crossSections.particle_!=NULL)
    {
        os<<"\t\tname:\t\t\t"<<crossSections.particle_->GetName()<<std::endl;
        os<<"\t\tenergy:\t\t\t"<<crossSections.particle_->GetEnergy()<<std::endl;
        os<<"\t\tdistance:\t\t\t"<<crossSections.particle_->GetPropagatedDistance()<<std::endl;
    }
    os<<std::endl;
    os<<"\tMedium:\t\t"<<crossSections.medium_<<std::endl;
    if(crossSections.medium_!=NULL)
    {
        os<<"\t\tname:\t\t\t"<<crossSections.medium_->GetName()<<std::endl;
        os<<"\t\trho:\t\t\t"<<crossSections.medium_->GetRho()<<std::endl;
    }
    os<<std::endl;
    os<<"\tCutSettings:\t"<<crossSections.cut_settings_<<std::endl;
    if(crossSections.cut_settings_!=NULL)
    {
        os<<"\t\tecut:\t\t\t"<<crossSections.cut_settings_->GetEcut()<<std::endl;
        os<<"\t\tvcut:\t\t\t"<<crossSections.cut_settings_->GetVcut()<<std::endl;
    }
    os<<std::endl;
    os<<"\tTemp. variables: " << std::endl;
    os<<"\t\tvMax:\t\t\t"<<crossSections.vMax_<<std::endl;
    os<<"\t\tvUp:\t\t\t"<<crossSections.vUp_<<std::endl;
    os<<"\t\tvMin:\t\t\t"<<crossSections.vMin_<<std::endl;
    os<<"\t\trnd:\t\t\t"<<crossSections.rnd_<<std::endl;
    os<<"\t\torder_of_interpolation:\t"<<crossSections.order_of_interpolation_<<std::endl;
    os<<"\t\tsum_of_rates:\t\t"<<crossSections.sum_of_rates_<<std::endl;
    os<<"\t\tinit_lpm_effect:\t\t"<<crossSections.init_lpm_effect_<<std::endl;
    os<<"\t\tlpm_effect_enabled:\t\t"<<crossSections.lpm_effect_enabled_<<std::endl;
    return os;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void CrossSections::swap(CrossSections &crossSections)
{
    using std::swap;

    swap( vMax_                     , crossSections.vMax_ );
    swap( vUp_                      , crossSections.vUp_ );
    swap( vMin_                     , crossSections.vMin_ );
    swap( ebig_                     , crossSections.ebig_ );
    swap( rnd_                      , crossSections.rnd_ );
    swap( do_dedx_Interpolation_    , crossSections.do_dedx_Interpolation_ );
    swap( do_dndx_Interpolation_    , crossSections.do_dndx_Interpolation_ );
    swap( multiplier_               , crossSections.multiplier_ );
    swap( parametrization_          , crossSections.parametrization_ );
    swap( lpm_effect_enabled_       , crossSections.lpm_effect_enabled_ );
    swap( init_lpm_effect_          , crossSections.init_lpm_effect_ );
    swap( order_of_interpolation_   , crossSections.order_of_interpolation_ );
    swap( sum_of_rates_             , crossSections.sum_of_rates_ );

    cut_settings_->swap(*crossSections.cut_settings_);
    particle_->swap(*crossSections.particle_);
    medium_->swap(*crossSections.medium_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void CrossSections::SetParametrizationLimit(double ebig){
    ebig_ = ebig;
}

void CrossSections::SetMultiplier(double multiplier){
    multiplier_ = multiplier;
}

void CrossSections::SetVMin(double vMin){
    vMin_ = vMin;
}

void CrossSections::SetVMax(double vMax){
    vMax_ = vMax;
}

void CrossSections::SetVUp(double vUp){
    vUp_ = vUp;
}

void CrossSections::SetParticle(Particle *particle){
    particle_ = particle;
}

void CrossSections::SetMedium(Medium *medium){
    medium_ = medium;
}

void CrossSections::SetEnergyCutSettings(EnergyCutSettings *cuts){
    cut_settings_ = cuts;
}

void CrossSections::EnableLpmEffect(bool lpm_effect_enabled){
    lpm_effect_enabled_ = lpm_effect_enabled;
}


