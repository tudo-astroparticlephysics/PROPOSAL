
// #include <algorithm>

#include <boost/bind.hpp>

#include "PROPOSAL/crossection/Photonuclear.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"

using namespace std;
using namespace PROPOSAL;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::CalculatedEdx(const PROPOSALParticle& particle)
{
    if(multiplier_<=0)
    {
        return 0;
    }

    if(do_dedx_Interpolation_)
    {
        return max(dedx_interpolant_->Interpolate(particle.GetEnergy()), 0.0);
    }

    double sum = 0;

    for(int i=0; i < medium_->GetNumComponents(); i++)
    {
        CrossSections::IntegralLimits limits = SetIntegralLimits(particle, i);
        sum +=  integral_for_dEdx_->Integrate(limits.vMin, limits.vUp, boost::bind(&Photonuclear::FunctionToDEdxIntegral, this, boost::cref(particle), _1),4);
    }

    return multiplier_*particle.GetEnergy()*sum;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::CalculatedNdx(const PROPOSALParticle& particle)
{
    if(multiplier_<=0)
    {
        return 0;
    }

    sum_of_rates_ = 0;

    for(int i=0; i<medium_->GetNumComponents(); i++)
    {
        if(do_dndx_Interpolation_)
        {
            prob_for_component_.at(i) = max(dndx_interpolant_1d_.at(i)->Interpolate(particle.GetEnergy()), 0.);
        }
        else
        {
            CrossSections::IntegralLimits limits = SetIntegralLimits(particle, i);
            prob_for_component_.at(i) = dndx_integral_.at(i)->Integrate(limits.vUp, limits.vMax, boost::bind(&Photonuclear::FunctionToDNdxIntegral, this, boost::cref(particle), _1),4);
        }
        sum_of_rates_ += prob_for_component_.at(i);
    }
    return sum_of_rates_;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::CalculatedNdx(const PROPOSALParticle& particle, double rnd)
{
    if(multiplier_<=0)
    {
        return 0.;
    }

    // The random number will be stored to be able
    // to check if dNdx is already calculated for this random number.
    // This avoids a second calculation in CalculateStochasticLoss

    rnd_ = rnd;
    sum_of_rates_ = 0;

    for(int i=0; i<medium_->GetNumComponents(); i++)
    {
        if(do_dndx_Interpolation_)
        {
            prob_for_component_.at(i) = max(dndx_interpolant_1d_.at(i)->Interpolate(particle.GetEnergy()), 0.);
        }
        else
        {
            CrossSections::IntegralLimits limits = SetIntegralLimits(particle, i);
            prob_for_component_.at(i) = dndx_integral_.at(i)->IntegrateWithRandomRatio(limits.vUp, limits.vMax, boost::bind(&Photonuclear::FunctionToDNdxIntegral, this, boost::cref(particle), _1),4,rnd);
        }
        sum_of_rates_ += prob_for_component_.at(i);
    }

    return sum_of_rates_;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::CalculateStochasticLoss(const PROPOSALParticle& particle, double rnd1, double rnd2)
{
    if(rnd1 != rnd_ )
    {
        CalculatedNdx(particle, rnd1);
    }

    return CalculateStochasticLoss(particle, rnd2);

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-----------------------Enable and disable interpolation---------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Photonuclear::EnableDNdxInterpolation(const PROPOSALParticle& particle, std::string path, bool raw)
{
    if(do_dndx_Interpolation_)return;

    EnablePhotoInterpolation(particle, path,raw);

    bool storing_failed =   false;
    bool reading_worked =   true;

    // charged anti leptons have the same cross sections like charged leptons
    // so they use the same interpolation tables

    if(!path.empty())
    {
        stringstream filename;
        filename<<path<<"/Photo_dNdx"
                <<"_particle_"<<particle.GetName()
                <<"_mass_"<<particle.GetMass()
                <<"_charge_"<<particle.GetCharge()
                <<"_lifetime_"<<particle.GetLifetime()
                <<"_para_"<<parametrization_
                <<"_med_"<<medium_->GetName()
                <<"_"<<medium_->GetMassDensity()
                <<"_ecut_"<<cut_settings_->GetEcut()
                <<"_vcut_"<<cut_settings_->GetVcut()
                <<"_multiplier_"<<multiplier_;

        if(!raw)
            filename<<".txt";

        dndx_interpolant_2d_.resize(medium_->GetNumComponents());
        dndx_interpolant_1d_.resize(medium_->GetNumComponents());

        if( FileExist(filename.str()) )
        {
            log_debug("Photonuclear parametrisation tables (dNdx) will be read from file:\t%s",filename.str().c_str());
            ifstream input;

            if(raw)
            {
                input.open(filename.str().c_str(), ios::binary);
            }
            else
            {
                input.open(filename.str().c_str());
            }

            for(int i=0; i<(medium_->GetNumComponents()); i++)
            {
                component_ = i;
                dndx_interpolant_2d_.at(i) = new Interpolant();
                dndx_interpolant_1d_.at(i) = new Interpolant();
                reading_worked = dndx_interpolant_2d_.at(i)->Load(input, raw);
                reading_worked = dndx_interpolant_1d_.at(i)->Load(input, raw);

            }
            input.close();
        }
        if(!FileExist(filename.str()) || !reading_worked )
        {
            if(!reading_worked)
            {
                log_info("File %s is corrupted! Write it again!",filename.str().c_str());
            }

            log_info("Photonuclear parametrisation tables (dNdx) will be saved to file:\t%s",filename.str().c_str());

            ofstream output;

            if(raw)
            {
                output.open(filename.str().c_str(), ios::binary);
            }
            else
            {
                output.open(filename.str().c_str());
            }

            if(output.good())
            {
                output.precision(16);

                for(int i=0; i<(medium_->GetNumComponents()); i++)
                {
                    component_ = i;

                    dndx_interpolant_2d_.at(i) = new Interpolant(NUM1
                                                                , particle.GetLow()
                                                                , BIGENERGY
                                                                , NUM1
                                                                , 0
                                                                , 1
                                                                , boost::bind(&Photonuclear::FunctionToBuildDNdxInterpolant2D, this, boost::cref(particle), _1 , _2)
                                                                , order_of_interpolation_
                                                                , false
                                                                , false
                                                                , true
                                                                , order_of_interpolation_
                                                                , false
                                                                , false
                                                                , false
                                                                , order_of_interpolation_
                                                                , true
                                                                , false
                                                                , false
                                                                );
                    dndx_interpolant_1d_.at(i) = new Interpolant(NUM1
                                                                , particle.GetLow()
                                                                , BIGENERGY
                                                                , boost::bind(&Photonuclear::FunctionToBuildDNdxInterpolant1D, this, _1)
                                                                , order_of_interpolation_
                                                                , false
                                                                , false
                                                                , true
                                                                , order_of_interpolation_
                                                                , true
                                                                , false
                                                                , false
                                                                );

                    dndx_interpolant_2d_.at(i)->Save(output, raw);
                    dndx_interpolant_1d_.at(i)->Save(output, raw);

                }
            }
            else
            {
                storing_failed  =   true;
                log_warn("Can not open file %s for writing! Table will not be stored!",filename.str().c_str());
            }

            output.close();
        }
    }
    if(path.empty() || storing_failed)
    {
        EnablePhotoInterpolation(particle, path);

        dndx_interpolant_1d_.resize(medium_->GetNumComponents());
        dndx_interpolant_2d_.resize(medium_->GetNumComponents());
        for(int i=0; i<(medium_->GetNumComponents()); i++)
        {
            component_ = i;
            dndx_interpolant_2d_.at(i) =    new Interpolant(NUM1
                                                        , particle.GetLow()
                                                        , BIGENERGY
                                                        , NUM1
                                                        , 0
                                                        , 1
                                                        , boost::bind(&Photonuclear::FunctionToBuildDNdxInterpolant2D, this, boost::cref(particle), _1 , _2)
                                                        , order_of_interpolation_
                                                        , false
                                                        , false
                                                        , true
                                                        , order_of_interpolation_
                                                        , false
                                                        , false
                                                        , false
                                                        , order_of_interpolation_
                                                        , true
                                                        , false
                                                        , false
                                                        );
            dndx_interpolant_1d_.at(i) =    new Interpolant(NUM1
                                                        , particle.GetLow()
                                                        , BIGENERGY
                                                        , boost::bind(&Photonuclear::FunctionToBuildDNdxInterpolant1D, this, _1)
                                                        , order_of_interpolation_
                                                        , false
                                                        , false
                                                        , true
                                                        , order_of_interpolation_
                                                        , true
                                                        , false
                                                        , false
                                                        );
        }
    }

    do_dndx_Interpolation_=true;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Photonuclear::EnableDEdxInterpolation(const PROPOSALParticle& particle, std::string path, bool raw)
{
    if(do_dedx_Interpolation_)return;

    EnablePhotoInterpolation(particle, path,raw);

    bool reading_worked =   true;
    bool storing_failed =   false;

    // charged anti leptons have the same cross sections like charged leptons
    // so they use the same interpolation tables

    if(!path.empty())
    {
        stringstream filename;
        filename<<path<<"/Photo_dEdx"
                <<"_particle_"<<particle.GetName()
                <<"_mass_"<<particle.GetMass()
                <<"_charge_"<<particle.GetCharge()
                <<"_lifetime_"<<particle.GetLifetime()
                <<"_para_"<<parametrization_
                <<"_med_"<<medium_->GetName()
                <<"_"<<medium_->GetMassDensity()
                <<"_ecut_"<<cut_settings_->GetEcut()
                <<"_vcut_"<<cut_settings_->GetVcut()
                <<"_multiplier_"<<multiplier_;

        if(!raw)
            filename<<".txt";

        if( FileExist(filename.str()) )
        {
            log_debug("Photonuclear parametrisation tables (dEdx) will be read from file:\t%s",filename.str().c_str());
            ifstream input;

            if(raw)
            {
                input.open(filename.str().c_str(), ios::binary);
            }
            else
            {
                input.open(filename.str().c_str());
            }

            dedx_interpolant_ = new Interpolant();
            reading_worked = dedx_interpolant_->Load(input, raw);

            input.close();
        }
        if(!FileExist(filename.str()) || !reading_worked )
        {
            if(!reading_worked)
            {
                log_info("File %s is corrupted! Write it again!",filename.str().c_str());
            }

            log_info("Photonuclear parametrisation tables (dEdx) will be saved to file:\t%s",filename.str().c_str());

            ofstream output;

            if(raw)
            {
                output.open(filename.str().c_str(), ios::binary);
            }
            else
            {
                output.open(filename.str().c_str());
            }

            if(output.good())
            {
                output.precision(16);

                dedx_interpolant_ = new Interpolant(NUM1, particle.GetLow(), BIGENERGY, boost::bind(&Photonuclear::FunctionToBuildDEdxInterpolant, this, boost::cref(particle), _1),
                                                    order_of_interpolation_, true, false, true, order_of_interpolation_, false, false, false);
                dedx_interpolant_->Save(output, raw);
            }
            else
            {
                storing_failed  =   true;
                log_warn("Can not open file %s for writing! Table will not be stored!",filename.str().c_str());
            }

            output.close();
        }
    }
    if(path.empty() || storing_failed)
    {
        dedx_interpolant_ = new Interpolant(NUM1, particle.GetLow(), BIGENERGY, boost::bind(&Photonuclear::FunctionToBuildDEdxInterpolant, this, boost::cref(particle), _1),
                                            order_of_interpolation_, true, false, true, order_of_interpolation_, false, false, false);

    }

    do_dedx_Interpolation_=true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Photonuclear::EnablePhotoInterpolation(const PROPOSALParticle& particle, std::string path, bool raw)
{
    if(do_photo_interpolation_)return;

    bool storing_failed =   false;
    bool reading_worked =   true;

    // charged anti leptons have the same cross sections like charged leptons
    // so they use the same interpolation tables

    if(!path.empty())
    {
        stringstream filename;
        filename<<path<<"/Photo"
                <<"_particle_"<< particle.GetName()
                <<"_mass_"<<particle.GetMass()
                <<"_charge_"<<particle.GetCharge()
                <<"_lifetime_"<<particle.GetLifetime()
                <<"_para_"<<parametrization_
                <<"_med_"<<medium_->GetName()
                <<"_"<<medium_->GetMassDensity()
                <<"_ecut_"<<cut_settings_->GetEcut()
                <<"_vcut_"<<cut_settings_->GetVcut()
                <<"_multiplier_"<<multiplier_;

        if(!raw)
            filename<<".txt";

        photo_interpolant_.resize( medium_->GetNumComponents() );

        if( FileExist(filename.str()) )
        {
            log_debug("Photonuclear parametrisation tables will be read from file:\t%s",filename.str().c_str());
            ifstream input;

            if(raw)
            {
                input.open(filename.str().c_str(), ios::binary);
            }
            else
            {
                input.open(filename.str().c_str());
            }

            for(int i=0; i<(medium_->GetNumComponents()); i++)
            {
                component_ = i;
                photo_interpolant_.at(i) = new Interpolant();
                reading_worked = photo_interpolant_.at(i)->Load(input, raw);

            }
            input.close();
        }
        if(!FileExist(filename.str()) || !reading_worked )
        {
            if(!reading_worked)
            {
                log_info("File %s is corrupted! Write it again!",filename.str().c_str());
            }

            log_info("Photonuclear parametrisation tables will be saved to file:\t%s",filename.str().c_str());

            ofstream output;

            if(raw)
            {
                output.open(filename.str().c_str(), ios::binary);
            }
            else
            {
                output.open(filename.str().c_str());
            }

            if(output.good())
            {
                output.precision(16);

                for(int i=0; i<(medium_->GetNumComponents()); i++)
                {
                    component_ = i;

                    photo_interpolant_.at(i)  = new Interpolant(NUM1
                                                            , particle.GetLow()
                                                            , BIGENERGY
                                                            , NUM1
                                                            , 0.
                                                            , 1.
                                                            , boost::bind(&Photonuclear::FunctionToBuildPhotoInterpolant, this, boost::cref(particle), _1, _2)
                                                            , order_of_interpolation_
                                                            , false
                                                            , false
                                                            , true
                                                            , order_of_interpolation_
                                                            , false
                                                            , false
                                                            , false
                                                            , order_of_interpolation_
                                                            , false
                                                            , false
                                                            , false
                                                            );

                    photo_interpolant_.at(i)->Save(output, raw);

                }
            }
            else
            {
                storing_failed  =   true;
                log_warn("Can not open file %s for writing! Table will not be stored!",filename.str().c_str());
            }

            output.close();
        }
    }
    if(path.empty() || storing_failed)
    {
        photo_interpolant_.resize( medium_->GetNumComponents() );

        for(int i=0; i<medium_->GetNumComponents(); i++)
        {
            component_ = i;
            photo_interpolant_.at(i)  = new Interpolant(NUM1
                                                    , particle.GetLow()
                                                    , BIGENERGY
                                                    , NUM1
                                                    , 0.
                                                    , 1.
                                                    , boost::bind(&Photonuclear::FunctionToBuildPhotoInterpolant, this, boost::cref(particle), _1, _2)
                                                    , order_of_interpolation_
                                                    , false
                                                    , false
                                                    , true
                                                    , order_of_interpolation_
                                                    , false
                                                    , false
                                                    , false
                                                    , order_of_interpolation_
                                                    , false
                                                    , false
                                                    , false
                                                    );

        }
    }

    do_photo_interpolation_ = true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Photonuclear::DisableDNdxInterpolation()
{
    for(unsigned int i = 0 ; i < dndx_interpolant_1d_.size() ; i++ )
    {
        delete dndx_interpolant_1d_.at(i);
    }

    for(unsigned int i = 0 ; i < dndx_interpolant_2d_.size() ; i++ )
    {
        delete dndx_interpolant_2d_.at(i);
    }

    dndx_interpolant_1d_.clear();
    dndx_interpolant_2d_.clear();

    do_dndx_Interpolation_  =   false;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Photonuclear::DisableDEdxInterpolation()
{
    delete dedx_interpolant_;
    dedx_interpolant_   = NULL;
    do_dedx_Interpolation_  =   false;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Photonuclear::DisablePhotoInterpolation()
{
    for(unsigned int i=0; i < photo_interpolant_.size(); i++)
    {
        delete photo_interpolant_.at(i);
    }
    photo_interpolant_.clear();
    do_photo_interpolation_  =   false;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------Set and validate options--------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// void Photonuclear::ValidateOptions()
// {
//     if(order_of_interpolation_ < 2)
//     {
//         order_of_interpolation_ = 5;
//         cerr<<"Photonuclear: Order of Interpolation is not a vaild number. Must be > 2.\t"<<"Set to 5"<<endl;
//     }
//     if(order_of_interpolation_ > 6)
//     {
//         cerr<<"Photonuclear: Order of Interpolation is set to "<<order_of_interpolation_
//             <<".\t Note a order of interpolation > 6 will slow down the program"<<endl;
//     }
//     switch (parametrization_)
//     {
//         case ParametrizationType::PhotoKokoulinShadowBezrukovSoft:
//         case ParametrizationType::PhotoKokoulinShadowBezrukovHard:
//         case ParametrizationType::PhotoRhodeShadowBezrukovSoft:
//         case ParametrizationType::PhotoRhodeShadowBezrukovHard:
//         case ParametrizationType::PhotoBezrukovBugaevShadowBezrukovSoft:
//         case ParametrizationType::PhotoBezrukovBugaevShadowBezrukovHard:
//         case ParametrizationType::PhotoZeusShadowBezrukovSoft:
//         case ParametrizationType::PhotoZeusShadowBezrukovHard:
//         case ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowDutta:
//         case ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowButkevich:
//         case ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowDutta:
//         case ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich:
//         case ParametrizationType::PhotoButkevichMikhailovShadowDutta:
//         case ParametrizationType::PhotoButkevichMikhailovShadowButkevich:
//             break;
//         default:
//             SetParametrization(ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich);
//             cerr<<"Photonuclear: Parametrization is not a vaild number. \t Set parametrization to ALLM97 with Butkevich Shadowing"<<endl;
//             break;
//     }
// }
//

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Photonuclear::Photonuclear()
    :component_             ( 0 )
    ,init_measured_         ( true )
    ,init_hardbb_           ( true )
    ,hmax_                  ( 8 )
    ,v_                     ( 0 )
    ,do_photo_interpolation_( false )
    ,hard_component_        ( false )
    ,shadow_                ( ShadowingType::ButkevichMikhailov )
    ,dndx_integral_         ( )
    ,interpolant_hardBB_    ( )
    ,dndx_interpolant_1d_   ( )
    ,dndx_interpolant_2d_   ( )
    ,photo_interpolant_     ( )
    ,prob_for_component_    ( )
    ,woodSaxonPotential_    ( )

{
    name_   =   "Photonuclear";
    type_   =   ParticleType::NuclInt;
    parametrization_ = ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich;

    integral_             = new Integral(IROMB, IMAXS, IPREC);
    integral_for_dEdx_    = new Integral(IROMB, IMAXS, IPREC);

    dndx_integral_.resize(medium_->GetNumComponents());

    for(int i =0 ; i<medium_->GetNumComponents();i++)
    {
        dndx_integral_.at(i) = new Integral(IROMB, IMAXS, IPREC);
    }

    prob_for_component_.resize(medium_->GetNumComponents());
    dedx_interpolant_     = NULL;
    interpolant_measured_ = NULL;

    CalculateWoodSaxonPotential();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Photonuclear::Photonuclear(const Photonuclear &photo)
    :CrossSections                      ( photo )
    ,component_                         ( photo.component_ )
    ,init_measured_                     ( photo.init_measured_ )
    ,init_hardbb_                       ( photo.init_hardbb_ )
    ,hmax_                              ( photo.hmax_ )
    ,v_                                 ( photo.v_ )
    ,do_photo_interpolation_            ( photo.do_photo_interpolation_ )
    ,hard_component_                    ( photo.hard_component_ )
    ,shadow_                            ( photo.shadow_ )
    ,integral_                          ( new Integral(*photo.integral_) )
    ,integral_for_dEdx_                 ( new Integral(*photo.integral_for_dEdx_) )
    ,prob_for_component_                ( photo.prob_for_component_ )
    ,woodSaxonPotential_                ( photo.woodSaxonPotential_ )
{

    if(photo.dedx_interpolant_ != NULL)
    {
        dedx_interpolant_ = new Interpolant(*photo.dedx_interpolant_) ;
    }
    else
    {
        dedx_interpolant_ = NULL;
    }

    if(photo.interpolant_measured_ != NULL)
    {
        interpolant_measured_ = new Interpolant(*photo.interpolant_measured_) ;
    }
    else
    {
        interpolant_measured_ = NULL;
    }

    dndx_integral_.resize( photo.dndx_integral_.size() );
    dndx_interpolant_1d_.resize( photo.dndx_interpolant_1d_.size() );
    dndx_interpolant_2d_.resize( photo.dndx_interpolant_2d_.size() );
    interpolant_hardBB_.resize( photo.interpolant_hardBB_.size() );
    photo_interpolant_.resize( photo.photo_interpolant_.size() );

    for(unsigned int i =0; i<photo.dndx_integral_.size(); i++)
    {
        dndx_integral_.at(i) = new Integral( *photo.dndx_integral_.at(i) );
    }
    for(unsigned int i =0; i<photo.dndx_interpolant_1d_.size(); i++)
    {
        dndx_interpolant_1d_.at(i) = new Interpolant( *photo.dndx_interpolant_1d_.at(i) );
    }
    for(unsigned int i =0; i<photo.dndx_interpolant_2d_.size(); i++)
    {
        dndx_interpolant_2d_.at(i) = new Interpolant( *photo.dndx_interpolant_2d_.at(i) );
    }
    for(unsigned int i =0; i<photo.interpolant_hardBB_.size(); i++)
    {
        interpolant_hardBB_.at(i) = new Interpolant( *photo.interpolant_hardBB_.at(i) );
    }
    for(unsigned int i =0; i<photo.photo_interpolant_.size(); i++)
    {
        photo_interpolant_.at(i) = new Interpolant( *photo.photo_interpolant_.at(i) );
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Photonuclear::Photonuclear(Medium* medium, EnergyCutSettings* cut_settings)
    : CrossSections(medium, cut_settings)
    , component_(0)
    , init_measured_(true)
    , init_hardbb_(true)
    , hmax_(8)
    , v_(0)
    , do_photo_interpolation_(false)
    , hard_component_(false)
    , shadow_(ShadowingType::ButkevichMikhailov)
    , dndx_integral_()
    , interpolant_hardBB_()
    , dndx_interpolant_1d_()
    , dndx_interpolant_2d_()
    , photo_interpolant_()
    , prob_for_component_()
    , woodSaxonPotential_()

{
    name_       = "Photonuclear";
    type_       = ParticleType::NuclInt;
    multiplier_ = 1.;
    parametrization_ = ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich;

    integral_             = new Integral(IROMB, IMAXS, IPREC);
    integral_for_dEdx_    = new Integral(IROMB, IMAXS, IPREC);

    dndx_integral_.resize(medium_->GetNumComponents());

    for(int i =0 ; i<medium_->GetNumComponents();i++){
        dndx_integral_.at(i) = new Integral(IROMB, IMAXS, IPREC);
    }

    prob_for_component_.resize(medium_->GetNumComponents());
    dedx_interpolant_     = NULL;
    interpolant_measured_ = NULL;

    CalculateWoodSaxonPotential();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Photonuclear& Photonuclear::operator=(const Photonuclear &photo)
{

    if (this != &photo)
    {
      Photonuclear tmp(photo);
      swap(tmp);
    }
    return *this;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Photonuclear::operator==(const Photonuclear &photo) const
{
    if( this->CrossSections::operator !=(photo) )                           return false;

    if( component_                  !=  photo.component_)                   return false;
    if( hmax_                       !=  photo.hmax_)                        return false;
    if( v_                          !=  photo.v_)                           return false;
    if( do_photo_interpolation_     !=  photo.do_photo_interpolation_ )     return false;
    if( hard_component_             !=  photo.hard_component_ )             return false;
    if( shadow_                     !=  photo.shadow_ )                     return false;
    if( parametrization_            !=  photo.parametrization_ )            return false;
    if( init_hardbb_                !=  photo.init_hardbb_)                 return false;
    if( init_measured_              !=  photo.init_measured_)               return false;
    if( *integral_                  != *photo.integral_)                    return false;
    if( *integral_for_dEdx_         != *photo.integral_for_dEdx_)           return false;
    if( prob_for_component_.size()  !=  photo.prob_for_component_.size())   return false;
    if( woodSaxonPotential_.size()  !=  photo.woodSaxonPotential_.size())   return false;
    if( dndx_integral_.size()       !=  photo.dndx_integral_.size())        return false;
    if( dndx_interpolant_1d_.size() !=  photo.dndx_interpolant_1d_.size())  return false;
    if( dndx_interpolant_2d_.size() !=  photo.dndx_interpolant_2d_.size())  return false;
    if( interpolant_hardBB_.size()  !=  photo.interpolant_hardBB_.size())   return false;
    if( photo_interpolant_.size()   !=  photo.photo_interpolant_.size())    return false;


    for(unsigned int i =0; i<photo.dndx_integral_.size(); i++)
    {
        if( *dndx_integral_.at(i)       != *photo.dndx_integral_.at(i) )        return false;
    }
    for(unsigned int i =0; i<photo.dndx_interpolant_1d_.size(); i++)
    {
        if( *dndx_interpolant_1d_.at(i) != *photo.dndx_interpolant_1d_.at(i) )  return false;
    }
    for(unsigned int i =0; i<photo.dndx_interpolant_2d_.size(); i++)
    {
        if( *dndx_interpolant_2d_.at(i) != *photo.dndx_interpolant_2d_.at(i) )  return false;
    }
    for(unsigned int i =0; i<photo.prob_for_component_.size(); i++)
    {
        if( prob_for_component_.at(i) != photo.prob_for_component_.at(i) )      return false;
    }
    for(unsigned int i =0; i<photo.woodSaxonPotential_.size(); i++)
    {
        if( woodSaxonPotential_.at(i) != photo.woodSaxonPotential_.at(i) )      return false;
    }
    for(unsigned int i =0; i<photo.interpolant_hardBB_.size(); i++)
    {
        if( *interpolant_hardBB_.at(i) != *photo.interpolant_hardBB_.at(i) )    return false;
    }
    for(unsigned int i =0; i<photo.photo_interpolant_.size(); i++)
    {
        if( *photo_interpolant_.at(i) != *photo.photo_interpolant_.at(i) )      return false;
    }

    if( dedx_interpolant_ != NULL && photo.dedx_interpolant_ != NULL)
    {
        if( *dedx_interpolant_   != *photo.dedx_interpolant_)                   return false;
    }
    else if( dedx_interpolant_ != photo.dedx_interpolant_)                      return false;

    if( interpolant_measured_ != NULL && photo.interpolant_measured_ != NULL)
    {
        if( *interpolant_measured_   != *photo.interpolant_measured_)           return false;
    }
    else if( interpolant_measured_ != photo.interpolant_measured_)              return false;

    //else
    return true;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Photonuclear::operator!=(const Photonuclear &photo) const
{
    return !(*this == photo);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Photonuclear::swap(Photonuclear &photo)
{
    using std::swap;

    this->CrossSections::swap(photo);

    swap(component_,photo.component_);
    swap(init_measured_,photo.init_measured_);
    swap(init_hardbb_,photo.init_hardbb_);
    swap(hmax_,photo.hmax_);
    swap(v_,photo.v_);
    swap(do_photo_interpolation_, photo.do_photo_interpolation_);
    swap(hard_component_ , photo.hard_component_ );
    swap(shadow_ , photo.shadow_ );

    integral_for_dEdx_->swap(*photo.integral_for_dEdx_);
    integral_->swap(*photo.integral_);

    prob_for_component_.swap(photo.prob_for_component_);
    woodSaxonPotential_.swap(photo.woodSaxonPotential_);
    dndx_integral_.swap(photo.dndx_integral_);
    dndx_interpolant_1d_.swap(photo.dndx_interpolant_1d_);
    dndx_interpolant_2d_.swap(photo.dndx_interpolant_2d_);
    interpolant_hardBB_.swap(photo.interpolant_hardBB_);
    photo_interpolant_.swap(photo.photo_interpolant_);

    if( dedx_interpolant_ != NULL && photo.dedx_interpolant_ != NULL)
    {
        dedx_interpolant_->swap(*photo.dedx_interpolant_);
    }
    else if( dedx_interpolant_ == NULL && photo.dedx_interpolant_ != NULL)
    {
        dedx_interpolant_ = new Interpolant(*photo.dedx_interpolant_);
        photo.dedx_interpolant_ = NULL;
    }
    else if( dedx_interpolant_ != NULL && photo.dedx_interpolant_ == NULL)
    {
        photo.dedx_interpolant_ = new Interpolant(*dedx_interpolant_);
        dedx_interpolant_ = NULL;
    }


    if( interpolant_measured_ != NULL && photo.interpolant_measured_ != NULL)
    {
        interpolant_measured_->swap(*photo.interpolant_measured_);
    }
    else if( interpolant_measured_ == NULL && photo.interpolant_measured_ != NULL)
    {
        interpolant_measured_ = new Interpolant(*photo.interpolant_measured_);
        photo.interpolant_measured_ = NULL;
    }
    else if( interpolant_measured_ != NULL && photo.interpolant_measured_ == NULL)
    {
        photo.interpolant_measured_ = new Interpolant(*interpolant_measured_);
        interpolant_measured_ = NULL;
    }


}

namespace PROPOSAL
{

ostream& operator<<(std::ostream& os, Photonuclear const &photo)
{
    os<<"---------------------------Photonuclear( "<<&photo<<" )---------------------------"<<std::endl;
    os<< static_cast <const CrossSections &>( photo ) << endl;
    os<< "------- Class Specific: " << endl;
    os<< "\tinit_measured:\t\t" << photo.init_measured_ << endl;
    os<< "\tinit_hardbb:\t\t" << photo.init_hardbb_ << endl;
    os<< "\tshadow:\t\t\t" << photo.shadow_ << endl;
    os<< "\thard_component:\t\t" << photo.hard_component_ << endl;
    os<< "\tparametrization_:\t" << photo.parametrization_ << endl;
    os<<endl;
    os<<"\tintegral:\t\t"<<photo.integral_ << endl;
    os<<"\tdedx_integral:\t"<<photo.integral_for_dEdx_ << endl;
    os<<"\tdndx_integral:\t"<<photo.dndx_integral_.size()<<endl;
    for(unsigned int i=0;i<photo.dndx_integral_.size();i++)
    {
        os<<"\t\tadress:\t\t"<<photo.dndx_integral_.at(i)<<endl;
    }
    os<<endl;
    os<<"\tdo_dedx_Interpolation:\t\t"<<photo.do_dedx_Interpolation_<<endl;
    os<<"\tdedx_interpolant:\t\t"<<photo.dedx_interpolant_<<endl;
    os<<endl;
    os<<"\tinterpolant_measured:\t\t"<<photo.interpolant_measured_<<endl;
    os<<"\tinterpolant_hardBB:\t\t"<<photo.interpolant_hardBB_.size()<<endl;
    for(unsigned int i=0;i<photo.interpolant_hardBB_.size();i++)
    {
        os<<"\t\tadress:\t\t"<<photo.interpolant_hardBB_.at(i)<<endl;
    }
    os<<endl;
    os<<"\tdo_dndx_Interpolation:\t\t"<<photo.do_dndx_Interpolation_<<endl;
    os<<"\tdndx_interpolant_1d:\t\t"<<photo.dndx_interpolant_1d_.size()<<endl;
    for(unsigned int i=0;i<photo.dndx_interpolant_1d_.size();i++)
    {
        os<<"\t\tadress:\t\t"<<photo.dndx_interpolant_1d_.at(i)<<endl;
    }
    os<<"\tdndx_interpolant_2d:\t\t"<<photo.dndx_interpolant_2d_.size()<<endl;
    for(unsigned int i=0;i<photo.dndx_interpolant_2d_.size();i++)
    {
        os<<"\t\tadress:\t\t"<<photo.dndx_interpolant_2d_.at(i)<<endl;
    }
    os<<endl;
    os<<"\tdo_photo_interpolation:\t"<<photo.do_photo_interpolation_<<endl;
    os<<"\tphoto_interpolant:\t\t"<<photo.photo_interpolant_.size()<<endl;
    for(unsigned int i=0;i<photo.photo_interpolant_.size();i++)
    {
        os<<"\t\tadress:\t\t"<<photo.photo_interpolant_.at(i)<<endl;
    }
    os<<endl;
    os<<"\tTemp. variables: " << endl;
    os<< "\t\tcomponent:\t\t" << photo.component_<<endl;
    os<< "\t\thmax:\t\t" << photo.hmax_<<endl;
    os<< "\t\tv:\t\t" << photo.v_<<endl;
    os<<"\t\tprob_for_component:\t"<<photo.prob_for_component_.size()<<endl;
    for(unsigned int i=0;i<photo.prob_for_component_.size();i++)
    {
        os<<"\t\t\tvalue:\t\t"<<photo.prob_for_component_.at(i)<<endl;
    }
    os<<"\t\twoodSaxonPotential:\t"<<photo.woodSaxonPotential_.size()<<endl;
    for(unsigned int i=0;i<photo.woodSaxonPotential_.size();i++)
    {
        os<<"\t\t\tvalue:\t\t"<<photo.woodSaxonPotential_.at(i)<<endl;
    }
    os<<"-----------------------------------------------------------------------------------------------";
    return os;
}

}  // namespace PROPOSAL

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-----------------Photonuclear Parametrizations-------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::ParametrizationOfRealPhotonAssumption(const PROPOSALParticle& particle, double v, int i)
{
    double aux, aum, k, G, t, sgn;
    double nu = v*particle.GetEnergy()*1.e-3;

    switch (parametrization_)
    {
        case ParametrizationType::PhotoKokoulinShadowBezrukovSoft:
        case ParametrizationType::PhotoKokoulinShadowBezrukovHard:
            sgn = PhotoNucleusCrossSectionKokoulin(nu);
            break;
        case ParametrizationType::PhotoRhodeShadowBezrukovSoft:
        case ParametrizationType::PhotoRhodeShadowBezrukovHard:
            sgn = PhotoNucleusCrossSectionRhode(nu);
            break;
        case ParametrizationType::PhotoBezrukovBugaevShadowBezrukovSoft:
        case ParametrizationType::PhotoBezrukovBugaevShadowBezrukovHard:
            sgn = PhotoNucleusCrossSectionBezrukovBugaev(nu);
            break;
        case ParametrizationType::PhotoZeusShadowBezrukovSoft:
        case ParametrizationType::PhotoZeusShadowBezrukovHard:
            sgn = PhotoNucleusCrossSectionZeus(nu, medium_->GetComponents().at(i)->GetAverageNucleonWeight());
            break;
        default:
            log_fatal("The photonuclear Parametrization %i is not supported.\n", parametrization_);
    }

    const double m1 =   0.54;
    const double m2 =   1.80;

    double particle_charge = particle.GetCharge();
    double particle_mass = particle.GetMass();
    double atomic_number = medium_->GetComponents().at(i)->GetAtomicNum();

    k   =   1 - 2/v + 2/(v*v);

    if(medium_->GetComponents().at(i)->GetNucCharge()==1)
    {
        G   =   1;
    }
    else
    {
        G   =   ShadowBezrukovBugaev(sgn, atomic_number);
    }

    G       *=  3;
    aux     =   v*particle_mass*1.e-3;
    t       =   aux*aux/(1-v);
    aum     =   particle_mass*1.e-3;
    aum     *=  aum;
    aux     =   2*aum/t;
    aux     =   G*((k + 4*aum/m1)*log(1 + m1/t) - (k*m1)/(m1 + t) - aux)
                + ((k + 2*aum/m2)*log(1 + m2/t) - aux)
                + aux*((G*(m1 - 4*t))/(m1 + t) + (m2/t)*log(1 + t/m2));

    aux     *=  ALPHA/(8*PI)*atomic_number*v*sgn*1.e-30;

    if(hard_component_)
    {
        if (particle.getHardBB() != NULL)
        {
                aux +=  atomic_number*1.e-30*HardBB(particle, v);
        }
        // switch (particle.getHardBB().empty())
        // {
        //     case ParticleType::MuMinus:
        //     case ParticleType::MuPlus:
        //     case ParticleType::TauMinus:
        //     case ParticleType::TauPlus:
        //         aux +=  atomic_number*1.e-30*HardBB(particle.GetEnergy(), v);
        //         break;
        //     default:
        //         break;
        // }
    }

    return medium_->GetMolDensity()*medium_->GetComponents().at(i)->GetAtomInMolecule()*particle_charge*particle_charge*aux;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Photonuclear::ParametrizationOfQ2Integration(const PROPOSALParticle& particle, double v, int i)
{
    double particle_energy = particle.GetEnergy();

    CrossSections::IntegralLimits limits = SetIntegralLimits(particle, i);

    if(do_photo_interpolation_)
    {
        if(v >= limits.vUp)
        {
            return max(photo_interpolant_.at(i)->Interpolate(particle_energy, log(v / limits.vUp)/log(limits.vMax / limits.vUp)), 0.0);
        }
    }

    double aux, min, max;
    double particle_mass = particle.GetMass();
    double particle_charge = particle.GetCharge();

    component_ =   i;
    v_         =   v;
    min        =   particle_mass*v;
    min        *=  min/(1-v);

    if(particle_mass < MPI)
    {
        aux     =   particle_mass*particle_mass/particle_energy;
        min     -=  (aux*aux)/(2*(1-v));
    }

    max =   2*medium_->GetComponents().at(i)->GetAverageNucleonWeight()*particle_energy*(v-limits.vMin);

    //  if(form==4) max=Math.min(max, 5.5e6);  // as requested in Butkevich and Mikheyev
    if(min > max)
    {
        return 0;
    }

    aux = medium_->GetMolDensity()*medium_->GetComponents().at(i)->GetAtomInMolecule()*particle_charge*particle_charge;

    switch (parametrization_)
    {
        case ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowDutta:
        case ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowButkevich:
            aux *= integral_->Integrate(min, max, boost::bind(&Photonuclear::FunctionToIntegralALLM91, this, boost::cref(particle), _1),4);
            break;
        case ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowDutta:
        case ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich:
            aux *= integral_->Integrate(min, max, boost::bind(&Photonuclear::FunctionToIntegralALLM97, this, boost::cref(particle), _1),4);
            break;
        case ParametrizationType::PhotoButkevichMikhailovShadowDutta:
        case ParametrizationType::PhotoButkevichMikhailovShadowButkevich:
            aux *= integral_->Integrate(min, max, boost::bind(&Photonuclear::FunctionToIntegralButMik, this, boost::cref(particle), _1),4);
            break;
        case ParametrizationType::PhotoRenoSarcevicSuShadowDutta:
        case ParametrizationType::PhotoRenoSarcevicSuShadowButkevich:
            aux *= integral_->Integrate(min, max, boost::bind(&Photonuclear::FunctionToIntegralRSS, this, boost::cref(particle), _1),4);
            break;
        default:
            log_fatal("The photonuclear Parametrization %i is not supported.\n", parametrization_);
    }

    return aux;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------Cross section / limit / private-----------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::CalculateStochasticLoss(const PROPOSALParticle& particle, double rnd)
{
    double rand, rsum, particle_energy;

    particle_energy = particle.GetEnergy();

    rand    =   rnd*sum_of_rates_;
    rsum    =   0;

    for(int i=0; i<(medium_->GetNumComponents()); i++)
    {
        rsum    += prob_for_component_.at(i);

        if(rsum > rand)
        {
            if(do_dndx_Interpolation_)
            {
                CrossSections::IntegralLimits limits = SetIntegralLimits(particle, i);

                if(limits.vUp == limits.vMax)
                {
                    return particle_energy * limits.vUp;
                }

                return particle_energy*(limits.vUp * exp(dndx_interpolant_2d_.at(i)->FindLimit(particle_energy, rnd_ * prob_for_component_.at(i))*log(limits.vMax / limits.vUp)));

            }
            else
            {
                component_ = i;
                return (particle_energy)*dndx_integral_.at(i)->GetUpperLimit();

            }
        }
    }

    //sometime everything is fine, just the probability for interaction is zero
    bool prob_for_all_comp_is_zero=true;
    for(int i=0; i<(medium_->GetNumComponents()); i++)
    {
        CrossSections::IntegralLimits limits = SetIntegralLimits(particle, i);
        if(limits.vUp != limits.vMax)
            prob_for_all_comp_is_zero=false;
    }
    if(prob_for_all_comp_is_zero)return 0;

    log_fatal("sum was not initialized correctly");

    return 0;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


CrossSections::IntegralLimits Photonuclear::SetIntegralLimits(const PROPOSALParticle& particle, int component)
{

    double aux;

    component_ = component;
    IntegralLimits limits;

    limits.vMin    =   (MPI + (MPI*MPI)/(2*medium_->GetComponents().at(component_)->GetAverageNucleonWeight()))/particle.GetEnergy();

    if(particle.GetMass() < MPI)
    {
        aux     =   particle.GetMass()/medium_->GetComponents().at(component_)->GetAverageNucleonWeight();
        limits.vMax    =   1 - medium_->GetComponents().at(component_)->GetAverageNucleonWeight()*(1 + aux*aux)/(2*particle.GetEnergy());
    }
    else
    {
        limits.vMax    =   1;
    }

    limits.vMax    =   min(limits.vMax, 1-particle.GetMass()/particle.GetEnergy());

    if(limits.vMax < limits.vMin)
    {
        limits.vMax    =   limits.vMin;
    }

    limits.vUp     =   min(limits.vMax, cut_settings_->GetCut(particle.GetEnergy()));

    if(limits.vUp < limits.vMin)
    {
        limits.vUp =   limits.vMin;
    }

    return limits;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::PhotoN(const PROPOSALParticle& particle, double v, int i)
{
    switch(parametrization_)
    {
        case ParametrizationType::PhotoKokoulinShadowBezrukovSoft:
        case ParametrizationType::PhotoKokoulinShadowBezrukovHard:
        case ParametrizationType::PhotoRhodeShadowBezrukovSoft:
        case ParametrizationType::PhotoRhodeShadowBezrukovHard:
        case ParametrizationType::PhotoBezrukovBugaevShadowBezrukovSoft:
        case ParametrizationType::PhotoBezrukovBugaevShadowBezrukovHard:
        case ParametrizationType::PhotoZeusShadowBezrukovSoft:
        case ParametrizationType::PhotoZeusShadowBezrukovHard:
            return ParametrizationOfRealPhotonAssumption(particle, v, i);
        case ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowDutta:
        case ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowButkevich:
        case ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowDutta:
        case ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich:
        case ParametrizationType::PhotoButkevichMikhailovShadowDutta:
        case ParametrizationType::PhotoButkevichMikhailovShadowButkevich:
        case ParametrizationType::PhotoRenoSarcevicSuShadowDutta:
        case ParametrizationType::PhotoRenoSarcevicSuShadowButkevich:
            return ParametrizationOfQ2Integration(particle, v, i);
        default:
            log_fatal("The photonuclear parametrization_ '%i' is not supported! Be careful 0 is returned. \n",parametrization_);
            return 0;
    }

}


//----------------------------------------------------------------------------//
//-----------Photon Nucleus Cross Section Parametrizations--------------------//
//----------------------------------------------------------------------------//


double Photonuclear::PhotoNucleusCrossSectionCaldwell(double nu)
{
    return 49.2 + 11.1*log(nu) + 151.8/sqrt(nu);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::PhotoNucleusCrossSectionKokoulin(double nu)
{
    if(nu<=200.)
    {
        if(nu<=17.)
        {
            return 96.1+82./sqrt(nu);
        }
        else
        {
            return PhotoNucleusCrossSectionBezrukovBugaev(nu);
        }
    }
    else
    {
        return PhotoNucleusCrossSectionCaldwell(nu);
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::PhotoNucleusCrossSectionRhode(double nu)
{
    if(nu<=200.)
    {
        return MeasuredSgN(nu);
    }
    else
    {
        return PhotoNucleusCrossSectionCaldwell(nu);
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::PhotoNucleusCrossSectionBezrukovBugaev(double nu)
{
    double aux;

    aux =   log(0.0213*nu);
    aux =   114.3 + 1.647*aux*aux;

    return aux;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::PhotoNucleusCrossSectionZeus(double nu, double medium_average_nucleon_weight)
{
    double aux;

    aux =   nu*2.e-3*medium_average_nucleon_weight;
    aux =   63.5*pow(aux, 0.097) + 145*pow(aux , -0.5);

    return aux;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Photonuclear::SetMeasured()
{
    if(init_measured_)
    {
        init_measured_=false;

        double x_aux[]  =   {0, 0.1, 0.144544, 0.20893, 0.301995,
                            0.436516, 0.630957, 0.912011, 1.31826, 1.90546,
                            2.75423, 3.98107, 5.7544, 8.31764, 12.0226,
                            17.378, 25.1189, 36.3078, 52.4807, 75.8577,
                            109.648, 158.489, 229.087, 331.131, 478.63,
                            691.831, 1000, 1445.44, 2089.3, 3019.95,
                            4365.16, 6309.58, 9120.12, 13182.6, 19054.6,
                            27542.3, 39810.8, 57544, 83176.4, 120226,
                            173780, 251188, 363078, 524807, 758576,
                            1.09648e+06, 1.58489e+06, 2.29086e+06, 3.3113e+06, 4.78628e+06,
                            6.91828e+06, 9.99996e+06};

        double y_aux[]  =   {0, 0.0666667, 0.0963626, 159.74, 508.103,
                            215.77, 236.403, 201.919, 151.381, 145.407,
                            132.096, 128.546, 125.046, 121.863, 119.16,
                            117.022, 115.496, 114.607, 114.368, 114.786,
                            115.864, 117.606, 120.011, 123.08, 126.815,
                            131.214, 136.278, 142.007, 148.401, 155.46,
                            163.185, 171.574, 180.628, 190.348, 200.732,
                            211.782, 223.497, 235.876, 248.921, 262.631,
                            277.006, 292.046, 307.751, 324.121, 341.157,
                            358.857, 377.222, 396.253, 415.948, 436.309,
                            457.334, 479.025};

        vector<double> x(x_aux, x_aux + sizeof(x_aux) / sizeof(double) );
        vector<double> y(y_aux, y_aux + sizeof(y_aux) / sizeof(double) );

        interpolant_measured_   =   new Interpolant(x, y, 4, false, false);
    }

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::MeasuredSgN(double e)
{
    SetMeasured();
    return interpolant_measured_->InterpolateArray(e);
}


//----------------------------------------------------------------------------//
//-------------------------Hard component-------------------------------------//
//----------------------------------------------------------------------------//


void Photonuclear::EnableHardBB(const PROPOSALParticle& particle)
{
    if(init_hardbb_)
    {
        init_hardbb_           =   false;
        double x_aux[]  =   {3, 4, 5, 6, 7, 8, 9};


        // double y_aux[][56][7]=
        // {
        //     {
        //         {7.174409e-4, 1.7132e-3, 4.082304e-3, 8.628455e-3, 0.01244159, 0.02204591, 0.03228755},
        //         {-0.2436045, -0.5756682, -1.553973, -3.251305, -5.976818, -9.495636, -13.92918},
        //         {-0.2942209, -0.68615, -2.004218, -3.999623, -6.855045, -10.05705, -14.37232},
        //         {-0.1658391, -0.3825223, -1.207777, -2.33175, -3.88775, -5.636636, -8.418409},
        //         {-0.05227727, -0.1196482, -0.4033373, -0.7614046, -1.270677, -1.883845, -2.948277},
        //         {-9.328318e-3, -0.02124577, -0.07555636, -0.1402496, -0.2370768, -0.3614146, -0.5819409},
        //         {-8.751909e-4, -1.987841e-3, -7.399682e-3, -0.01354059, -0.02325118, -0.03629659, -0.059275},
        //         {-3.343145e-5, -7.584046e-5, -2.943396e-4, -5.3155e-4, -9.265136e-4, -1.473118e-3, -2.419946e-3}
        //     },
        //     {
        //         {-1.269205e-4, -2.843877e-4, -5.761546e-4, -1.195445e-3, -1.317386e-3, -9.689228e-15, -6.4595e-15},
        //         {-0.01563032, -0.03589573, -0.07768545, -0.157375, -0.2720009, -0.4186136, -0.8045046},
        //         {0.04693954, 0.1162945, 0.3064255, 0.7041273, 1.440518, 2.533355, 3.217832},
        //         {0.05338546, 0.130975, 0.3410341, 0.7529364, 1.425927, 2.284968, 2.5487},
        //         {0.02240132, 0.05496, 0.144945, 0.3119032, 0.5576727, 0.8360727, 0.8085682},
        //         {4.658909e-3, 0.01146659, 0.03090286, 0.06514455, 0.1109868, 0.1589677, 0.1344223},
        //         {4.822364e-4, 1.193018e-3, 3.302773e-3, 6.843364e-3, 0.011191, 0.015614, 0.01173827},
        //         {1.9837e-5, 4.940182e-5, 1.409573e-4, 2.877909e-4, 4.544877e-4, 6.280818e-4, 4.281932e-4}
        //     }
        // };

        //TODO(mario): Must be cleared Tue 2017/08/08
        interpolant_hardBB_.resize(hmax_);

        const HardBBTables::VecType* y = particle.getHardBB();

        if (y != NULL)
        {
            for(int i=0; i < hmax_; i++)
            {
                vector<double> x(x_aux, x_aux + sizeof(x_aux) / sizeof(double) );
                interpolant_hardBB_.at(i) = new Interpolant(x, y->at(i), 4, false, false);
            }
        }

        // for(int i=0; i < hmax_; i++)
        // {
        //     vector<double> x(x_aux, x_aux + sizeof(x_aux) / sizeof(double) );
        //
        //     if (particle_->GetType() == ParticleType::MuMinus || particle_->GetType() == ParticleType::MuPlus)
        //     {
        //         vector<double> y(y_aux[0][i], y_aux[0][i] + sizeof(y_aux[0][i]) / sizeof(double) );
        //         interpolant_hardBB_.at(i)    =   new Interpolant(x, y, 4, false, false) ;
        //     }
        //     else if (particle_->GetType() == ParticleType::TauMinus || particle_->GetType() == ParticleType::TauPlus)
        //     {
        //         vector<double> y(y_aux[1][i], y_aux[1][i] + sizeof(y_aux[1][i]) / sizeof(double) );
        //         interpolant_hardBB_.at(i)    =   new Interpolant(x, y, 4, false, false) ;
        //     }
        //     else
        //     {
        //         vector<double> y(y_aux[0][i], y_aux[0][i] + sizeof(y_aux[0][i]) / sizeof(double) );
        //         interpolant_hardBB_.at(i)    =   new Interpolant(x, y, 4, false, false) ;
        //         log_warn("No hard component paraemtrization found for particle %s! Parametrization for muons are used. Be careful!", particle_->GetName().c_str());
        //     }
        // }
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::HardBB(const PROPOSALParticle& particle, double v)
{
    EnableHardBB(particle);
    double e = particle.GetEnergy();

    if(e<1.e5 || v<1.e-7)
    {
        return 0;
    }

    double aux, sum, lov, loe;

    sum =   0;
    aux =   1;
    lov =   log(v)/LOG10;
    loe =   log(e)/LOG10-3;

    for(int i=0; i < hmax_; i++)
    {
        if(i>0)
        {
            aux *=  lov;
        }

        sum +=  aux*interpolant_hardBB_.at(i)->InterpolateArray(loe);
    }
    return sum/v;

}


//----------------------------------------------------------------------------//
//-------------------------Shadow Parametrizations----------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::ShadowBezrukovBugaev(double sgn, double atomic_number)
{
    double G, aux;

    aux =   0.00282*pow(atomic_number, 1./3)*sgn;
    G   =   (3/aux)*(0.5 + ((1 + aux)*exp(-aux) - 1)/(aux*aux));

    return G;
}

double Photonuclear::ShadowEffect(double x , double nu)
{
    if(medium_->GetComponents().at(component_)->GetNucCharge()==1) return 1;

    double G, atomic_number;

    atomic_number = medium_->GetComponents().at(component_)->GetAtomicNum();

    if(shadow_ == ShadowingType::Dutta)
    {
        if(x<0.0014)
        {
            G   =   pow(atomic_number, -0.1);
        }
        else if(x<0.04)
        {
            G   =   pow(atomic_number, 0.069*log(x)/LOG10+0.097);
        }
        else
        {
            G   =   1;
        }
    }
    else if (shadow_ == ShadowingType::ButkevichMikhailov)
    {
        if(x>0.3)
        {
            const double Mb =   0.437;
            const double la =   0.5;
            const double x2 =   0.278;

            double mb, Aosc, mu, au, ac;

            mb      =   Mb*woodSaxonPotential_.at(component_);
            au      =   1/(1 - x);
            ac      =   1/(1 - x2);
            mu      =   MPI/medium_->GetComponents().at(component_)->GetAverageNucleonWeight();
            Aosc    =   (1 - la*x)*((au - ac)-mu*(au*au - ac*ac));
            G       =   1 - mb*Aosc;
        }
        else
        {
            const double M1 =   0.129;
            const double M2 =   0.456;
            const double M3 =   0.553;

            double m1, m2, m3, x0, sgn;

            m1  =   M1*woodSaxonPotential_.at(component_);
            m2  =   M2*woodSaxonPotential_.at(component_);
            m3  =   M3*woodSaxonPotential_.at(component_);
            nu  *=  1.e-3;
            sgn =   112.2*(0.609*pow(nu, 0.0988) + 1.037*pow(nu, -0.5944));
            G   =   ShadowBezrukovBugaev(sgn, atomic_number);
            G   =   0.75*G + 0.25;
            x0  =   pow(G/(1+m2), 1/m1);

            if(x>=x0)
            {
                G   =   pow(x, m1)*(1+m2)*(1-m3*x);
            }
        }
    }
    else
    {
        log_warn("shadow type '%i' is not valid must be Dutta or Butkevich, other values are not supported! \
            Be careful shadow effect factor is set to 1!", shadow_);

        G   =   1.;
    }

    return G;
}


//----------------------------------------------------------------------------//
//--------------Wood-Saxon (for Butkevich-Mikheyev Shadowing)-----------------//
//----------------------------------------------------------------------------//


void Photonuclear::CalculateWoodSaxonPotential()
{
    woodSaxonPotential_.resize( medium_->GetNumComponents() );

    for(int i=0; i<(medium_->GetNumComponents()); i++)
    {
        if (medium_->GetComponents().at(i)->GetNucCharge() != 1.0)
        {
            Integral integral(IROMB, IMAXS, IPREC);

            double r0;
            r0 = pow(medium_->GetComponents().at(i)->GetAtomicNum(), 1.0 / 3.0);
            r0 = 1.12 * r0 - 0.86 / r0;

            woodSaxonPotential_.at(i) = 1.0 -
                  4.0 * PI * 0.17 *
                      integral.Integrate(
                          r0, -1.0,
                          boost::bind(&Photonuclear::FunctionToWoodSaxonPotentialIntegral, this, r0, _1),
                          3, 2.0) /
                      medium_->GetComponents().at(i)->GetAtomicNum();
        }
        else
        {
            woodSaxonPotential_.at(i) = 0.;
        }
    }
}

// ------------------------------------------------------------------------- //

double Photonuclear::FunctionToWoodSaxonPotentialIntegral(const double r0, double r)
{
    const double a = 0.54;

    return r * r / (1 + exp((r - r0) / a));
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------Functions to interpolate---------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::FunctionToBuildDEdxInterpolant(const PROPOSALParticle& particle,  double energy)
{
    PROPOSALParticle tmp_particle(particle);
    tmp_particle.SetEnergy(energy);

    return CalculatedEdx(tmp_particle);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::FunctionToBuildDNdxInterpolant1D(double energy)
{
    return dndx_interpolant_2d_.at(component_)->Interpolate(energy, 1.);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::FunctionToBuildDNdxInterpolant2D(const PROPOSALParticle& particle, double energy, double v)
{
    PROPOSALParticle tmp_particle(particle);
    tmp_particle.SetEnergy(energy);

    CrossSections::IntegralLimits limits = SetIntegralLimits(tmp_particle, component_);
    if(limits.vUp == limits.vMax)
    {
        return 0;
    }

    v   =   limits.vUp*exp(v*log(limits.vMax / limits.vUp));

    return dndx_integral_.at(component_)->Integrate(limits.vUp, v, boost::bind(&Photonuclear::FunctionToDNdxIntegral, this, boost::cref(tmp_particle), _1) ,4);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::FunctionToBuildPhotoInterpolant(const PROPOSALParticle& particle, double energy, double v)
{
    PROPOSALParticle tmp_particle(particle);
    tmp_particle.SetEnergy(energy);

    CrossSections::IntegralLimits limits = SetIntegralLimits(particle, component_);

    if(limits.vUp == limits.vMax)
    {
        return 0;
    }

    v   =   limits.vUp*exp(v*log(limits.vMax / limits.vUp));

    return PhotoN(particle, v, component_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------Functions to integrate----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::FunctionToDEdxIntegral(const PROPOSALParticle& particle, double variable)
{
    return variable*PhotoN(particle, variable, component_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::FunctionToDNdxIntegral(const PROPOSALParticle& particle, double variable)
{
    return multiplier_*PhotoN(particle, variable, component_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::FunctionToIntegralALLM91(const PROPOSALParticle& particle, double Q2)
{
    Components::Component* component = medium_->GetComponents().at(component_);

    double x, aux, nu, G, F2, R2;

    nu  =   v_*particle.GetEnergy();
    x   =   Q2/(2*medium_->GetComponents().at(component_)->GetAverageNucleonWeight()*nu);

    G = ShadowEffect(x , nu);


    double P, W2;

    aux =   x*x;
    P   =   1 - 1.85*x + 2.45*aux - 2.35*aux*x + aux*aux;
    G   *=  (component->GetNucCharge() + (component->GetAtomicNum() - component->GetNucCharge())*P);
    W2  =   component->GetAverageNucleonWeight()*component->GetAverageNucleonWeight()
            - Q2 + 2*component->GetAverageNucleonWeight()*particle.GetEnergy()*v_;

    double cp1, cp2, cp3, cr1, cr2, cr3, ap1, ap2, ap3, ar1, ar2, ar3;
    double bp1, bp2, bp3, br1, br2, br3, m2o, m2r, L2, m2p, Q2o;


    cp1     =   0.26550;
    cp2     =   0.04856;
    cp3     =   1.04682;
    cr1     =   0.67639;
    cr2     =   0.49027;
    cr3     =   2.66275;
    ap1     =   -0.04503;
    ap2     =   -0.36407;
    ap3     =   8.17091;
    ar1     =   0.60408;
    ar2     =   0.17353;
    ar3     =   1.61812;
    bp1     =   0.49222;
    bp2     =   0.52116;
    bp3     =   3.55115;
    br1     =   1.26066;
    br2     =   1.83624;
    br3     =   0.81141;
    m2o     =   0.30508;
    m2r     =   0.20623;
    L2      =   0.06527;
    m2p     =   10.67564;
    Q2o     =   0.27799;


    // GeV -> MeV conversion
    m2o     *=  1e6;
    m2r     *=  1e6;
    L2      *=  1e6;
    m2p     *=  1e6;
    Q2o     *=  1e6;

    // these values are corrected according to the file f2allm.f from Halina Abramowicz
    bp1     *=  bp1;
    bp2     *=  bp2;
    br1     *=  br1;
    br2     *=  br2;
    Q2o     +=  L2;

    const double R  =   0;

    double cr, ar, cp, ap, br, bp, t;

    t   =   log(log((Q2 + Q2o)/L2)/log(Q2o/L2));

    if(t<0)
    {
        t=0;
    }

    cr  =   cr1 + cr2*pow(t, cr3);
    ar  =   ar1 + ar2*pow(t, ar3);
    cp  =   cp1 + (cp1 - cp2)*(1/(1 + pow(t, cp3)) - 1);
    ap  =   ap1 + (ap1 - ap2)*(1/(1 + pow(t, ap3)) - 1);
    br  =   br1 + br2*pow(t, br3);
    bp  =   bp1 + bp2*pow(t, bp3);

    double xp, xr, F2p, F2r;

    xp  =   (Q2 + m2p)/(Q2 + m2p + W2 - component->GetAverageNucleonWeight()*component->GetAverageNucleonWeight());
    xr  =   (Q2 + m2r)/(Q2 + m2r + W2 - component->GetAverageNucleonWeight()*component->GetAverageNucleonWeight());
    F2p =   cp*pow(xp, ap)*pow(1 - x, bp);
    F2r =   cr*pow(xr, ar)*pow(1 - x, br);
    F2  =   (Q2/(Q2 + m2o))*(F2p + F2r)*G;
    R2  =   (2*(1 + R));



    aux =   ME*RE/Q2;
    aux *=  aux*(1 - v_ - component->GetAverageNucleonWeight()*x*v_/(2*particle.GetEnergy()) +
                 (1 - 2*particle.GetMass()*particle.GetMass()/Q2)
                 *v_*v_*(1 + 4*component->GetAverageNucleonWeight()*component->GetAverageNucleonWeight()*x*x/Q2)/R2);

    return (4*PI*F2/v_)*aux;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::FunctionToIntegralALLM97(const PROPOSALParticle& particle, double Q2)
{
    Components::Component* component = medium_->GetComponents().at(component_);

    double x, aux, nu, G, F2, R2;

    nu  =   v_*particle.GetEnergy();
    x   =   Q2/(2*component->GetAverageNucleonWeight()*nu);

    G = ShadowEffect(x , nu);

    double P, W2;

    aux =   x*x;
    P   =   1 - 1.85*x + 2.45*aux - 2.35*aux*x + aux*aux;
    G   *=  (component->GetNucCharge() + (component->GetAtomicNum() - component->GetNucCharge())*P);
    W2  =   component->GetAverageNucleonWeight()*component->GetAverageNucleonWeight()
            - Q2 + 2*component->GetAverageNucleonWeight()*particle.GetEnergy()*v_;

    double cp1, cp2, cp3, cr1, cr2, cr3, ap1, ap2, ap3, ar1, ar2, ar3;
    double bp1, bp2, bp3, br1, br2, br3, m2o, m2r, L2, m2p, Q2o;


    cp1     =   0.28067;
    cp2     =   0.22291;
    cp3     =   2.1979;
    cr1     =   0.80107;
    cr2     =   0.97307;
    cr3     =   3.4942;
    ap1     =   -0.0808;
    ap2     =   -0.44812;
    ap3     =   1.1709;
    ar1     =   0.58400;
    ar2     =   0.37888;
    ar3     =   2.6063;
    bp1     =   0.60243;
    bp2     =   1.3754;
    bp3     =   1.8439;
    br1     =   0.10711;
    br2     =   1.9386;
    br3     =   0.49338;
    m2o     =   0.31985;
    m2r     =   0.15052;
    L2      =   0.06527;
    m2p     =   49.457;
    Q2o     =   0.46017;


    // GeV -> MeV conversion
    m2o     *=  1e6;
    m2r     *=  1e6;
    L2      *=  1e6;
    m2p     *=  1e6;
    Q2o     *=  1e6;

    // these values are corrected according to the file f2allm.f from Halina Abramowicz
    bp1     *=  bp1;
    bp2     *=  bp2;
    br1     *=  br1;
    br2     *=  br2;
    Q2o     +=  L2;

    const double R  =   0;

    double cr, ar, cp, ap, br, bp, t;

    t   =   log(log((Q2 + Q2o)/L2)/log(Q2o/L2));

    if(t<0)
    {
        t=0;
    }

    cr  =   cr1 + cr2*pow(t, cr3);
    ar  =   ar1 + ar2*pow(t, ar3);
    cp  =   cp1 + (cp1 - cp2)*(1/(1 + pow(t, cp3)) - 1);
    ap  =   ap1 + (ap1 - ap2)*(1/(1 + pow(t, ap3)) - 1);
    br  =   br1 + br2*pow(t, br3);
    bp  =   bp1 + bp2*pow(t, bp3);

    double xp, xr, F2p, F2r;

    xp  =   (Q2 + m2p)/(Q2 + m2p + W2 - component->GetAverageNucleonWeight()*component->GetAverageNucleonWeight());
    xr  =   (Q2 + m2r)/(Q2 + m2r + W2 - component->GetAverageNucleonWeight()*component->GetAverageNucleonWeight());
    F2p =   cp*pow(xp, ap)*pow(1 - x, bp);
    F2r =   cr*pow(xr, ar)*pow(1 - x, br);
    F2  =   (Q2/(Q2 + m2o))*(F2p + F2r)*G;
    R2  =   (2*(1 + R));

    aux =   ME*RE/Q2;
    aux *=  aux*(1 - v_ - component->GetAverageNucleonWeight()*x*v_/(2*particle.GetEnergy()) +
                 (1 - 2*particle.GetMass()*particle.GetMass()/Q2)
                 *v_*v_*(1 + 4*component->GetAverageNucleonWeight()*component->GetAverageNucleonWeight()*x*x/Q2)/R2);

    return (4*PI*F2/v_)*aux;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



double Photonuclear::FunctionToIntegralButMik(const PROPOSALParticle& particle, double Q2)
{
    Components::Component* component = medium_->GetComponents().at(component_);

    double x, aux, nu, G, F2, R2;

    nu  =   v_*particle.GetEnergy();
    x   =   Q2/(2*component->GetAverageNucleonWeight()*nu);

    G = ShadowEffect(x , nu);

    const double a  =   0.2513e6;
    const double b  =   0.6186e6;
    const double c  =   3.0292e6;
    const double d  =   1.4817e6;
    const double d0 =   0.0988;
    const double ar =   0.4056;
    const double t  =   1.8152;
    const double As =   0.12;
    const double Bu =   1.2437;
    const double Bd =   0.1853;
    const double R  =   0.25;

    double F2p, F2n, FSp, FNp, FSn, FNn, n, dl, xUv, xDv;

    n   =   1.5*(1 + Q2/(Q2 + c));
    dl  =   d0*(1 + 2*Q2/(Q2 + d));
    aux =   As*pow(x, -dl)*pow(Q2/(Q2 + a), 1 + dl);
    FSp =   aux*pow(1 - x, n + 4);
    FSn =   aux*pow(1 - x, n + t);
    aux =   pow(x, 1 - ar)*pow(1 - x, n)*pow(Q2/(Q2 + b), ar);
    xUv =   Bu*aux;
    xDv =   Bd*aux*(1 - x);
    FNp =   xUv + xDv;
    FNn =   xUv/4 + xDv*4;
    F2p =   FSp + FNp;
    F2n =   FSn + FNn;
    F2  =   G*(component->GetNucCharge()*F2p + (component->GetAtomicNum() - component->GetNucCharge())*F2n);
    R2  =   (2*(1 + R));


    aux =   ME*RE/Q2;
    aux *=  aux*(1 - v_ - component->GetAverageNucleonWeight()*x*v_/(2*particle.GetEnergy()) +
                 (1 - 2*particle.GetMass()*particle.GetMass()/Q2)
                 *v_*v_*(1 + 4*component->GetAverageNucleonWeight()*component->GetAverageNucleonWeight()*x*x/Q2)/R2);

    return (4*PI*F2/v_)*aux;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::FunctionToIntegralRSS(const PROPOSALParticle& particle, double Q2)
{
    Components::Component* component = medium_->GetComponents().at(component_);

    double x, aux, nu;

    nu  =   v_*particle.GetEnergy();
    x   =   Q2/(2*component->GetAverageNucleonWeight()*nu);

    // -------------[ Evaluate shadowfactor ]---------------- //

    double a;

    if(component->GetNucCharge()==1)
    {
        a   =   1;
    }
    else
    {
        a = ShadowEffect(x , nu);
    }

    double P;

    aux =   x*x;
    P   =   1 - 1.85*x + 2.45*aux - 2.35*aux*x + aux*aux;

    // ---------[ Evaluate ALLM form factor F_2 ]------------ //
    //
    // F_2 = c_i(t) * x_i^{a_i(t)} * (1 - x)^{b_i(t)}; i = P,R
    // ------------------------------------------------------ //

    double cp1, cp2, cp3;
    double cr1, cr2, cr3;

    double ap1, ap2, ap3;
    double ar1, ar2, ar3;

    double bp1, bp2, bp3;
    double br1, br2, br3;

    double m2o, m2r, L2, m2p, Q2o;

    cp1     =   0.28067;
    cp2     =   0.22291;
    cp3     =   2.1979;

    cr1     =   0.80107;
    cr2     =   0.97307;
    cr3     =   3.4942;

    ap1     =   -0.0808;
    ap2     =   -0.44812;
    ap3     =   1.1709;

    ar1     =   0.58400;
    ar2     =   0.37888;
    ar3     =   2.6063;

    bp1     =   0.60243;
    bp2     =   1.3754;
    bp3     =   1.8439;

    br1     =   0.10711;
    br2     =   1.9386;
    br3     =   0.49338;

    m2o     =   0.31985;
    m2r     =   0.15052;
    L2      =   0.06527;
    m2p     =   49.457;
    Q2o     =   0.46017;


    // GeV -> MeV conversion
    m2o     *=  1e6;
    m2r     *=  1e6;
    L2      *=  1e6;
    m2p     *=  1e6;
    Q2o     *=  1e6;

    // these values are corrected according to the file f2allm.f from Halina Abramowicz
    bp1     *=  bp1;
    bp2     *=  bp2;
    br1     *=  br1;
    br2     *=  br2;
    Q2o     +=  L2;

    // R(x, Q^2) is approximated to 0
    const double R = 0;

    double cr, ar, cp, ap, br, bp, t;

    t   =   log(log((Q2 + Q2o)/L2)/log(Q2o/L2));

    if(t<0)
    {
        t=0;
    }

    cr  =   cr1 + cr2*pow(t, cr3);
    ar  =   ar1 + ar2*pow(t, ar3);
    cp  =   cp1 + (cp1 - cp2)*(1/(1 + pow(t, cp3)) - 1);
    ap  =   ap1 + (ap1 - ap2)*(1/(1 + pow(t, ap3)) - 1);
    br  =   br1 + br2*pow(t, br3);
    bp  =   bp1 + bp2*pow(t, bp3);

    double xp, xr, F2p, F2A, F2P, F2R, W2;

    W2  =   component->GetAverageNucleonWeight()*component->GetAverageNucleonWeight()
            - Q2 + 2*component->GetAverageNucleonWeight()*particle.GetEnergy()*v_;
    xp  =   (Q2 + m2p)/(Q2 + m2p + W2 - component->GetAverageNucleonWeight()*component->GetAverageNucleonWeight());
    xr  =   (Q2 + m2r)/(Q2 + m2r + W2 - component->GetAverageNucleonWeight()*component->GetAverageNucleonWeight());
    F2P =   cp*pow(xp, ap)*pow(1 - x, bp);
    F2R =   cr*pow(xr, ar)*pow(1 - x, br);
    F2p =   (Q2/(Q2 + m2o))*(F2P + F2R);
    F2A =   a*(component->GetNucCharge() + (component->GetAtomicNum() - component->GetNucCharge())) *P*F2p;

    // ---------[ Write together cross section ]------------- //

    aux =   ME*RE/Q2;
    aux *=  aux*(1 - v_ + 0.25*v_*v_ -
                (1 + 4*particle.GetMass()*particle.GetMass()/Q2)*0.25*v_*v_ *
                (1 + 4*component->GetAverageNucleonWeight()*component->GetAverageNucleonWeight()*x*x/Q2) /
                (1 + R));
                 // *v_*v_*(1 + 4*component->GetAverageNucleonWeight()*component->GetAverageNucleonWeight()*x*x/Q2)/R2);

    return (4*PI*F2A/v_)*aux;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// void Photonuclear::SetParametrization(ParametrizationType::Enum parametrization)
// {
//     parametrization_    =   parametrization;
//
//     switch(parametrization_)
//     {
//         case ParametrizationType::PhotoKokoulinShadowBezrukovSoft:
//             hard_component_         =   false;
//             break;
//         case ParametrizationType::PhotoKokoulinShadowBezrukovHard:
//             hard_component_         =   true;
//             break;
//         case ParametrizationType::PhotoRhodeShadowBezrukovSoft:
//             hard_component_         =   false;
//             break;
//         case ParametrizationType::PhotoRhodeShadowBezrukovHard:
//             hard_component_         =   true;
//             break;
//         case ParametrizationType::PhotoBezrukovBugaevShadowBezrukovSoft:
//             hard_component_         =   false;
//             break;
//         case ParametrizationType::PhotoBezrukovBugaevShadowBezrukovHard:
//             hard_component_         =   true;
//             break;
//         case ParametrizationType::PhotoZeusShadowBezrukovSoft:
//             hard_component_         =   false;
//             break;
//         case ParametrizationType::PhotoZeusShadowBezrukovHard:
//             hard_component_         =   true;
//             break;
//         case ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowDutta:
//             shadow_                 =   ShadowingType::Dutta;
//             break;
//         case ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowButkevich:
//             shadow_                 =   ShadowingType::ButkevichMikhailov;
//             break;
//         case ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowDutta:
//             shadow_                 =   ShadowingType::Dutta;
//             break;
//         case ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich:
//             shadow_                 =   ShadowingType::ButkevichMikhailov;
//             break;
//         case ParametrizationType::PhotoButkevichMikhailovShadowDutta:
//             shadow_                 =   ShadowingType::Dutta;
//             break;
//         case ParametrizationType::PhotoButkevichMikhailovShadowButkevich:
//             shadow_                 =   ShadowingType::ButkevichMikhailov;
//             break;
//         case ParametrizationType::PhotoRenoSarcevicSuShadowDutta:
//             shadow_                 =   ShadowingType::Dutta;
//             break;
//         case ParametrizationType::PhotoRenoSarcevicSuShadowButkevich:
//             shadow_                 =   ShadowingType::ButkevichMikhailov;
//             break;
//         default:
//             log_warn("Photonuclear Parametrization '%i' not supported. Set to default '%i' "
//                 ,parametrization_, ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich);
//             parametrization_ =  ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich;
//             shadow_          =  ShadowingType::ButkevichMikhailov;
//     }
//
//     if(do_photo_interpolation_)
//     {
//         log_warn("photo-interpolation enabled before choosing the parametrization. Building the tables again");
//         DisablePhotoInterpolation();
//         EnablePhotoInterpolation();
//     }
//
//     if(do_dedx_Interpolation_)
//     {
//         log_warn("dEdx-interpolation enabled before choosing the parametrization. Building the tables again");
//         DisableDEdxInterpolation();
//         EnableDEdxInterpolation();
//     }
//     if(do_dndx_Interpolation_)
//     {
//         log_warn("dNdx-interpolation enabled before choosing the parametrization. Building the tables again");
//         DisableDNdxInterpolation();
//         EnableDNdxInterpolation();
//     }
// }



