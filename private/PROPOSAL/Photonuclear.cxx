#include "PROPOSAL/Photonuclear.h"
#include <algorithm>
#include "PROPOSAL/Output.h"

using namespace std;

namespace po	= boost::program_options;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::CalculatedEdx()
{
    if(multiplier_<=0)
    {
        return 0;
    }

    if(do_dedx_Interpolation_)
    {
        return max(dedx_interpolant_->Interpolate(particle_->GetEnergy()), 0.0);
    }

    double sum  =   0;

    for(int i=0; i < medium_->GetNumComponents(); i++)
    {
        SetIntegralLimits(i);
        sum +=  integral_for_dEdx_->Integrate(vMin_, vUp_, boost::bind(&Photonuclear::FunctionToDEdxIntegral, this, _1),4);
    }

    return multiplier_*particle_->GetEnergy()*sum;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::CalculatedNdx()
{
    if(multiplier_<=0)
    {
        return 0;
    }

    sum_of_rates_    =   0;

    for(int i=0; i<medium_->GetNumComponents(); i++)
    {

        if(do_dndx_Interpolation_)
        {
            sum_of_rates_    +=  max(dndx_interpolant_1d_.at(i)->Interpolate(particle_->GetEnergy()), 0.0);
        }
        else
        {
            SetIntegralLimits(i);
            sum_of_rates_    +=  dndx_integral_.at(i)->Integrate(vUp_, vMax_, boost::bind(&Photonuclear::FunctionToDNdxIntegral, this, _1),4);
        }
    }

    return sum_of_rates_;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::CalculatedNdx(double rnd)
{
    if(multiplier_<=0)
    {
        return 0.;
    }

    // The random number will be stored to be able
    // to check if dNdx is already calculated for this random number.
    // This avoids a second calculation in CalculateStochaticLoss

    rnd_    =   rnd;

    sum_of_rates_   =   0;

    for(int i=0; i<medium_->GetNumComponents(); i++)
    {

        if(do_dndx_Interpolation_)
        {
            prob_for_component_.at(i) = max(dndx_interpolant_1d_.at(i)->Interpolate(particle_->GetEnergy()), 0.0);
        }
        else
        {
            SetIntegralLimits(i);
            prob_for_component_.at(i) = dndx_integral_.at(i)->IntegrateWithLog(vUp_,vMax_, boost::bind(&Photonuclear::FunctionToDNdxIntegral, this, _1),rnd);
        }

        sum_of_rates_    +=  prob_for_component_.at(i);
    }

    return sum_of_rates_;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::CalculateStochasticLoss(double rnd1, double rnd2)
{
    if(rnd1 != rnd_ )
    {
        CalculatedNdx(rnd1);
    }

    return CalculateStochasticLoss(rnd2);

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-----------------------Enable and disable interpolation---------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Photonuclear::EnableDNdxInterpolation(std::string path, bool raw)
{
    if(do_dndx_Interpolation_)return;

    EnablePhotoInterpolation(path,raw);

    bool storing_failed =   false;
    bool reading_worked =   true;

    if(!path.empty())
    {
        stringstream filename;
        filename<<path<<"/Photo_dNdx_particle_"<<particle_->GetName()
                <<"_para_"<<parametrization_
                <<"_med_"<<medium_->GetName()
                <<medium_->GetMassDensity()
                <<"_ecut_"<<cut_settings_->GetEcut()
                <<"_vcut_"<<cut_settings_->GetVcut()
                <<"_multiplier_"<<multiplier_;

        if(!raw)
            filename<<".txt";

        dndx_interpolant_2d_.resize(medium_->GetNumComponents());
        dndx_interpolant_1d_.resize(medium_->GetNumComponents());

        if( FileExist(filename.str()) )
        {
            log_info("Photonuclear parametrisation tables (dNdx) will be read from file:\t%s",filename.str().c_str());
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

            double energy = particle_->GetEnergy();

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

                    dndx_interpolant_2d_.at(i) =    new Interpolant(NUM1, particle_->GetLow(), BIGENERGY,  NUM1, 0, 1, boost::bind(&Photonuclear::FunctionToBuildDNdxInterpolant2D, this, _1 , _2), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false, order_of_interpolation_, true, false, false);
                    dndx_interpolant_1d_.at(i) =    new Interpolant(NUM1, particle_->GetLow(), BIGENERGY,  boost::bind(&Photonuclear::FunctionToBuildDNdxInterpolant1D, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, true, false, false);

                    dndx_interpolant_2d_.at(i)->Save(output, raw);
                    dndx_interpolant_1d_.at(i)->Save(output, raw);

                }
            }
            else
            {
                storing_failed  =   true;
                log_warn("Can not open file %s for writing! Table will not be stored!",filename.str().c_str());
            }
            particle_->SetEnergy(energy);

            output.close();
        }
    }
    if(path.empty() || storing_failed)
    {
        EnablePhotoInterpolation(path);

        double energy = particle_->GetEnergy();
        dndx_interpolant_1d_.resize(medium_->GetNumComponents());
        dndx_interpolant_2d_.resize(medium_->GetNumComponents());
        for(int i=0; i<(medium_->GetNumComponents()); i++)
        {
            component_ = i;
            dndx_interpolant_2d_.at(i) =    new Interpolant(NUM1, particle_->GetLow(), BIGENERGY,  NUM1, 0, 1, boost::bind(&Photonuclear::FunctionToBuildDNdxInterpolant2D, this, _1 , _2), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false, order_of_interpolation_, true, false, false);
            dndx_interpolant_1d_.at(i) =    new Interpolant(NUM1, particle_->GetLow(), BIGENERGY,  boost::bind(&Photonuclear::FunctionToBuildDNdxInterpolant1D, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, true, false, false);
        }
        particle_->SetEnergy(energy);
    }

    do_dndx_Interpolation_=true;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Photonuclear::EnableDEdxInterpolation(std::string path, bool raw)
{
    if(do_dedx_Interpolation_)return;

    EnablePhotoInterpolation(path,raw);

    bool reading_worked =   true;
    bool storing_failed =   false;

    if(!path.empty())
    {
        stringstream filename;
        filename<<path<<"/Photo_dEdx_particle_"<<particle_->GetName()
                <<"_para_"<<parametrization_
                <<"_med_"<<medium_->GetName()
                <<medium_->GetMassDensity()
                <<"_ecut_"<<cut_settings_->GetEcut()
                <<"_vcut_"<<cut_settings_->GetVcut()
                <<"_multiplier_"<<multiplier_;

        if(!raw)
            filename<<".txt";

        if( FileExist(filename.str()) )
        {
            log_info("Photonuclear parametrisation tables (dEdx) will be read from file:\t%s",filename.str().c_str());
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

            double energy = particle_->GetEnergy();

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

                dedx_interpolant_ = new Interpolant(NUM1, particle_->GetLow(), BIGENERGY, boost::bind(&Photonuclear::FunctionToBuildDEdxInterpolant, this, _1),
                                                    order_of_interpolation_, true, false, true, order_of_interpolation_, false, false, false);
                dedx_interpolant_->Save(output, raw);
            }
            else
            {
                storing_failed  =   true;
                log_warn("Can not open file %s for writing! Table will not be stored!",filename.str().c_str());
            }
            particle_->SetEnergy(energy);

            output.close();
        }
    }
    if(path.empty() || storing_failed)
    {
        double energy = particle_->GetEnergy();
        dedx_interpolant_ = new Interpolant(NUM1, particle_->GetLow(), BIGENERGY, boost::bind(&Photonuclear::FunctionToBuildDEdxInterpolant, this, _1),
                                            order_of_interpolation_, true, false, true, order_of_interpolation_, false, false, false);
        particle_->SetEnergy(energy);

    }

    do_dedx_Interpolation_=true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Photonuclear::EnablePhotoInterpolation(std::string path, bool raw)
{
    if(do_photo_interpolation_)return;

    bool storing_failed =   false;
    bool reading_worked =   true;

    if(!path.empty())
    {
        stringstream filename;
        filename<<path<<"/Photo_particle_"<<particle_->GetName()
                <<"_para_"<<parametrization_
                <<"_med_"<<medium_->GetName()
                <<medium_->GetMassDensity()
                <<"_ecut_"<<cut_settings_->GetEcut()
                <<"_vcut_"<<cut_settings_->GetVcut()
                <<"_multiplier_"<<multiplier_;

        if(!raw)
            filename<<".txt";

        photo_interpolant_.resize( medium_->GetNumComponents() );

        if( FileExist(filename.str()) )
        {
            log_info("Photonuclear parametrisation tables will be read from file:\t%s",filename.str().c_str());
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

            double energy = particle_->GetEnergy();

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

                    photo_interpolant_.at(i)  = new Interpolant(NUM1, particle_->GetLow(), BIGENERGY, NUM1, 0., 1., boost::bind(&Photonuclear::FunctionToBuildPhotoInterpolant, this, _1, _2), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false, order_of_interpolation_, false, false, false);

                    photo_interpolant_.at(i)->Save(output, raw);

                }
            }
            else
            {
                storing_failed  =   true;
                log_warn("Can not open file %s for writing! Table will not be stored!",filename.str().c_str());
            }
            particle_->SetEnergy(energy);

            output.close();
        }
    }
    if(path.empty() || storing_failed)
    {
        double energy = particle_->GetEnergy();

        photo_interpolant_.resize( medium_->GetNumComponents() );

        for(int i=0; i<medium_->GetNumComponents(); i++)
        {
            component_ = i;
            photo_interpolant_.at(i)  = new Interpolant(NUM1, particle_->GetLow(), BIGENERGY, NUM1, 0., 1., boost::bind(&Photonuclear::FunctionToBuildPhotoInterpolant, this, _1, _2), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false, order_of_interpolation_, false, false, false);

        }
        particle_->SetEnergy(energy);
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


boost::program_options::options_description Photonuclear::CreateOptions()
{
    po::options_description photonuclear("Photonuclear options");
    photonuclear.add_options()
        ("photonuclear.para",             po::value<int>()->notifier(boost::bind(&Photonuclear::SetParametrization, this, _1)),              "1 = Kakoulin \n2 = Rhode \n3 = Bezrukov-Bugaev \n4 = Zeus \n5 = ALLM 91 \n6 = ALLM 97 \n7 = Butkevich-Mikhailov")
        ("photonuclear.interpol_dedx",    po::value<bool>(&do_dedx_Interpolation_)->implicit_value(false),  "Enables interpolation for dEdx")
        ("photonuclear.interpol_dndx",    po::value<bool>(&do_dndx_Interpolation_)->implicit_value(false),  "Enables interpolation for dNdx")
        ("photonuclear.multiplier",       po::value<double>(&multiplier_)->default_value(1.),               "modify the cross section by this factor")
        ("photonuclear.interpol_order",   po::value<int>(&order_of_interpolation_)->default_value(5),       "number of interpolation points");

   return photonuclear;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Photonuclear::ValidateOptions()
{
    if(order_of_interpolation_ < 2)
    {
        order_of_interpolation_ = 5;
        cerr<<"Photonuclear: Order of Interpolation is not a vaild number. Must be > 2.\t"<<"Set to 5"<<endl;
    }
    if(order_of_interpolation_ > 6)
    {
        cerr<<"Photonuclear: Order of Interpolation is set to "<<order_of_interpolation_
            <<".\t Note a order of interpolation > 6 will slow down the program"<<endl;
    }
    if(parametrization_ < 1 || parametrization_ >14)
    {
        SetParametrization(12);
        cerr<<"Photonuclear: Parametrization is not a vaild number. Must be 1-14.\tSet parametrization to 12 (icecube default)"<<endl;
    }

}


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
    ,shadow_                ( 1 )
    ,hard_component_        ( false )
    ,parametrization_family_( 1 )
    ,dndx_integral_         ( )
    ,interpolant_hardBB_    ( )
    ,dndx_interpolant_1d_   ( )
    ,dndx_interpolant_2d_   ( )
    ,photo_interpolant_     ( )
    ,prob_for_component_    ( )

{
    name_   =   "Photonuclear";

    integral_             = new Integral(IROMB, IMAXS, IPREC);
    integral_for_dEdx_    = new Integral(IROMB, IMAXS, IPREC);

    dndx_integral_.resize(medium_->GetNumComponents());

    for(int i =0 ; i<medium_->GetNumComponents();i++){
        dndx_integral_.at(i) = new Integral(IROMB, IMAXS, IPREC);
    }

    prob_for_component_.resize(medium_->GetNumComponents());
    dedx_interpolant_     = NULL;
    interpolant_measured_ = NULL;
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
    ,shadow_                            ( photo.shadow_ )
    ,hard_component_                    ( photo.hard_component_ )
    ,parametrization_family_            ( photo.parametrization_family_ )
    ,integral_                          ( new Integral(*photo.integral_) )
    ,integral_for_dEdx_                 ( new Integral(*photo.integral_for_dEdx_) )
    ,prob_for_component_                ( photo.prob_for_component_ )
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


Photonuclear::Photonuclear(Particle* particle,
                             Medium* medium,
                             EnergyCutSettings* cut_settings)
    :CrossSections(particle, medium, cut_settings)
    ,component_             ( 0 )
    ,init_measured_         ( true )
    ,init_hardbb_           ( true )
    ,hmax_                  ( 8 )
    ,v_                     ( 0 )
    ,do_photo_interpolation_( false )
    ,shadow_                ( 2 )
    ,hard_component_        ( false )
    ,parametrization_family_( 1 )
    ,dndx_integral_         ( )
    ,interpolant_hardBB_    ( )
    ,dndx_interpolant_1d_   ( )
    ,dndx_interpolant_2d_   ( )
    ,photo_interpolant_     ( )
    ,prob_for_component_    ( )


{
    name_       = "Photonuclear";
    multiplier_ = 1.;

    integral_             = new Integral(IROMB, IMAXS, IPREC);
    integral_for_dEdx_    = new Integral(IROMB, IMAXS, IPREC);

    dndx_integral_.resize(medium_->GetNumComponents());

    for(int i =0 ; i<medium_->GetNumComponents();i++){
        dndx_integral_.at(i) = new Integral(IROMB, IMAXS, IPREC);
    }

    prob_for_component_.resize(medium_->GetNumComponents());
    dedx_interpolant_     = NULL;
    interpolant_measured_ = NULL;


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
    if( shadow_                     !=  photo.shadow_ )                     return false;
    if( hard_component_             !=  photo.hard_component_ )             return false;
    if( parametrization_family_     !=  photo.parametrization_family_ )     return false;
    if( init_hardbb_                !=  photo.init_hardbb_)                 return false;
    if( init_measured_              !=  photo.init_measured_)               return false;
    if( *integral_                  != *photo.integral_)                    return false;
    if( *integral_for_dEdx_         != *photo.integral_for_dEdx_)           return false;
    if( prob_for_component_.size()  !=  photo.prob_for_component_.size())   return false;
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
    swap(shadow_ , photo.shadow_ );
    swap(hard_component_ , photo.hard_component_ );
    swap(parametrization_family_    ,photo.parametrization_family_ );

    integral_for_dEdx_->swap(*photo.integral_for_dEdx_);
    integral_->swap(*photo.integral_);

    prob_for_component_.swap(photo.prob_for_component_);
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

ostream& operator<<(std::ostream& os, Photonuclear const &photo)
{
    os<<"---------------------------Photonuclear( "<<&photo<<" )---------------------------"<<std::endl;
    os<< static_cast <const CrossSections &>( photo ) << endl;
    os<< "------- Class Specific: " << endl;
    os<< "\tinit_measured:\t\t" << photo.init_measured_ << endl;
    os<< "\tinit_hardbb:\t\t" << photo.init_hardbb_ << endl;
    os<< "\tshadow:\t\t\t" << photo.shadow_ << endl;
    os<< "\thard_component:\t\t" << photo.hard_component_ << endl;
    os<< "\tparametrization_family_:\t" << photo.parametrization_family_ << endl;
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
    os<<"-----------------------------------------------------------------------------------------------";
    return os;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------Parametrizations--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::KokoulinParametrization(double v, int i)
{
    double nu, sgn, aux, aum, k, G, t;

    const double m1 =   0.54;
    const double m2 =   1.80;
    nu              =   v*particle_->GetEnergy();
    nu              *=  1.e-3;

    if(nu<=200.)
    {
        if(nu<=17.)
        {
            sgn =   96.1+82./sqrt(nu);
        }
        else
        {
            aux =   log(0.0213*nu);
            sgn =   114.3 + 1.647*aux*aux;
        }
    }
    else
    {
        sgn =   49.2 + 11.1*log(nu) + 151.8/sqrt(nu);
    }


    k   =   1 - 2/v + 2/(v*v);

    if(medium_->GetNucCharge().at(i)==1)
    {
        G   =   1;
    }
    else
    {
        aux =   0.00282*pow(medium_->GetAtomicNum().at(i), 1./3)*sgn;
        G   =   (3/aux)*(0.5 + ((1 + aux)*exp(-aux) - 1)/(aux*aux));
    }

    G       *=  3;
    aux     =   v*particle_->GetMass()*1.e-3;
    t       =   aux*aux/(1-v);
    sgn     *=  1.e-30;
    aum     =   particle_->GetMass()*1.e-3;
    aum     *=  aum;
    aux     =   2*aum/t;
    aux     =   G*((k + 4*aum/m1)*log(1 + m1/t) - (k*m1)/(m1 + t) - aux) +
                ((k + 2*aum/m2)*log(1 + m2/t) - aux) + 0.25*aux*((G*(m1 - 4*t))/(m1 + t) + (m2/t)*log(1 + t/m2));

    aux     *=  ALPHA/(8*PI)*medium_->GetAtomicNum().at(i)*sgn*v;

    if(hard_component_)
    {
        if(particle_->GetType()==1 || particle_->GetType()==2)
        {
            aux +=  medium_->GetAtomicNum().at(i)*1.e-30*HardBB(particle_->GetEnergy(), v);
        }

    }

    return medium_->GetMolDensity()*medium_->GetAtomInMolecule().at(i)*particle_->GetCharge()*particle_->GetCharge()*aux;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::RhodeParametrization(double v, int i)
{
    double nu, sgn, aux, aum, k, G, t;

    const double m1 =   0.54;
    const double m2 =   1.80;
    nu              =   v*particle_->GetEnergy();
    nu              *=  1.e-3;

    if(nu<=200.)
    {
        sgn =   MeasuredSgN(nu);
    }
    else
    {
        sgn =   49.2 + 11.1*log(nu) + 151.8/sqrt(nu);
    }


    k   =   1 - 2/v + 2/(v*v);

    if(medium_->GetNucCharge().at(i)==1)
    {
        G   =   1;
    }
    else
    {
        aux =   0.00282*pow(medium_->GetAtomicNum().at(i), 1./3)*sgn;
        G   =   (3/aux)*(0.5 + ((1 + aux)*exp(-aux) - 1)/(aux*aux));
    }

    G       *=  3;
    aux     =   v*particle_->GetMass()*1.e-3;
    t       =   aux*aux/(1-v);
    sgn     *=  1.e-30;
    aum     =   particle_->GetMass()*1.e-3;
    aum     *=  aum;
    aux     =   2*aum/t;
    aux     =   G*((k + 4*aum/m1)*log(1 + m1/t) - (k*m1)/(m1 + t) - aux) +
                ((k + 2*aum/m2)*log(1 + m2/t) - aux) + 0.25*aux*((G*(m1 - 4*t))/(m1 + t) + (m2/t)*log(1 + t/m2));

    aux     *=  ALPHA/(8*PI)*medium_->GetAtomicNum().at(i)*sgn*v;

    if(hard_component_)
    {
        if(particle_->GetType()==1 || particle_->GetType()==2)
        {
            aux +=  medium_->GetAtomicNum().at(i)*1.e-30*HardBB(particle_->GetEnergy(), v);
        }

    }

    return medium_->GetMolDensity()*medium_->GetAtomInMolecule().at(i)*particle_->GetCharge()*particle_->GetCharge()*aux;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::BezrukovBugaevParametrization(double v, int i)
{
    double nu, sgn, aux, aum, k, G, t;

    const double m1 =   0.54;
    const double m2 =   1.80;
    nu              =   v*particle_->GetEnergy();
    nu              *=  1.e-3;


    aux =   log(0.0213*nu);
    sgn =   114.3 + 1.647*aux*aux;

    k   =   1 - 2/v + 2/(v*v);

    if(medium_->GetNucCharge().at(i)==1)
    {
        G   =   1;
    }
    else
    {
        aux =   0.00282*pow(medium_->GetAtomicNum().at(i), 1./3)*sgn;
        G   =   (3/aux)*(0.5 + ((1 + aux)*exp(-aux) - 1)/(aux*aux));
    }

    G       *=  3;
    aux     =   v*particle_->GetMass()*1.e-3;
    t       =   aux*aux/(1-v);
    sgn     *=  1.e-30;
    aum     =   particle_->GetMass()*1.e-3;
    aum     *=  aum;
    aux     =   2*aum/t;
    aux     =   G*((k + 4*aum/m1)*log(1 + m1/t) - (k*m1)/(m1 + t) - aux) +
                ((k + 2*aum/m2)*log(1 + m2/t) - aux) + 0.25*aux*((G*(m1 - 4*t))/(m1 + t) + (m2/t)*log(1 + t/m2));

    aux     *=  ALPHA/(8*PI)*medium_->GetAtomicNum().at(i)*sgn*v;

    if(hard_component_)
    {
        if(particle_->GetType()==1 || particle_->GetType()==2)
        {
            aux +=  medium_->GetAtomicNum().at(i)*1.e-30*HardBB(particle_->GetEnergy(), v);
        }

    }

    return medium_->GetMolDensity()*medium_->GetAtomInMolecule().at(i)*particle_->GetCharge()*particle_->GetCharge()*aux;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::ZeusParametrization(double v, int i)
{
    double nu, sgn, aux, aum, k, G, t;

    const double m1 =   0.54;
    const double m2 =   1.80;
    nu              =   v*particle_->GetEnergy();
    nu              *=  1.e-3;

    aux =   nu*2.e-3*medium_->GetAverageNucleonWeight().at(i);
    sgn =   63.5*pow(aux, 0.097) + 145*pow(aux , -0.5);


    k   =   1 - 2/v + 2/(v*v);

    if(medium_->GetNucCharge().at(i)==1)
    {
        G   =   1;
    }
    else
    {
        aux =   0.00282*pow(medium_->GetAtomicNum().at(i), 1./3)*sgn;
        G   =   (3/aux)*(0.5 + ((1 + aux)*exp(-aux) - 1)/(aux*aux));
    }

    G       *=  3;
    aux     =   v*particle_->GetMass()*1.e-3;
    t       =   aux*aux/(1-v);
    sgn     *=  1.e-30;
    aum     =   particle_->GetMass()*1.e-3;
    aum     *=  aum;
    aux     =   2*aum/t;
    aux     =   G*((k + 4*aum/m1)*log(1 + m1/t) - (k*m1)/(m1 + t) - aux) +
                ((k + 2*aum/m2)*log(1 + m2/t) - aux) + 0.25*aux*((G*(m1 - 4*t))/(m1 + t) + (m2/t)*log(1 + t/m2));

    aux     *=  ALPHA/(8*PI)*medium_->GetAtomicNum().at(i)*sgn*v;

    if(hard_component_)
    {
        if(particle_->GetType()==1 || particle_->GetType()==2)
        {
            aux +=  medium_->GetAtomicNum().at(i)*1.e-30*HardBB(particle_->GetEnergy(), v);
        }

    }

    return medium_->GetMolDensity()*medium_->GetAtomInMolecule().at(i)*particle_->GetCharge()*particle_->GetCharge()*aux;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::ALLM91Parametrization(double v, int i)
{
    if(do_photo_interpolation_)
    {
        SetIntegralLimits(i);
        if(v >= vUp_)
        {
            return max(photo_interpolant_.at(i)->Interpolate(particle_->GetEnergy(), log(v/vUp_)/log(vMax_/vUp_)), 0.0);
        }
    }

    double aux, min, max;

    component_ =   i;
    v_         =   v;
    min        =   particle_->GetMass()*v;
    min        *=  min/(1-v);

    if(particle_->GetMass() < MPI)
    {
        aux     =   particle_->GetMass()*particle_->GetMass()/particle_->GetEnergy();
        min     -=  (aux*aux)/(2*(1-v));
    }

    max =   2*medium_->GetAverageNucleonWeight().at(i)*particle_->GetEnergy()*(v-vMin_);

    //  if(form==4) max=Math.min(max, 5.5e6);  // as requested in Butkevich and Mikheyev
    if(min > max)
    {
        return 0;
    }

    return medium_->GetMolDensity()*medium_->GetAtomInMolecule().at(i)*particle_->GetCharge()*particle_->GetCharge()*integral_->Integrate(min, max, boost::bind(&Photonuclear::FunctionToIntegralALLM91, this, _1),4);

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::ALLM97Parametrization(double v, int i)
{
    if(do_photo_interpolation_)
    {
        SetIntegralLimits(i);
        if(v >= vUp_)
        {
            return max(photo_interpolant_.at(i)->Interpolate(particle_->GetEnergy(), log(v/vUp_)/log(vMax_/vUp_)), 0.0);
        }
    }

    double aux, min, max;

    component_ =   i;
    v_         =   v;
    min        =   particle_->GetMass()*v;
    min        *=  min/(1-v);

    if(particle_->GetMass() < MPI)
    {
        aux     =   particle_->GetMass()*particle_->GetMass()/particle_->GetEnergy();
        min     -=  (aux*aux)/(2*(1-v));
    }

    max =   2*medium_->GetAverageNucleonWeight().at(i)*particle_->GetEnergy()*(v-vMin_);

    //  if(form==4) max=Math.min(max, 5.5e6);  // as requested in Butkevich and Mikheyev
    if(min>max)
    {
        return 0;
    }

    return medium_->GetMolDensity()*medium_->GetAtomInMolecule().at(i)*particle_->GetCharge()*particle_->GetCharge()*integral_->Integrate(min, max, boost::bind(&Photonuclear::FunctionToIntegralALLM97, this, _1),4);

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::ButkevichMikhailovParametrization(double v, int i)
{
    if(do_photo_interpolation_)
    {
        SetIntegralLimits(i);
        if(v >= vUp_)
        {
            return max(photo_interpolant_.at(i)->Interpolate(particle_->GetEnergy(), log(v/vUp_)/log(vMax_/vUp_)), 0.0);
        }
    }

    double aux, min, max;

    component_ =   i;
    v_         =   v;
    min        =   particle_->GetMass()*v;
    min        *=  min/(1-v);

    if(particle_->GetMass() < MPI)
    {
        aux     =   particle_->GetMass()*particle_->GetMass()/particle_->GetEnergy();
        min     -=  (aux*aux)/(2*(1-v));
    }

    max =   2*medium_->GetAverageNucleonWeight().at(i)*particle_->GetEnergy()*(v-vMin_);

    //  if(form==4) max=Math.min(max, 5.5e6);  // as requested in Butkevich and Mikheyev
    if(min>max)
    {
        return 0;
    }

    return medium_->GetMolDensity()*medium_->GetAtomInMolecule().at(i)*particle_->GetCharge()*particle_->GetCharge()*integral_->Integrate(min, max, boost::bind(&Photonuclear::FunctionToIntegralButMik, this, _1),4);

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------Cross section / limit / private-----------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::CalculateStochasticLoss(double rnd)
{
    double rand;
    double rsum;

    rand    =   rnd*sum_of_rates_;
    rsum    =   0;


    for(int i=0; i<(medium_->GetNumComponents()); i++)
    {
        rsum    += prob_for_component_.at(i);
        if(rsum > rand)
        {
            if(do_dndx_Interpolation_)
            {
                SetIntegralLimits(i);

                if(vUp_==vMax_)
                {
                    return (particle_->GetEnergy())*vUp_;
                }

                return particle_->GetEnergy()*(vUp_*exp(dndx_interpolant_2d_.at(i)->FindLimit((particle_->GetEnergy()), (rnd_)*prob_for_component_.at(i))*log(vMax_/vUp_)));

            }
            else
            {
                component_ = i;
                return (particle_->GetEnergy())*dndx_integral_.at(i)->GetUpperLimit();

            }
        }
    }

    //sometime everything is fine, just the probability for interaction is zero
    bool prob_for_all_comp_is_zero=true;
    for(int i=0; i<(medium_->GetNumComponents()); i++)
    {
        SetIntegralLimits(i);
        if(vUp_!=vMax_)prob_for_all_comp_is_zero=false;
    }
    if(prob_for_all_comp_is_zero)return 0;

    log_fatal("sum was not initialized correctly");

    return 0;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Photonuclear::SetIntegralLimits(int component)
{

    double aux;

    component_ = component;
    vMin_    =   (MPI + (MPI*MPI)/(2*medium_->GetAverageNucleonWeight().at(component_)))/particle_->GetEnergy();

    if(particle_->GetMass() < MPI)
    {
        aux     =   particle_->GetMass()/medium_->GetAverageNucleonWeight().at(component_);
        vMax_    =   1 - medium_->GetAverageNucleonWeight().at(component_)*(1 + aux*aux)/(2*particle_->GetEnergy());
    }
    else
    {
        vMax_    =   1;
    }

    vMax_    =   min(vMax_, 1-particle_->GetMass()/particle_->GetEnergy());

    if(vMax_ < vMin_)
    {
        vMax_    =   vMin_;
    }

    vUp_     =   min(vMax_, cut_settings_->GetCut(particle_->GetEnergy()));

    if(vUp_ < vMin_)
    {
        vUp_ =   vMin_;
    }

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::PhotoN(double v, int i)
{
    switch(parametrization_family_)
    {
        case 1: return KokoulinParametrization(v, i);

        case 2: return RhodeParametrization(v, i);

        case 3: return BezrukovBugaevParametrization(v, i);

        case 4: return ZeusParametrization(v, i);

        case 5: return ALLM91Parametrization(v, i);

        case 6: return ALLM97Parametrization(v, i);

        case 7: return ButkevichMikhailovParametrization(v, i);

        default:
            log_error("parametrization_family_ %i is not supported! Be careful 0 is returned",parametrization_family_);
            return 0;

    }

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
//----------------------------------------------------------------------------//


void Photonuclear::EnableHardBB()
{
    if(init_hardbb_)
    {
        init_hardbb_           =   false;
        double x_aux[]  =   {3, 4, 5, 6, 7, 8, 9};


        double y_aux[][56][7]=
        {
            {
                {7.174409e-4, 1.7132e-3, 4.082304e-3, 8.628455e-3, 0.01244159, 0.02204591, 0.03228755},
                {-0.2436045, -0.5756682, -1.553973, -3.251305, -5.976818, -9.495636, -13.92918},
                {-0.2942209, -0.68615, -2.004218, -3.999623, -6.855045, -10.05705, -14.37232},
                {-0.1658391, -0.3825223, -1.207777, -2.33175, -3.88775, -5.636636, -8.418409},
                {-0.05227727, -0.1196482, -0.4033373, -0.7614046, -1.270677, -1.883845, -2.948277},
                {-9.328318e-3, -0.02124577, -0.07555636, -0.1402496, -0.2370768, -0.3614146, -0.5819409},
                {-8.751909e-4, -1.987841e-3, -7.399682e-3, -0.01354059, -0.02325118, -0.03629659, -0.059275},
                {-3.343145e-5, -7.584046e-5, -2.943396e-4, -5.3155e-4, -9.265136e-4, -1.473118e-3, -2.419946e-3}
            },
            {
                {-1.269205e-4, -2.843877e-4, -5.761546e-4, -1.195445e-3, -1.317386e-3, -9.689228e-15, -6.4595e-15},
                {-0.01563032, -0.03589573, -0.07768545, -0.157375, -0.2720009, -0.4186136, -0.8045046},
                {0.04693954, 0.1162945, 0.3064255, 0.7041273, 1.440518, 2.533355, 3.217832},
                {0.05338546, 0.130975, 0.3410341, 0.7529364, 1.425927, 2.284968, 2.5487},
                {0.02240132, 0.05496, 0.144945, 0.3119032, 0.5576727, 0.8360727, 0.8085682},
                {4.658909e-3, 0.01146659, 0.03090286, 0.06514455, 0.1109868, 0.1589677, 0.1344223},
                {4.822364e-4, 1.193018e-3, 3.302773e-3, 6.843364e-3, 0.011191, 0.015614, 0.01173827},
                {1.9837e-5, 4.940182e-5, 1.409573e-4, 2.877909e-4, 4.544877e-4, 6.280818e-4, 4.281932e-4}
            }
        };

        interpolant_hardBB_.resize(hmax_);

        for(int i=0; i < hmax_; i++)
        {
            vector<double> x(x_aux, x_aux + sizeof(x_aux) / sizeof(double) );
            vector<double> y(y_aux[particle_->GetType()-1][i], y_aux[particle_->GetType()-1][i] + sizeof(y_aux[particle_->GetType()-1][i]) / sizeof(double) );

            interpolant_hardBB_.at(i)    =   new Interpolant(x, y, 4, false, false) ;
        }
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::HardBB(double e, double v)
{
    EnableHardBB();

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
//----------------------------------------------------------------------------//


double Photonuclear::ShadowEffect(double x , double nu)
{
    double G , aux;

    if(shadow_ == 1)
    {
        if(x<0.0014)
        {
            G   =   pow(medium_->GetAtomicNum().at(component_), -0.1);
        }
        else if(x<0.04)
        {
            G   =   pow(medium_->GetAtomicNum().at(component_), 0.069*log(x)/LOG10+0.097);
        }
        else
        {
            G   =   1;
        }
    }
    else if (shadow_ == 2)
    {
        if(x>0.3)
        {
            const double Mb =   0.437;
            const double la =   0.5;
            const double x2 =   0.278;

            double mb, Aosc, mu, au, ac;

            mb      =   Mb*medium_->GetMN().at(component_);
            au      =   1/(1 - x);
            ac      =   1/(1 - x2);
            mu      =   MPI/medium_->GetAverageNucleonWeight().at(component_);
            Aosc    =   (1 - la*x)*((au - ac)-mu*(au*au - ac*ac));
            G       =   1 - mb*Aosc;
        }
        else
        {
            const double M1 =   0.129;
            const double M2 =   0.456;
            const double M3 =   0.553;

            double m1, m2, m3, x0, sgn;

            m1  =   M1*medium_->GetMN().at(component_);
            m2  =   M2*medium_->GetMN().at(component_);
            m3  =   M3*medium_->GetMN().at(component_);
            nu  *=  1.e-3;
            sgn =   112.2*(0.609*pow(nu, 0.0988) + 1.037*pow(nu, -0.5944));
            aux =   0.00282*pow(medium_->GetAtomicNum().at(component_), 1./3)*sgn;
            G   =   (3/aux)*(0.5 + ((1 + aux)*exp(-aux)-1)/(aux*aux));
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
        log_warn("shadow_ must be 1 or 2, other values are not supported! Be careful shadow effect factor is set to 1!");

        G   =   1.;
    }

    return G;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------Functions to interpolate---------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::FunctionToBuildDEdxInterpolant(double energy)
{
    particle_->SetEnergy(energy);
    return CalculatedEdx();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::FunctionToBuildDNdxInterpolant1D(double energy)
{
    return dndx_interpolant_2d_.at(component_)->Interpolate(energy, 1.);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::FunctionToBuildDNdxInterpolant2D(double energy, double v)
{

    particle_->SetEnergy(energy);
    SetIntegralLimits(component_);
    if(vUp_==vMax_)
    {
        return 0;
    }

    v   =   vUp_*exp(v*log(vMax_/vUp_));

    return dndx_integral_.at(component_)->Integrate(vUp_, v, boost::bind(&Photonuclear::FunctionToDNdxIntegral, this, _1) ,4);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::FunctionToBuildPhotoInterpolant(double energy, double v)
{
    particle_->SetEnergy(energy);
    SetIntegralLimits(component_);

    if(vUp_==vMax_)
    {
        return 0;
    }

    v   =   vUp_*exp(v*log(vMax_/vUp_));

    return PhotoN(v, component_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------Functions to integrate----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::FunctionToDEdxIntegral(double variable)
{
    return variable*PhotoN(variable, component_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::FunctionToDNdxIntegral(double variable)
{
    return multiplier_*PhotoN(variable, component_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::FunctionToIntegralALLM91(double Q2)
{
    double x, aux, nu, G, F2, R2;

    nu  =   v_*particle_->GetEnergy();
    x   =   Q2/(2*medium_->GetAverageNucleonWeight().at(component_)*nu);

    if(medium_->GetNucCharge().at(component_)==1)
    {
        G   =   1;
    }
    else G = ShadowEffect(x , nu);


    double P, W2;

    aux =   x*x;
    P   =   1 - 1.85*x + 2.45*aux - 2.35*aux*x + aux*aux;
    G   *=  (medium_->GetNucCharge().at(component_) + (medium_->GetAtomicNum().at(component_) - medium_->GetNucCharge().at(component_))*P);
    W2  =   medium_->GetAverageNucleonWeight().at(component_)*medium_->GetAverageNucleonWeight().at(component_)
            - Q2 + 2*medium_->GetAverageNucleonWeight().at(component_)*particle_->GetEnergy()*v_;

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

    xp  =   (Q2 + m2p)/(Q2 + m2p + W2 - medium_->GetAverageNucleonWeight().at(component_)*medium_->GetAverageNucleonWeight().at(component_));
    xr  =   (Q2 + m2r)/(Q2 + m2r + W2 - medium_->GetAverageNucleonWeight().at(component_)*medium_->GetAverageNucleonWeight().at(component_));
    F2p =   cp*pow(xp, ap)*pow(1 - x, bp);
    F2r =   cr*pow(xr, ar)*pow(1 - x, br);
    F2  =   (Q2/(Q2 + m2o))*(F2p + F2r)*G;
    R2  =   (2*(1 + R));



    aux =   ME*RE/Q2;
    aux *=  aux*(1 - v_ - medium_->GetAverageNucleonWeight().at(component_)*x*v_/(2*particle_->GetEnergy()) +
                 (1 - 2*particle_->GetMass()*particle_->GetMass()/Q2)
                 *v_*v_*(1 + 4*medium_->GetAverageNucleonWeight().at(component_)*medium_->GetAverageNucleonWeight().at(component_)*x*x/Q2)/R2);

    return (4*PI*F2/v_)*aux;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Photonuclear::FunctionToIntegralALLM97(double Q2)
{
    double x, aux, nu, G, F2, R2;

    nu  =   v_*particle_->GetEnergy();
    x   =   Q2/(2*medium_->GetAverageNucleonWeight().at(component_)*nu);

    if(medium_->GetNucCharge().at(component_)==1)
    {
        G   =   1;
    }
    else G = ShadowEffect(x , nu);

    double P, W2;

    aux =   x*x;
    P   =   1 - 1.85*x + 2.45*aux - 2.35*aux*x + aux*aux;
    G   *=  (medium_->GetNucCharge().at(component_) + (medium_->GetAtomicNum().at(component_) - medium_->GetNucCharge().at(component_))*P);
    W2  =   medium_->GetAverageNucleonWeight().at(component_)*medium_->GetAverageNucleonWeight().at(component_)
            - Q2 + 2*medium_->GetAverageNucleonWeight().at(component_)*particle_->GetEnergy()*v_;

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

    xp  =   (Q2 + m2p)/(Q2 + m2p + W2 - medium_->GetAverageNucleonWeight().at(component_)*medium_->GetAverageNucleonWeight().at(component_));
    xr  =   (Q2 + m2r)/(Q2 + m2r + W2 - medium_->GetAverageNucleonWeight().at(component_)*medium_->GetAverageNucleonWeight().at(component_));
    F2p =   cp*pow(xp, ap)*pow(1 - x, bp);
    F2r =   cr*pow(xr, ar)*pow(1 - x, br);
    F2  =   (Q2/(Q2 + m2o))*(F2p + F2r)*G;
    R2  =   (2*(1 + R));

    aux =   ME*RE/Q2;
    aux *=  aux*(1 - v_ - medium_->GetAverageNucleonWeight().at(component_)*x*v_/(2*particle_->GetEnergy()) +
                 (1 - 2*particle_->GetMass()*particle_->GetMass()/Q2)
                 *v_*v_*(1 + 4*medium_->GetAverageNucleonWeight().at(component_)*medium_->GetAverageNucleonWeight().at(component_)*x*x/Q2)/R2);

    return (4*PI*F2/v_)*aux;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



double Photonuclear::FunctionToIntegralButMik(double Q2)
{
    double x, aux, nu, G, F2, R2;

    nu  =   v_*particle_->GetEnergy();
    x   =   Q2/(2*medium_->GetAverageNucleonWeight().at(component_)*nu);

    if(medium_->GetNucCharge().at(component_)==1)
    {
        G   =   1;
    }
    else G = ShadowEffect(x , nu);

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
    F2  =   G*(medium_->GetNucCharge().at(component_)*F2p + (medium_->GetAtomicNum().at(component_) - medium_->GetNucCharge().at(component_))*F2n);
    R2  =   (2*(1 + R));


    aux =   ME*RE/Q2;
    aux *=  aux*(1 - v_ - medium_->GetAverageNucleonWeight().at(component_)*x*v_/(2*particle_->GetEnergy()) +
                 (1 - 2*particle_->GetMass()*particle_->GetMass()/Q2)
                 *v_*v_*(1 + 4*medium_->GetAverageNucleonWeight().at(component_)*medium_->GetAverageNucleonWeight().at(component_)*x*x/Q2)/R2);

    return (4*PI*F2/v_)*aux;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
void Photonuclear::SetParametrization(int parametrization)
{
    // Now: parametrization_ = 1,  Former: form=1 and bb=1 Kokoulin
    // Now: parametrization_ = 2,  Former: form=2 and bb=1 Kokoulin + hard component
    // Now: parametrization_ = 3,  Former: form=1 and bb=2 Rhode
    // Now: parametrization_ = 4,  Former: form=2 and bb=2 Rhode + hard component
    // Now: parametrization_ = 5,  Former: form=1 and bb=3 Bezrukov/Bugaev
    // Now: parametrization_ = 6,  Former: form=2 and bb=3 Bezrukov/Bugaev + hard component
    // Now: parametrization_ = 7,  Former: form=1 and bb=4 Zeus
    // Now: parametrization_ = 8,  Former: form=2 and bb=4 Zeus + hard component
    // Now: parametrization_ = 9,  Former: form=3 and bb=1 shadow=1 ALLM 91
    // Now: parametrization_ = 10, Former: form=3 and bb=1 shadow=2 ALLM 91
    // Now: parametrization_ = 11, Former: form=3 and bb=2 shadow=1 ALLM 97
    // Now: parametrization_ = 12, Former: form=3 and bb=2 shadow=2 ALLM 97
    // Now: parametrization_ = 13, Former: form=4 and bb=1 shadow=1 Butkevich/Mikhailov
    // Now: parametrization_ = 14, Former: form=4 and bb=1 shadow=2 Butkevich/Mikhailov

    parametrization_    =   parametrization;

    switch(parametrization)
    {
        case 1:     //Kokoulin
            parametrization_family_ =   1;
            hard_component_         =   false;
            break;
        case 2:     //Kokoulin + hard component
            parametrization_family_ =   1;
            hard_component_         =   true;
            break;
        case 3:     //Rhode
            parametrization_family_ =   2;
            hard_component_         =   false;
            break;
        case 4:     //Rhode + hard component
            parametrization_family_ =   2;
            hard_component_         =   true;
            break;
        case 5:     //Bezrukov/Bugaev
            parametrization_family_ =   3;
            hard_component_         =   false;
            break;
        case 6:     //Bezrukov/Bugaev + hard component
            parametrization_family_ =   3;
            hard_component_         =   true;
            break;
        case 7:     //Zeus
            parametrization_family_ =   4;
            hard_component_         =   false;
            break;
        case 8:     //Zeus + hard component
            parametrization_family_ =   4;
            hard_component_         =   true;
            break;
        case 9:     //ALLM 91 shadow 1
            parametrization_family_ =   5;
            shadow_                 =   1;
            break;
        case 10:    //ALLM 91 shadow 2
            parametrization_family_ =   5;
            shadow_                 =   2;
            break;
        case 11:    //ALLM 97 shadow 1
            parametrization_family_ =   6;
            shadow_                 =   1;
            break;
        case 12:    //ALLM 97 shadow 2
            parametrization_family_ =   6;
            shadow_                 =   2;
            break;
        case 13:    //Butkevich/Mikhailov shadow 1
            parametrization_family_ =   7;
            shadow_                 =   1;
            break;
        case 14:    //Butkevich/Mikhailov shadow 2
            parametrization_family_ =   7;
            shadow_                 =   2;
            break;
        default:
            log_warn("Parametrization %i not supported. Set to 12 (icecube default)",parametrization_);
            parametrization_family_ =   6;
            shadow_                 =   2;
    }

    if(do_photo_interpolation_)
    {
        log_warn("photo-interpolation enabled before choosing the parametrization. Building the tables again");
        DisablePhotoInterpolation();
        EnablePhotoInterpolation();
    }

    if(do_dedx_Interpolation_)
    {
        log_warn("dEdx-interpolation enabled before choosing the parametrization. Building the tables again");
        DisableDEdxInterpolation();
        EnableDEdxInterpolation();
    }
    if(do_dndx_Interpolation_)
    {
        log_warn("dNdx-interpolation enabled before choosing the parametrization. Building the tables again");
        DisableDNdxInterpolation();
        EnableDNdxInterpolation();
    }
}

void Photonuclear::SetComponent(int component) {
	component_ = component;
}

void Photonuclear::SetDedxInterpolant(Interpolant* dedxInterpolant) {
	dedx_interpolant_ = dedxInterpolant;
}

void Photonuclear::SetDndxIntegral(std::vector<Integral*> dndxIntegral) {
	dndx_integral_ = dndxIntegral;
}

void Photonuclear::SetDndxInterpolant1d(
		std::vector<Interpolant*> dndxInterpolant1d) {
	dndx_interpolant_1d_ = dndxInterpolant1d;
}

void Photonuclear::SetDndxInterpolant2d(
		std::vector<Interpolant*> dndxInterpolant2d) {
	dndx_interpolant_2d_ = dndxInterpolant2d;
}

void Photonuclear::SetHmax(int hmax) {
	hmax_ = hmax;
}

void Photonuclear::SetInitHardbb(bool initHardbb) {
	init_hardbb_ = initHardbb;
}

void Photonuclear::SetInitMeasured(bool initMeasured) {
	init_measured_ = initMeasured;
}

void Photonuclear::SetIntegral(Integral* integral) {
	integral_ = integral;
}

void Photonuclear::SetIntegralForDEdx(Integral* integralForDEdx) {
	integral_for_dEdx_ = integralForDEdx;
}

void Photonuclear::SetInterpolantHardBb(
		std::vector<Interpolant*> interpolantHardBb) {
	interpolant_hardBB_ = interpolantHardBb;
}

void Photonuclear::SetInterpolantMeasured(
		Interpolant* interpolantMeasured) {
	interpolant_measured_ = interpolantMeasured;
}

void Photonuclear::SetProbForComponent(std::vector<double> probForComponent) {
	prob_for_component_ = probForComponent;
}

void Photonuclear::SetV(double v) {
	v_ = v;
}

void Photonuclear::SetShadow(int shadow){
    shadow_  =   shadow;
}

void Photonuclear::SetHardComponent(bool hard){
    hard_component_  =   hard;
}
