#include "PROPOSAL/Epairproduction.h"
#include <algorithm>
#include "boost/bind.hpp"
#include "PROPOSAL/Output.h"


using namespace std;

namespace po	= boost::program_options;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Epairproduction::CalculatedEdx()
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

    for(int i=0; i<medium_->GetNumComponents(); i++)
    {
        SetIntegralLimits(i);
        double r1   =   0.8;
        double rUp  =   vUp_*(1-HALF_PRECISION);
        bool rflag  =   false;

        if(r1<rUp)
        {
            if(2*FunctionToDEdxIntegral(r1)<FunctionToDEdxIntegral(rUp))
            {
                rflag   =   true;
            }
        }

        if(rflag)
        {
            if(r1>vUp_)
            {
                r1  =   vUp_;
            }

            if(r1<vMin_)
            {
                r1  =   vMin_;
            }

            sum         +=  integral_for_dEdx_->Integrate(vMin_, r1, boost::bind(&Epairproduction::FunctionToDEdxIntegral, this, _1),4);
            reverse_    =   true;
            double r2   =   max(1-vUp_, COMPUTER_PRECISION);

            if(r2>1-r1)
            {
                r2  =   1-r1;
            }

            sum         +=  integral_for_dEdx_->Integrate(1-vUp_, r2, boost::bind(&Epairproduction::FunctionToDEdxIntegral, this, _1),2)
                        +   integral_for_dEdx_->Integrate(r2, 1-r1, boost::bind(&Epairproduction::FunctionToDEdxIntegral, this, _1),4);

            reverse_    =   false;
        }

        else
        {
            sum +=  integral_for_dEdx_->Integrate(vMin_, vUp_, boost::bind(&Epairproduction::FunctionToDEdxIntegral, this, _1),4);
        }
    }

    return multiplier_*particle_->GetEnergy()*sum;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Epairproduction::CalculatedNdx()
{
    if(multiplier_<=0)
    {
        return 0;
    }


    sum_of_rates_ = 0;

    for(int i=0; i<medium_->GetNumComponents(); i++){
        if(do_dndx_Interpolation_)
        {
            sum_of_rates_ += max(dndx_interpolant_1d_.at(i)->Interpolate(particle_->GetEnergy()), 0.0);
        }
        else
        {
            SetIntegralLimits(i);
            sum_of_rates_ += dndx_integral_.at(i)->Integrate(vUp_,vMax_, boost::bind(&Epairproduction::FunctionToDNdxIntegral, this, _1),4);
        }
    }

    return sum_of_rates_;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Epairproduction::CalculatedNdx(double rnd)
{
    if(multiplier_<=0)
    {
        return 0;
    }

    // The random number will be stored to be able
    // to check if dNdx is already calculated for this random number.
    // This avoids a second calculation in CalculateStochaticLoss

    rnd_    =   rnd;

    sum_of_rates_  =   0;

    for(int i=0; i<medium_->GetNumComponents(); i++){
        if(do_dndx_Interpolation_)
        {
            prob_for_component_.at(i) = max(dndx_interpolant_1d_.at(i)->Interpolate(particle_->GetEnergy()), 0.0);
        }
        else
        {
            SetIntegralLimits(i);
            prob_for_component_.at(i) = dndx_integral_.at(i)->IntegrateWithLog(vUp_,vMax_, boost::bind(&Epairproduction::FunctionToDNdxIntegral, this, _1),rnd);
        }

        sum_of_rates_ += prob_for_component_.at(i);
    }

    return sum_of_rates_;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Epairproduction::CalculateStochasticLoss(double rnd1, double rnd2)
{
    if(rnd1 != rnd_ )
    {
        CalculatedNdx(rnd1);
        log_warn("CalculatedNdx was not called! rnd1 and rnd_ don´t match! \n Calculationg again with rnd1=%f \t rnd2=%f",rnd1,rnd2);
    }

    return CalculateStochasticLoss(rnd2);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-----------------------Enable and disable interpolation---------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Epairproduction::EnableDNdxInterpolation(std::string path, bool raw)
{

    if(do_dndx_Interpolation_)return;

    EnableEpairInterpolation(path,raw);

    bool storing_failed =   false;
    bool reading_worked =   true;

    if(!path.empty())
    {
        stringstream filename;
        filename<<path<<"/Epair_dNdx_particle_"<<particle_->GetName()
                <<"_med_"<<medium_->GetName()
                <<medium_->GetMassDensity()
                <<"_ecut_"<<cut_settings_->GetEcut()
                <<"_vcut_"<<cut_settings_->GetVcut()
                <<"_lpm_"<<lpm_effect_enabled_
                <<"_multiplier_"<<multiplier_;

        if(!raw)
            filename<<".txt";

        dndx_interpolant_1d_.resize(medium_->GetNumComponents());
        dndx_interpolant_2d_.resize(medium_->GetNumComponents());

        if( FileExist(filename.str()) )
        {
            log_info("Epairproduction parametrisation tables (dNdx) will be read from file:\t%s",filename.str().c_str());
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
                reading_worked = dndx_interpolant_2d_.at(i)->Load(input,raw);
                reading_worked = dndx_interpolant_1d_.at(i)->Load(input,raw);

            }
            input.close();
        }
        if(!FileExist(filename.str()) || !reading_worked )
        {
            if(!reading_worked)
            {
                log_info("File %s is corrupted! Write it again!",filename.str().c_str());
            }

            log_info("Epairproduction parametrisation tables (dNdx) will be saved to file:\t%s",filename.str().c_str());

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

                    dndx_interpolant_2d_.at(i) =    new Interpolant(NUM1, particle_->GetLow(), BIGENERGY,  NUM1, 0, 1, boost::bind(&Epairproduction::FunctionToBuildDNdxInterpolant2D, this, _1 , _2), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false, order_of_interpolation_, true, false, false);
                    dndx_interpolant_1d_.at(i) =    new Interpolant(NUM1, particle_->GetLow(), BIGENERGY,  boost::bind(&Epairproduction::FunctionToBuildDNdxInterpolant1D, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, true, false, false);

                    dndx_interpolant_2d_.at(i)->Save(output,raw);
                    dndx_interpolant_1d_.at(i)->Save(output,raw);

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
        dndx_interpolant_1d_.resize( medium_->GetNumComponents() );
        dndx_interpolant_2d_.resize( medium_->GetNumComponents() );
        for(int i=0; i<(medium_->GetNumComponents()); i++)
        {
            component_ = i;
            dndx_interpolant_2d_.at(i) =    new Interpolant(NUM1, particle_->GetLow(), BIGENERGY,  NUM1, 0, 1, boost::bind(&Epairproduction::FunctionToBuildDNdxInterpolant2D, this, _1 , _2), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false, order_of_interpolation_, true, false, false);
            dndx_interpolant_1d_.at(i) =    new Interpolant(NUM1, particle_->GetLow(), BIGENERGY,  boost::bind(&Epairproduction::FunctionToBuildDNdxInterpolant1D, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, true, false, false);

        }
        particle_->SetEnergy(energy);
    }

    do_dndx_Interpolation_=true;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Epairproduction::EnableDEdxInterpolation(std::string path, bool raw)
{

    if(do_dedx_Interpolation_)return;

    EnableEpairInterpolation(path,raw);

    bool reading_worked =   true;
    bool storing_failed =   false;

    if(!path.empty())
    {
        stringstream filename;
        filename<<path<<"/Epair_dEdx_particle_"<<particle_->GetName()
                <<"_med_"<<medium_->GetName()
                <<medium_->GetMassDensity()
                <<"_ecut_"<<cut_settings_->GetEcut()
                <<"_vcut_"<<cut_settings_->GetVcut()
                <<"_lpm_"<<lpm_effect_enabled_
                <<"_multiplier_"<<multiplier_;

        if(!raw)
            filename<<".txt";

        if( FileExist(filename.str()) )
        {
            log_info("Epairproduction parametrisation tables (dEdx) will be read from file:\t%s",filename.str().c_str());
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
            reading_worked = dedx_interpolant_->Load(input,raw);

            input.close();
        }
        if(!FileExist(filename.str()) || !reading_worked )
        {
            if(!reading_worked)
            {
                log_info("File %s is corrupted! Write it again!",filename.str().c_str());
            }

            log_info("Epairproduction parametrisation tables (dEdx) will be saved to file:\t%s",filename.str().c_str());

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

                dedx_interpolant_ = new Interpolant(NUM1, particle_->GetLow(), BIGENERGY, boost::bind(&Epairproduction::FunctionToBuildDEdxInterpolant, this, _1), order_of_interpolation_, true, false, true, order_of_interpolation_, false, false, false);

                dedx_interpolant_->Save(output,raw);
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

        dedx_interpolant_ = new Interpolant(NUM1, particle_->GetLow(), BIGENERGY, boost::bind(&Epairproduction::FunctionToBuildDEdxInterpolant, this, _1), order_of_interpolation_, true, false, true, order_of_interpolation_, false, false, false);

        particle_->SetEnergy(energy);
    }

    do_dedx_Interpolation_=true;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Epairproduction::EnableEpairInterpolation(std::string path, bool raw)
{

    if(do_epair_interpolation_)return;

    bool storing_failed =   false;
    bool reading_worked =   true;

    if(!path.empty())
    {
        stringstream filename;
        filename<<path<<"/Epair_particle_"<<particle_->GetName()
                <<"_med_"<<medium_->GetName()
                <<medium_->GetMassDensity()
                <<"_ecut_"<<cut_settings_->GetEcut()
                <<"_vcut_"<<cut_settings_->GetVcut()
                <<"_lpm_"<<lpm_effect_enabled_
                <<"_multiplier_"<<multiplier_;

        if(!raw)
            filename<<".txt";

        epair_interpolant_.resize( medium_->GetNumComponents() );

        if( FileExist(filename.str()) )
        {
            log_info("Epairproduction parametrisation tables will be read from file:\t%s",filename.str().c_str());
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
                epair_interpolant_.at(i) = new Interpolant();
                reading_worked = epair_interpolant_.at(i)->Load(input,raw);

            }
            input.close();
        }
        if(!FileExist(filename.str()) || !reading_worked )
        {
            if(!reading_worked)
            {
                log_info("File %s is corrupted! Write it again!",filename.str().c_str());
            }

            log_info("Epairproduction parametrisation tables will be saved to file:\t%s",filename.str().c_str());

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

                    epair_interpolant_.at(i)   = new Interpolant(NUM1, particle_->GetLow(), BIGENERGY, NUM1, 0., 1.,boost::bind(&Epairproduction::FunctionToBuildEpairInterpolant, this, _1 , _2) , order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false, order_of_interpolation_, false, false, false);

                    epair_interpolant_.at(i)->Save(output,raw);

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

        epair_interpolant_.resize( medium_->GetNumComponents() );

        for(int i=0; i < medium_->GetNumComponents() ; i++)
        {
            component_ = i;
            epair_interpolant_.at(i)   = new Interpolant(NUM1, particle_->GetLow(), BIGENERGY, NUM1, 0., 1.,boost::bind(&Epairproduction::FunctionToBuildEpairInterpolant, this, _1 , _2) , order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false, order_of_interpolation_, false, false, false);

        }
        particle_->SetEnergy(energy);
    }

    do_epair_interpolation_=true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Epairproduction::DisableDNdxInterpolation()
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

void Epairproduction::DisableDEdxInterpolation()
{
    delete dedx_interpolant_;
    dedx_interpolant_   =   NULL;
    do_dedx_Interpolation_  =   false;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

void Epairproduction::DisableEpairInterpolation()
{
    for(unsigned int i=0; i < epair_interpolant_.size(); i++)
    {
        delete epair_interpolant_.at(i);
    }
    epair_interpolant_.clear();
    do_epair_interpolation_  =   false;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------Set and validate options--------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


boost::program_options::options_description Epairproduction::CreateOptions()
{
    po::options_description epairproduction("Epairproduction options");
    epairproduction.add_options()
        ("epairproduction.lpm",              po::value<bool>(&lpm_effect_enabled_)->implicit_value(false),     "Enables   Landau-Pomeranchuk-Migdal supression")
        ("epairproduction.interpol_dedx",    po::value<bool>(&do_dedx_Interpolation_)->implicit_value(false),  "Enables interpolation for dEdx")
        ("epairproduction.interpol_dndx",    po::value<bool>(&do_dndx_Interpolation_)->implicit_value(false),  "Enables interpolation for dNdx")
        ("epairproduction.multiplier",       po::value<double>(&multiplier_)->default_value(1.),               "modify the cross section by this factor")
        ("epairproduction.interpol_order",   po::value<int>(&order_of_interpolation_)->default_value(5),       "number of interpolation points");

   return epairproduction;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Epairproduction::ValidateOptions()
{
    if(order_of_interpolation_ < 2)
    {
        order_of_interpolation_ = 5;
        cerr<<"Epairproduction: Order of Interpolation is not a vaild number\t"<<"Set to 5"<<endl;
    }
    if(order_of_interpolation_ > 6)
    {
        cerr<<"Epairproduction: Order of Interpolation is set to "<<order_of_interpolation_
            <<".\t Note a order of interpolation > 6 will slow down the program"<<endl;
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Epairproduction::Epairproduction()
    :component_             ( 0 )
    ,v_                     ( 0 )
    ,reverse_               ( false )
    ,eLpm_                  ( 0 )
    ,do_epair_interpolation_( false )
    ,dndx_integral_         ( )
    ,dndx_interpolant_1d_   ( )
    ,dndx_interpolant_2d_   ( )
    ,epair_interpolant_     ( )
    ,prob_for_component_    ( )
{
    integral_             = new Integral(IROMB, IMAXS, IPREC);
    integral_for_dEdx_    = new Integral(IROMB, IMAXS, IPREC);
    dedx_interpolant_     = NULL;
    name_                 = "Epairproduction";

}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Epairproduction::Epairproduction(const Epairproduction &epair)
    :CrossSections                      ( epair )
    ,component_                         ( epair.component_ )
    ,v_                                 ( epair.v_ )
    ,reverse_                           ( epair.reverse_ )
    ,eLpm_                              ( epair.eLpm_)
    ,do_epair_interpolation_            ( epair.do_epair_interpolation_ )
    ,integral_                          ( new Integral(*epair.integral_) )
    ,integral_for_dEdx_                 ( new Integral(*epair.integral_for_dEdx_) )
    ,prob_for_component_                ( epair.prob_for_component_)

{    
    if(epair.dedx_interpolant_ != NULL)
    {
        dedx_interpolant_ = new Interpolant(*epair.dedx_interpolant_) ;
    }
    else
    {
        dedx_interpolant_ = NULL;
    }

    dndx_integral_.resize( epair.dndx_integral_.size() );
    dndx_interpolant_1d_.resize( epair.dndx_interpolant_1d_.size() );
    dndx_interpolant_2d_.resize( epair.dndx_interpolant_2d_.size() );
    epair_interpolant_.resize( epair.epair_interpolant_.size() );

    for(unsigned int i =0; i<epair.dndx_integral_.size(); i++)
    {
        dndx_integral_.at(i) = new Integral( *epair.dndx_integral_.at(i) );
    }
    for(unsigned int i =0; i<epair.dndx_interpolant_1d_.size(); i++)
    {
        dndx_interpolant_1d_.at(i) = new Interpolant( *epair.dndx_interpolant_1d_.at(i) );
    }
    for(unsigned int i =0; i<epair.dndx_interpolant_2d_.size(); i++)
    {
        dndx_interpolant_2d_.at(i) = new Interpolant( *epair.dndx_interpolant_2d_.at(i) );
    }
    for(unsigned int i =0; i<epair.epair_interpolant_.size(); i++)
    {
        epair_interpolant_.at(i) = new Interpolant( *epair.epair_interpolant_.at(i) );
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Epairproduction::Epairproduction(Particle* particle,
                             Medium* medium,
                             EnergyCutSettings* cut_settings)
    :CrossSections(particle, medium, cut_settings)
    ,v_                     ( 0 )
    ,reverse_               ( false )
    ,eLpm_                  ( 0 )
    ,do_epair_interpolation_( false )
    ,dndx_integral_         ( )
    ,dndx_interpolant_1d_   ( )
    ,dndx_interpolant_2d_   ( )
    ,epair_interpolant_     ( )
    ,prob_for_component_    ( )
{
    name_                       = "Epairproduction";
    vMax_                       = 0;
    vUp_                        = 0;
    vMin_                       = 0;
    ebig_                       = BIGENERGY;
    do_dedx_Interpolation_      = false;
    do_dndx_Interpolation_      = false;
    multiplier_                 = 1.;
    parametrization_            = 1;
    lpm_effect_enabled_         = false;
    init_lpm_effect_            = true;
    component_                  = 0;


    integral_             = new Integral(IROMB, IMAXS, IPREC);
    integral_for_dEdx_    = new Integral(IROMB, IMAXS, IPREC);

    dndx_integral_.resize(medium_->GetNumComponents());


    for(int i =0 ; i<medium_->GetNumComponents();i++){
        dndx_integral_.at(i) = new Integral(IROMB, IMAXS, IPREC);
    }

    prob_for_component_.resize(medium_->GetNumComponents());
    dedx_interpolant_     = NULL;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Epairproduction& Epairproduction::operator=(const Epairproduction &epair)
{

    if (this != &epair)
    {
      Epairproduction tmp(epair);
      swap(tmp);
    }
    return *this;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Epairproduction::operator==(const Epairproduction &epair) const
{
    if( this->CrossSections::operator !=(epair) )                           return false;
    if( component_                  !=  epair.component_)                   return false;
    if( v_                          !=  epair.v_)                           return false;
    if( reverse_                    !=  epair.reverse_)                     return false;
    if( eLpm_                       !=  epair.eLpm_)                        return false;
    if( do_epair_interpolation_     !=  epair.do_epair_interpolation_ )     return false;
    if( *integral_                  != *epair.integral_)                    return false;
    if( *integral_for_dEdx_         != *epair.integral_for_dEdx_)           return false;
    if( prob_for_component_.size()  !=  epair.prob_for_component_.size())   return false;
    if( dndx_integral_.size()       !=  epair.dndx_integral_.size())        return false;
    if( dndx_interpolant_1d_.size() !=  epair.dndx_interpolant_1d_.size())  return false;
    if( dndx_interpolant_2d_.size() !=  epair.dndx_interpolant_2d_.size())  return false;
    if( epair_interpolant_.size()   !=  epair.epair_interpolant_.size())    return false;


    for(unsigned int i =0; i<epair.dndx_integral_.size(); i++)
    {
        if( *dndx_integral_.at(i)       != *epair.dndx_integral_.at(i) )        return false;
    }
    for(unsigned int i =0; i<epair.dndx_interpolant_1d_.size(); i++)
    {
        if( *dndx_interpolant_1d_.at(i) != *epair.dndx_interpolant_1d_.at(i) )  return false;
    }
    for(unsigned int i =0; i<epair.dndx_interpolant_2d_.size(); i++)
    {
        if( *dndx_interpolant_2d_.at(i) != *epair.dndx_interpolant_2d_.at(i) )  return false;
    }
    for(unsigned int i =0; i<epair.prob_for_component_.size(); i++)
    {
        if( prob_for_component_.at(i) != epair.prob_for_component_.at(i) )      return false;
    }
    for(unsigned int i =0; i<epair.epair_interpolant_.size(); i++)
    {
        if( *epair_interpolant_.at(i) != *epair.epair_interpolant_.at(i) )      return false;
    }

    if( dedx_interpolant_ != NULL && epair.dedx_interpolant_ != NULL)
    {
        if( *dedx_interpolant_ != *epair.dedx_interpolant_)                     return false;
    }
    else if( dedx_interpolant_ != epair.dedx_interpolant_)                      return false;

    //else
    return true;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Epairproduction::operator!=(const Epairproduction &epair) const
{
    return !(*this == epair);

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Epairproduction::swap(Epairproduction &epair)
{
    using std::swap;

    this->CrossSections::swap(epair);

    swap(component_,epair.component_);
    swap(v_,epair.v_);
    swap(reverse_,epair.reverse_);
    swap(eLpm_,epair.eLpm_);
    swap(do_epair_interpolation_,epair.do_epair_interpolation_);
    integral_for_dEdx_->swap(*epair.integral_for_dEdx_);
    integral_->swap(*epair.integral_);

    if( dedx_interpolant_ != NULL && epair.dedx_interpolant_ != NULL)
    {
        dedx_interpolant_->swap(*epair.dedx_interpolant_);
    }
    else if( dedx_interpolant_ == NULL && epair.dedx_interpolant_ != NULL)
    {
        dedx_interpolant_ = new Interpolant(*epair.dedx_interpolant_);
        epair.dedx_interpolant_ = NULL;
    }
    else if( dedx_interpolant_ != NULL && epair.dedx_interpolant_ == NULL)
    {
        epair.dedx_interpolant_ = new Interpolant(*dedx_interpolant_);
        dedx_interpolant_ = NULL;
    }

    prob_for_component_.swap(epair.prob_for_component_);
    dndx_integral_.swap(epair.dndx_integral_);
    dndx_interpolant_1d_.swap(epair.dndx_interpolant_1d_);
    dndx_interpolant_2d_.swap(epair.dndx_interpolant_2d_);
    epair_interpolant_.swap(epair.epair_interpolant_);

}


ostream& operator<<(std::ostream& os, Epairproduction const &epair)
{
    os<<"---------------------------Epairproduction( "<<&epair<<" )---------------------------"<<std::endl;
    os<< static_cast <const CrossSections &>( epair ) << endl;
    os<< "------- Class Specific: " << endl;
    os<<"\tintegral:\t\t"<<epair.integral_ << endl;
    os<<"\tdedx_integral:\t"<<epair.integral_for_dEdx_ << endl;
    for(unsigned int i=0;i<epair.dndx_integral_.size();i++)
    {
        os<<"\t\tadress:\t\t"<<epair.dndx_integral_.at(i)<<endl;
    }
    os<<endl;
    os<<"\tdo_dedx_Interpolation:\t\t"<<epair.do_dedx_Interpolation_<<endl;
    os<<"\tdedx_interpolant:\t\t"<<epair.dedx_interpolant_<<endl;
    os<<endl;
    os<<"\tdo_dndx_Interpolation:\t\t"<<epair.do_dndx_Interpolation_<<endl;
    os<<"\tdndx_interpolant_1d:\t\t"<<epair.dndx_interpolant_1d_.size()<<endl;
    for(unsigned int i=0;i<epair.dndx_interpolant_1d_.size();i++)
    {
        os<<"\t\tadress:\t\t"<<epair.dndx_interpolant_1d_.at(i)<<endl;
    }
    os<<"\tdndx_interpolant_2d:\t\t"<<epair.dndx_interpolant_2d_.size()<<endl;
    for(unsigned int i=0;i<epair.dndx_interpolant_2d_.size();i++)
    {
        os<<"\t\tadress:\t\t"<<epair.dndx_interpolant_2d_.at(i)<<endl;
    }
    os<<endl;
    os<<"\tdo_epair_interpolation:\t"<<epair.do_epair_interpolation_<<endl;
    os<<"\tepair_interpolant:\t\t"<<epair.epair_interpolant_.size()<<endl;
    for(unsigned int i=0;i<epair.epair_interpolant_.size();i++)
    {
        os<<"\t\tadress:\t\t"<<epair.epair_interpolant_.at(i)<<endl;
    }
    os<<endl;
    os<<"\tTemp. variables: " << endl;
    os<< "\t\tcomponent:\t\t" << epair.component_<<endl;
    os<< "\t\tv:\t\t" << epair.v_<<endl;
    os<< "\t\treverse:\t\t" << epair.reverse_<<endl;
    os<< "\t\teLpm:\t\t" << epair.eLpm_<<endl;
    os<<"\t\tprob_for_component:\t"<<epair.prob_for_component_.size()<<endl;
    for(unsigned int i=0;i<epair.prob_for_component_.size();i++)
    {
        os<<"\t\t\tvalue:\t\t"<<epair.prob_for_component_.at(i)<<endl;
    }
    os<<"-----------------------------------------------------------------------------------------------";
    return os;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------Cross section / limit / lpm---------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Epairproduction::CalculateStochasticLoss(double rnd)
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

                return (particle_->GetEnergy())*(vUp_*exp(dndx_interpolant_2d_.at(i)->FindLimit((particle_->GetEnergy()), (rnd_)*prob_for_component_.at(i))*log(vMax_/vUp_)));
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


void Epairproduction::SetIntegralLimits(int component)
{

    component_ = component;
    double aux;

    vMin_    =   4*ME/particle_->GetEnergy();
    vMax_    =   1 - (3./4)*SQRTE*(particle_->GetMass()/particle_->GetEnergy())
                * pow(medium_->GetNucCharge().at(component) , 1./3);
    aux      =   particle_->GetMass()/particle_->GetEnergy();
    aux      =   1-6*aux*aux;
    vMax_    =   min(vMax_, aux);
    vMax_    =   min(vMax_, 1-particle_->GetMass()/particle_->GetEnergy());

    if(vMax_ < vMin_)
    {
        vMax_    =   vMin_;
    }

    vUp_     =   min(vMax_, cut_settings_->GetCut(particle_->GetEnergy()));

    if(vUp_ < vMin_)
    {
        vUp_     =   vMin_;
    }

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Epairproduction::lpm(double r2, double b, double x)
{


    if(init_lpm_effect_)
    {
        init_lpm_effect_        =   false;
        double sum  =   0;

        for(int i=0; i<medium_->GetNumComponents(); i++)
        {
            sum +=  medium_->GetNucCharge().at(i)*medium_->GetNucCharge().at(i)
                    *log(3.25*medium_->GetLogConstant().at(i)*pow(medium_->GetNucCharge().at(i), -1./3));
        }

        eLpm_    =   particle_->GetMass()/(ME*RE);
        eLpm_    *=  (eLpm_*eLpm_)*ALPHA*particle_->GetMass()
                /(2*PI*medium_->GetMolDensity()*particle_->GetCharge()*particle_->GetCharge()*sum);
    }

    double A, B, C, D, E, s;
    double s2, s36, s6, d1, d2, atan_, log1, log2;

    s       =   sqrt(eLpm_/(particle_->GetEnergy()*v_*(1 - r2)))/4;
    s6      =   6*s;
    atan_   =   s6*(x + 1);

    if(atan_>1/COMPUTER_PRECISION)
    {
        return 1;
    }

    s2      =   s*s;
    s36     =   36*s2;
    d1      =   s6/(s6 + 1);
    d2      =   s36/(s36 + 1);
    atan_   =   atan(atan_) - PI/2;
    log1    =   log((s36*(1 + x)*(1 + x) + 1)/(s36*x*x));
    log2    =   log((s6*(1 + x) + 1)/(s6*x));
    A       =   0.5*d2*(1 + 2*d2*x)*log1 - d2 + 6*d2*s*(1 + ((s36 - 1)/(s36 + 1))*x)*atan_;
    B       =   d1*(1 + d1*x)*log2 - d1;
    C       =   -d2*d2*x*log1 + d2 - (d2*d2*(s36 - 1)/(6*s))*x*atan_;
    D       =   d1-d1*d1*x*log2;
    E       =   -s6*atan_;

    return ((1 + b)*(A + (1 + r2)*B) + b*(C + (1 + r2)*D) + (1 - r2)*E)
            /(((2 + r2)*(1 +b ) + x*(3 + r2))*log(1 +1 /x) + (1 - r2 - b)/(1 + x) - (3 + r2));


}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Epairproduction::EPair(double v, int component)
{

    if(do_epair_interpolation_)
    {
        SetIntegralLimits(component);

        if(v>=vUp_)
        {
            return max(epair_interpolant_.at(component)->Interpolate(particle_->GetEnergy(), log(v/vUp_)/log(vMax_/vUp_)), 0.0);
        }
    }

    double rMax, aux, aux2;

    component_  =   component;
    v_          =   v;
    aux         =   1 - (4*ME)/(particle_->GetEnergy()*v_);
    aux2        =   1 - (6*particle_->GetMass()*particle_->GetMass())/(particle_->GetEnergy()*particle_->GetEnergy()*(1 - v_));

    if(aux>0 && aux2>0)
    {
        rMax    =   sqrt(aux)*aux2;
    }
    else
    {
        rMax    =   0;
    }

    aux =   max(1 - rMax , COMPUTER_PRECISION);

    return medium_->GetMolDensity()*medium_->GetAtomInMolecule().at(component_)
           *particle_->GetCharge()*particle_->GetCharge()
           *(integral_->Integrate(1 - rMax, aux, boost::bind(&Epairproduction::FunctionToIntegral, this, _1),2)
                + integral_->Integrate(aux, 1,  boost::bind(&Epairproduction::FunctionToIntegral, this, _1),4));

}
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------Functions to interpolate---------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Epairproduction::FunctionToBuildDNdxInterpolant1D(double energy)
{
    return dndx_interpolant_2d_.at(component_)->Interpolate(energy,1.);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Epairproduction::FunctionToBuildDNdxInterpolant2D(double energy, double v)
{
    particle_->SetEnergy(energy);
    SetIntegralLimits(component_);

    if(vUp_==vMax_)
    {
    return 0.;
    }

    v   =   vUp_*exp(v*log(vMax_/vUp_));

    return dndx_integral_.at(component_)->Integrate(vUp_,v, boost::bind(&Epairproduction::FunctionToDNdxIntegral, this, _1),4);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Epairproduction::FunctionToBuildDEdxInterpolant(double energy)
{
    particle_->SetEnergy(energy);
    return CalculatedEdx();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Epairproduction::FunctionToBuildEpairInterpolant(double energy , double v)
{
    particle_->SetEnergy(energy);
    SetIntegralLimits(component_);

    if(vUp_==vMax_)
    {
    return 0;
    }

    v   =   vUp_*exp(v*log(vMax_/vUp_));

    return EPair(v, component_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------Functions to integrate----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Epairproduction::FunctionToDEdxIntegral(double variable)
{

    if(reverse_)
    {
        variable   =   1-variable;
    }

    return variable*EPair(variable, component_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Epairproduction::FunctionToDNdxIntegral(double variable)
{
    return  multiplier_ * EPair(variable, component_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Epairproduction::FunctionToIntegral(double r)
{

    double Fe, Fm, Le, Lm, Ye, Ym, s, b, k, g1, g2;
    double aux, aux1, aux2, r2, Z3;

    r       =   1-r; // only for integral optimization - do not forget to swap integration limits!
    r2      =   r*r;
    Z3      =   pow(medium_->GetNucCharge().at(component_) , -1./3);
    aux     =   (particle_->GetMass()*v_)/(2*ME);
    aux     *=  aux;
    s       =   aux*(1 - r2)/(1 - v_);
    b       =   (v_*v_)/(2*(1 - v_));
    Ye      =   (5 - r2 + 4*b*(1 + r2))/(2*(1 + 3*b)*log(3 + 1/s) - r2 - 2*b*(2 - r2));
    Ym      =   (4 + r2 + 3*b*(1 + r2))/((1 + r2)*(1.5 + 2*b)*log(3 + s) + 1 - 1.5*r2);
    aux     =   (1.5*ME)/(particle_->GetMass()*Z3);
    aux     *=  aux;
    aux1    =   (1 + s)*(1 + Ye);
    aux2    =   (2*ME*SQRTE*medium_->GetLogConstant().at(component_)*Z3) / (particle_->GetEnergy()*v_*(1 - r2));
    Le      =   log((medium_->GetLogConstant().at(component_)*Z3*sqrt(aux1)) / (1 + aux2*aux1)) - 0.5*log(1 + aux*aux1);
    Lm      =   log(((particle_->GetMass()/(1.5*ME))*medium_->GetLogConstant().at(component_)*Z3*Z3)/(1 + aux2*(1 + s)*(1 + Ym)));

    if ( Le > 0 )
    {
        if (1/s < HALF_PRECISION)
        {
            Fe = (1.5 - r2/2 + b*(1 + r2))/s*Le;
        }
        else
        {
            Fe = (((2 + r2)*(1 + b) + s*(3 + r2))*log(1 + 1/s) + (1 - r2 - b)/(1 + s) - (3 + r2))*Le;
        }
    }
    else
    {
        Fe = 0;
    }

    if ( Le > 0)
    {
        Fm = (((1 + r2)*(1 + 1.5*b) - (1 + 2*b)*(1 - r2)/s)*log(1 + s) + s*(1 - r2 - b)/(1 + s) + (1 + 2*b)*(1 - r2))*Lm;
    }

    else
    {
        Fm = 0;
    }

    if(medium_->GetNucCharge().at(component_)==1)
    {
        g1  =   4.4e-5;
        g2  =   4.8e-5;
    }
    else
    {
        g1  =   1.95e-5;
        g2  =   5.3e-5;
    }

    aux     =   particle_->GetEnergy()/particle_->GetMass();
    k       =   0.058*log(aux/(1 + g2*aux/Z3)) - 0.14;

    if(k<=0)
    {
        k   =   0;
    }

    else
    {
        k   =   (0.073*log(aux/(1 + g1*aux/(Z3*Z3))) - 0.26)/k;
    }

    if(k<0)
    {
        k   =   0;
    }

    aux     =   ALPHA*RE;
    aux     *=  aux/(1.5*PI);
    aux1    =   ME/particle_->GetMass();
    aux1    *=  aux1;

    if(lpm_effect_enabled_)
    {
        aux     *=  2*medium_->GetNucCharge().at(component_)
                *  (medium_->GetNucCharge().at(component_) + k)*((1 - v_)/v_)*lpm(r2, b, s)*(Fe + aux1*Fm);
    }
    else
    {
        aux     *=  2*medium_->GetNucCharge().at(component_)
                *  (medium_->GetNucCharge().at(component_) + k)*((1 - v_)/v_)*(Fe + aux1*Fm);
    }
    if(aux<0)
    {
        aux =   0;
    }

    return aux;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

void Epairproduction::SetParametrization(int parametrization){
    parametrization_ = parametrization;
    log_warn("This has no effect. Till now only one parametrization for Epairproduction implemented");
}

void Epairproduction::SetComponent(int component) {
	component_ = component;
}

void Epairproduction::SetDedxInterpolant(Interpolant* dedxInterpolant) {
	dedx_interpolant_ = dedxInterpolant;
}

void Epairproduction::SetDndxIntegral(std::vector<Integral*> dndxIntegral) {
	dndx_integral_ = dndxIntegral;
}

void Epairproduction::SetDndxInterpolant1d(
		std::vector<Interpolant*> dndxInterpolant1d) {
	dndx_interpolant_1d_ = dndxInterpolant1d;
}

void Epairproduction::SetDndxInterpolant2d(
		std::vector<Interpolant*> dndxInterpolant2d) {
	dndx_interpolant_2d_ = dndxInterpolant2d;
}

void Epairproduction::SetLpm(double lpm) {
	eLpm_ = lpm;
}

void Epairproduction::SetIntegral(Integral* integral) {
	integral_ = integral;
}

void Epairproduction::SetIntegralForDEdx(Integral* integralForDEdx) {
	integral_for_dEdx_ = integralForDEdx;
}

void Epairproduction::SetProbForComponent(
		std::vector<double> probForComponent) {
	prob_for_component_ = probForComponent;
}

void Epairproduction::SetReverse(bool reverse) {
	reverse_ = reverse;
}

void Epairproduction::SetV(double v) {
	v_ = v;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Destructor---------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Epairproduction::~Epairproduction()
{
    delete integral_for_dEdx_;
    for(unsigned int i = 0 ; i < dndx_integral_.size() ; i++ ){
        delete dndx_integral_[i];
    }

    dndx_integral_.clear();
}


