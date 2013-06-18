/*
 * Propagator.cxx
 *
 *  Created on: 23.04.2013
 *      Author: koehne
 */

#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/Constants.h"
#include <cmath>
#include <utility>
#include <boost/program_options.hpp>
#include "boost/lexical_cast.hpp"

using namespace std;

namespace po	= boost::program_options;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


std::pair<double,double> Propagator::CalculateEnergyTillStochastic( double initial_energy )
{
    double rndd    =-  log(RandomDouble());
    double rndi    =-  log(RandomDouble());

    double rndiMin = 0;
    double rnddMin = 0;

    pair<double,double> final;

    //solving the tracking integral
    if(particle_->GetLifetime() < 0)
    {
        rnddMin =   0;
    }
    else
    {
        rnddMin =   collection_->CalculateTrackingIntegal(initial_energy, rndd, false)/collection_->GetDensityCorrection();
    }

    rndiMin =   collection_->CalculateTrackingIntegal(initial_energy, rndi, true);

    //evaluating the energy loss
    if(rndd >= rnddMin || rnddMin<=0)
    {
        final.second =   particle_->GetLow();
    }
    else
    {
        final.second =   collection_->CalculateFinalEnergy( initial_energy, rndd*collection_->GetDensityCorrection(), false );
    }

    if(rndi >= rndiMin || rndiMin <= 0)
    {
        final.first =   particle_->GetLow();
    }
    else
    {
        final.first =   collection_->CalculateFinalEnergy( initial_energy, rndi, true );
    }

    return final;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Propagator::Propagate( double distance )
{
    bool    flag;
    double  displacement;

    double  initial_energy  =   particle_->GetEnergy();
    double  final_energy    =   particle_->GetEnergy();

    pair<double,string> decay;
    pair<double,string> energy_loss;

    //first: final energy befor first interaction second: energy at which the
    // particle decay
    //first and second are compared to decide if interaction happens or decay
    pair<double,double> energy_till_stochastic_;


    if(distance < 0)
    {
        distance   =   0;
    }

    if(initial_energy <= particle_->GetLow() || distance==0)
    {
        flag    =   false;
    }
    else
    {
        flag    =   true;
    }

    while(flag)
    {
        energy_till_stochastic_ = CalculateEnergyTillStochastic( initial_energy );

        if(energy_till_stochastic_.first > energy_till_stochastic_.second)
        {
            particle_interaction_   =   true;
            final_energy            =   energy_till_stochastic_.first;
        }
        else
        {
            particle_interaction_   =   false;
            final_energy            =   energy_till_stochastic_.second;

        }
        cout<<"efi "<<energy_till_stochastic_.first<<"\t"<<energy_till_stochastic_.second<<"\t";

        cout<<final_energy<<"\t";

        //Calculate the displacement according to initial energy initial_energy and final_energy
        displacement  =   collection_->CalculateDisplacement(
                    initial_energy,
                    final_energy,
                    collection_->GetDensityCorrection()*(distance - particle_->GetPropagatedDistance())) / collection_->GetDensityCorrection();
        cout<<particle_->GetT()<<"\t";
        // The first interaction or decay happens behind the distance we want to propagate
        // So we calculate the final energy using only continuous losses
        if( displacement > distance - particle_->GetPropagatedDistance() )
        {
            displacement  =   distance - particle_->GetPropagatedDistance();

            final_energy  =   collection_->CalculateFinalEnergy(initial_energy, collection_->GetDensityCorrection()*displacement);

        }

        //Advance the Particle according to the displacement
        //Initial energy and final energy are needed if Molier Scattering is enabled
        AdvanceParticle(displacement, initial_energy, final_energy);

        if(abs(distance - particle_->GetPropagatedDistance()) < abs(distance)*COMPUTER_PRECISION)
        {
            particle_->SetPropagatedDistance( distance );  // computer precision control
        }

        //Randomize the continuous energy loss if this option is enabled
        if( collection_->GetDoRandomization() )
        {
            if(final_energy != particle_->GetLow())
            {
                final_energy  = collection_->Randomize( initial_energy, final_energy );
            }

        }

        // Lower limit of particle energy is reached or
        // or complete particle is propagated the whole distance
        if( final_energy == particle_->GetLow() || particle_->GetPropagatedDistance() == distance)
        {
            break;
        }

        //Set the particle energy to the current energy before making
        //stochatic losses or decay
        particle_->SetEnergy( final_energy );

        if(particle_interaction_)
        {
            energy_loss     =   collection_->MakeStochasticLoss();
            final_energy    -=  energy_loss.first;

            cout<<energy_loss.first<<"\t"<<energy_loss.second<<endl;
        }
        else
        {
            decay           =   collection_->MakeDecay();
            final_energy    =   0;

            cout<<decay.first<<"\t"<<decay.second<<endl;
        }

        //break if the lower limit of particle energy is reached
        if(final_energy <= particle_->GetLow())
        {

            break;
        }

        //Next round: update the inital energy
        initial_energy  =   final_energy;

    }

//    if(sdec)
//    {
//        if(particle_->r!=r && ef!=0 && particle_->l>=0)
//        {
//            particle_->setEnergy(particle_->m);

//            particle_->t    +=  -particle_->l*log(RandomDouble());

//            if(particle_->type==2)
//            {
//                aux =   cros->get_decay()->e(RandomDouble(), 0.5, RandomDouble(), o);
//            }
//            else
//            {
//                aux =   cros->get_decay()->e(RandomDouble(), 0.5, 0.5, o);
//            }

//            ef  =   0;

//            o->output(1, cros->get_decay()->get_out(), aux, ef);
//        }
//    }

    particle_->SetEnergy(final_energy);

    //Particle reached the border, final energy is returned
    if(particle_->GetPropagatedDistance()==distance)
    {
        return final_energy;
    }
    //The particle stopped/decayed, the propageted distance is return with a minus sign
    else
    {
        return -particle_->GetPropagatedDistance();
    }
    //Should never be here
    return 0;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Propagator::CalculateParticleTime(double ei, double ef)
{
    if(do_time_interpolation_)
    {
        if(abs(ei-ef) > abs(ei)*HALF_PRECISION)
        {
            double aux  =   interpol_time_particle_->Interpolate(ei);
            double aux2 =   aux - interpol_time_particle_->Interpolate(ef);

            if(abs(aux2) > abs(aux)*HALF_PRECISION)
            {
                return aux2;
            }
        }

        return interpol_time_particle_diff_->Interpolate( (ei+ef)/2 )*(ef-ei);
    }
    else
    {
        return time_particle_->Integrate(ei, ef, boost::bind(&Propagator::FunctionToTimeIntegral, this, _1),4);
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::AdvanceParticle(double dr, double ei, double ef)
{

    double dist = particle_->GetPropagatedDistance();
    double time = particle_->GetT();
    double x    = particle_->GetX();
    double y    = particle_->GetY();
    double z    = particle_->GetZ();

    dist   +=  dr;

    if(do_exact_time_calulation_)
    {
        time   +=  CalculateParticleTime(ei, ef)/collection_->GetDensityCorrection();
    }
    else
    {
        time   +=  dr/SPEED;
    }


//    if(propagate_->get_molieScat())
//    Implement the Molie Scattering here see PROPOSALParticle::advance of old version

//    else
//    {
    x   +=  particle_->GetSinTheta() * particle_->GetCosPhi() * dr;
    y   +=  particle_->GetSinTheta() * particle_->GetSinPhi() * dr;
    z   +=  particle_->GetCosTheta() * dr;

    particle_->SetPropagatedDistance(dist);
    particle_->SetT(time);
    particle_->SetX(x);
    particle_->SetY(y);
    particle_->SetZ(z);

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::ReadConfigFile(string config_file)
{
    bool found_detector         =   false;

    //global



    if(!FileExist(config_file))
    {
        cerr<<"Error: config file "<<config_file<<" does not exist!"<<endl;
        exit(1);
    }
    ifstream file;
    file.open(config_file.c_str());

   // int mediacur;
    char buf[256];
    string str_buf;
    deque<string> *token;
    string taux;

    while(file.good())
    {
        file.getline(buf,256);

        str_buf =   string(buf);
        token   =   SplitString(str_buf," \t");

        // Ignore empty lines
        if(token->empty())
        {
            continue;
        }

        taux =   NextToken(token);

        // Ignore lines starting with # (comments)
        if(taux.at(0)=='#')
        {
            continue;
        }

        // Reading the global options
        else if(ToLowerCase(taux).compare("global")==0)
        {
            cerr<<"Reading the global options"<<endl;
            continue;
        }
        // seed
        else if(ToLowerCase(taux).compare("seed")==0)
        {
            taux    =   NextToken(token);

            try {
                seed_ = boost::lexical_cast<int>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: The seed must be an integer! Set to 1"<<endl;
                seed_ = 1;
            }
        }
        // bremsstrahlungs parametrization
        else if(ToLowerCase(taux).compare("brems")==0)
        {
            taux    =   NextToken(token);

            try {
                brems_ = boost::lexical_cast<int>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: The bremsstrahlungs parametrization indentifier must be an integer! Set to 1"<<endl;
                brems_ = 1;
            }
        }
        // photonuclear parametrization
        else if(ToLowerCase(taux).compare("photo")==0)
        {
            taux    =   NextToken(token);

            try {
                photo_ = boost::lexical_cast<int>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: The photonuclear parametrization indentifier must be an integer! Set to 12"<<endl;
                photo_ = 12;
            }
        }
        // bremsstahlungs mulitpiler
        else if(ToLowerCase(taux).compare("brems_multiplier")==0)
        {
            taux    =   NextToken(token);

            try {
                brems_multiplier_  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: The bremsstrahlungs multiplier must be a double! Set to 1."<<endl;
                brems_multiplier_ = 1.;
            }
        }
        // photonuclear multiplier
        else if(ToLowerCase(taux).compare("photo_multiplier")==0)
        {
            taux    =   NextToken(token);

            try {
                photo_multiplier_  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: The photonuclear multiplier must be a double! Set to 1."<<endl;
                photo_multiplier_ = 1.;
            }
        }
        // epairproduction multiplier
        else if(ToLowerCase(taux).compare("epair_multiplier")==0)
        {
            taux    =   NextToken(token);

            try {
                epair_multiplier_  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: The epairproduction multiplier must be a double! Set to 1."<<endl;
                epair_multiplier_ = 1.;
            }
        }
        // ionization multiplier
        else if(ToLowerCase(taux).compare("ioniz_multiplier")==0)
        {
            taux    =   NextToken(token);

            try {
                ioniz_multiplier_  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: The ionization multiplier must be a double! Set to 1."<<endl;
                ioniz_multiplier_ = 1.;
            }
        }
        // global ecut inside the detector
        else if(ToLowerCase(taux).compare("ecut_inside")==0)
        {
            taux    =   NextToken(token);

            try {
                global_ecut_inside_  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: The ecut for inside the detector must be a double! Set to 500."<<endl;
                global_ecut_inside_ = 500.;
            }
        }
        // global ecut behind the detector
        else if(ToLowerCase(taux).compare("ecut_behind")==0)
        {
            taux    =   NextToken(token);

            try {
                global_ecut_behind_  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: The ecut for behind the detector must be a double! Set to -1."<<endl;
                global_ecut_behind_ = -1.;
            }
        }
        // global ecut infront of the detector
        else if(ToLowerCase(taux).compare("ecut_infront")==0)
        {
            taux    =   NextToken(token);

            try {
                global_ecut_infront_  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: The ecut for infront of the detector must be a double! Set to -1."<<endl;
                global_ecut_infront_ = -1.;
            }
        }
        // global vcut inside the detector
        else if(ToLowerCase(taux).compare("vcut_inside")==0)
        {
            taux    =   NextToken(token);

            try {
                global_vcut_inside_  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: The vcut for inside the detector must be a double! Set to -1."<<endl;
                global_vcut_inside_ = -1.;
            }
        }
        // global vcut behind the detector
        else if(ToLowerCase(taux).compare("vcut_behind")==0)
        {
            taux    =   NextToken(token);

            try {
                global_vcut_behind_  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: The vcut for behind the detector must be a double! Set to -1."<<endl;
                global_vcut_behind_ = -1.;
            }
        }
        // global vcut infront of the detector
        else if(ToLowerCase(taux).compare("vcut_infront")==0)
        {
            taux    =   NextToken(token);

            try {
                global_vcut_infront_  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: The vcut for infront of the detector must be a double! Set to 0.001"<<endl;
                global_vcut_infront_ = 0.001;
            }
        }
        // global continuous randominzation inside the detector
        else if(ToLowerCase(taux).compare("cont_inside")==0)
        {
            taux    =   NextToken(token);

            try {
                global_cont_inside_  = boost::lexical_cast<bool>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: cont for inside the detector must be a bool! Set to false"<<endl;
                global_cont_inside_ = false;
            }
        }
        // global continuous randominzation behind the detector
        else if(ToLowerCase(taux).compare("cont_behind")==0)
        {
            taux    =   NextToken(token);

            try {
                global_cont_behind_  = boost::lexical_cast<bool>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: cont for behind the detector must be a double! Set to false."<<endl;
                global_vcut_behind_ = false;
            }
        }
        // global continuous randominzation infront of the detector
        else if(ToLowerCase(taux).compare("cont_infront")==0)
        {
            taux    =   NextToken(token);

            try {
                global_cont_infront_  = boost::lexical_cast<bool>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: cont for infront of the detector must be a double! Set to true"<<endl;
                global_vcut_infront_ = true;
            }
        }
        // lpm effect
        else if(ToLowerCase(taux).compare("lpm")==0)
        {
            lpm_ =   true;
        }
        // moliere scattering
        else if(ToLowerCase(taux).compare("moliere")==0)
        {
            moliere_ =   true;
        }
        // exact location time
        else if(ToLowerCase(taux).compare("exact_time")==0)
        {
            do_exact_time_calulation_ =   true;
        }
        // do not interpolate: intergrate everything
        else if(ToLowerCase(taux).compare("integrate")==0)
        {
            integrate_ =   true;
        }
        // path to interpolation tables
        else if(ToLowerCase(taux).compare("path_to_tables")==0)
        {
            taux    =   NextToken(token);
            path_to_tables_ =   taux;
        }

        //Builing the detector geometry
        else if(ToLowerCase(taux).compare("detector")==0)
        {
            if(found_detector)
            {
                cerr<<"Warning: Detector already specified. There can be only one"<<endl;
                cerr<<"This one will be ignored"<<endl;
                continue;
            }
            found_detector  =   true;

            //find the first not empty line not starting with #
            while(file.good())
            {
                file.getline(buf,256);
                str_buf =   string(buf);
                token   =   SplitString(str_buf," \t");
                if(!token->empty())
                {
                    taux =   NextToken(token);
                    if(taux.at(0)!='#')
                    {
                        break;
                    }
                }
            }

            detector_   =   new Geometry();

            InitGeometry(detector_,token,taux);

        }

        // a sector consists of geometry, a medium and cut settings
        // here the ProccessCollections are initialized
        else if(ToLowerCase(taux).compare("sector")==0)
        {
            InitProcessCollections(file);
        }
        else
        {
            cerr<<"Unrecognized option: "<<taux<<endl;
            continue;
        }
    }
    cout<<photo_<<"\t"<<brems_<<"\t"<<seed_<<"\t"<<brems_multiplier_<<"\t"<<photo_multiplier_<<"\t"<<ioniz_multiplier_<<"\t"<<epair_multiplier_<<endl;
    cout<<lpm_<<"\t"<<moliere_<<"\t"<<do_exact_time_calulation_<<"\t"<<integrate_<<"\t"<<path_to_tables_<<endl;
    cout<<global_ecut_inside_<<"\t"<<global_ecut_infront_<<"\t"<<global_ecut_behind_<<endl;
    cout<<global_vcut_inside_<<"\t"<<global_vcut_infront_<<"\t"<<global_vcut_behind_<<endl;
    cout<<global_cont_inside_<<"\t"<<global_cont_infront_<<"\t"<<global_cont_behind_<<endl;


}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-----------------------Enable and disable interpolation---------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::EnableParticleTimeInterpolation()
{
    if(do_time_interpolation_)return;

    double energy = particle_->GetEnergy();

    interpol_time_particle_         =   new Interpolant(NUM3, particle_->GetLow(), BIGENERGY, boost::bind(&Propagator::InterpolTimeParticle, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);
    interpol_time_particle_diff_    =   new Interpolant(NUM3, particle_->GetLow(), BIGENERGY, boost::bind(&Propagator::InterpolTimeParticleDiff, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);

    particle_->SetEnergy(energy);
    do_time_interpolation_ =true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::DisableParticleTimeInterpolation()
{
    delete interpol_time_particle_;
    delete interpol_time_particle_diff_;

    interpol_time_particle_         = NULL;
    interpol_time_particle_diff_    = NULL;

    do_time_interpolation_ =false;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::EnableInterpolation()
{
    collection_->EnableInterpolation();
    EnableParticleTimeInterpolation();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::DisableInterpolation()
{
    collection_->DisableInterpolation();
    DisableParticleTimeInterpolation();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------Set and validate options--------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


boost::program_options::options_description Propagator::CreateOptions()
{

    po::options_description general("General options");
    general.add_options()
        ("help,h",		"shows this message")
        ("version,v",	"shows the version of the program");


    po::options_description propagator("Propagator options");
    propagator.add_options()
        ("propagator.interpol_time",        po::value<bool>(&do_time_interpolation_)->implicit_value(false),    "Enables interpolation for particle time calculation")
        ("propagator.exact_time",           po::value<bool>(&do_exact_time_calulation_)->implicit_value(false), "Do exact particle time calculations")
        ("propagator.interpol_order",       po::value<int>(&order_of_interpolation_)->default_value(5),         "number of interpolation points");


    po::options_description all("All options");
        all.add(general);
        all.add(propagator);

    for(unsigned int i =0 ; i < collection_->GetCrosssections().size(); i++)
    {
        all.add(collection_->GetCrosssections().at(i)->CreateOptions());
    }

    return all;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::Setup(int argc, char** argv)
{
    vector <string> test = po::split_unix("a b c");
    cout<<test.at(2)<<endl;
    po::options_description all = CreateOptions();

    //parse cmd line
    po::variables_map vm;
    po::store( po::command_line_parser(argc, argv).options(all).run(), vm);

    //print help message if wanted
    if(vm.count("help")) {
        std::cout<< all;
        exit(0);
    }
    //notifies globalVar
    try {
        //set the variables
        vm.notify();
    }
    catch (po::invalid_command_line_syntax &e) {
        std::cerr<<"Error: "<<e.what()<<"\n";
        exit(1);
    }

    for(unsigned int i = 0 ; i < collection_->GetCrosssections().size() ; i++)
    {
        collection_->GetCrosssections().at(i)->ValidateOptions();
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Propagator::Propagator()
    :order_of_interpolation_    ( 5 )
    ,debug_                     ( false )
    ,particle_interaction_      ( false )
    ,do_time_interpolation_     ( false )
    ,seed_                      ( 1 )
    ,brems_                     ( 1 )
    ,photo_                     ( 12 )
    ,lpm_                       ( false )
    ,moliere_                   ( false )
    ,do_exact_time_calulation_  ( false )
    ,integrate_                 ( false )
    ,brems_multiplier_          ( 1 )
    ,photo_multiplier_          ( 1 )
    ,ioniz_multiplier_          ( 1 )
    ,epair_multiplier_          ( 1 )
    ,global_ecut_inside_        ( 500 )
    ,global_ecut_infront_       ( -1 )
    ,global_ecut_behind_        ( -1 )
    ,global_vcut_inside_        ( -1 )
    ,global_vcut_infront_       ( 0.001 )
    ,global_vcut_behind_        ( -1 )
    ,global_cont_inside_        ( false )
    ,global_cont_infront_       ( true )
    ,global_cont_behind_        ( false )
    ,path_to_tables_            ( "" )
{
    particle_              = new Particle("mu",0,0,0,0,0,1e6,0);
    time_particle_         = new Integral();

    interpol_time_particle_         = NULL;
    interpol_time_particle_diff_    = NULL;

    InitDefaultCollection();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Propagator::Propagator(const Propagator &propagator)
    :order_of_interpolation_    ( propagator.order_of_interpolation_ )
    ,debug_                     ( propagator.debug_ )
    ,particle_interaction_      ( propagator.particle_interaction_ )
    ,do_time_interpolation_     ( propagator.do_time_interpolation_ )
    ,seed_                      ( propagator.seed_ )
    ,brems_                     ( propagator.brems_ )
    ,photo_                     ( propagator.photo_ )
    ,lpm_                       ( propagator.lpm_ )
    ,moliere_                   ( propagator.moliere_ )
    ,do_exact_time_calulation_  ( propagator.do_exact_time_calulation_ )
    ,integrate_                 ( propagator.integrate_ )
    ,brems_multiplier_          ( propagator.brems_multiplier_ )
    ,photo_multiplier_          ( propagator.photo_multiplier_ )
    ,ioniz_multiplier_          ( propagator.ioniz_multiplier_ )
    ,epair_multiplier_          ( propagator.epair_multiplier_ )
    ,global_ecut_inside_        ( propagator.global_ecut_inside_ )
    ,global_ecut_infront_       ( propagator.global_ecut_infront_ )
    ,global_ecut_behind_        ( propagator.global_ecut_behind_ )
    ,global_vcut_inside_        ( propagator.global_vcut_inside_ )
    ,global_vcut_infront_       ( propagator.global_vcut_infront_ )
    ,global_vcut_behind_        ( propagator.global_vcut_behind_ )
    ,global_cont_inside_        ( propagator.global_cont_inside_ )
    ,global_cont_infront_       ( propagator.global_cont_infront_ )
    ,global_cont_behind_        ( propagator.global_cont_behind_ )
    ,path_to_tables_            ( propagator.path_to_tables_ )
    ,particle_                  ( propagator.particle_ )
    ,collection_                ( new ProcessCollection(*propagator.collection_) )
    ,time_particle_             ( new Integral(*propagator.time_particle_) )

{
    if(propagator.interpol_time_particle_ != NULL)
    {
        interpol_time_particle_ = new Interpolant(*propagator.interpol_time_particle_) ;
    }
    else
    {
        interpol_time_particle_ = NULL;
    }

    if(propagator.interpol_time_particle_diff_ != NULL)
    {
        interpol_time_particle_diff_ = new Interpolant(*propagator.interpol_time_particle_diff_) ;
    }
    else
    {
        interpol_time_particle_diff_ = NULL;
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Propagator& Propagator::operator=(const Propagator &propagator)
{
    if (this != &propagator)
    {
      Propagator tmp(propagator);
      swap(tmp);
    }
    return *this;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Propagator::operator==(const Propagator &propagator) const
{
    if( order_of_interpolation_   != propagator.order_of_interpolation_ ) return false;
    if( debug_                    != propagator.debug_ )                  return false;
    if( particle_                 != propagator.particle_ )               return false;
    if( particle_interaction_     != propagator.particle_interaction_ )   return false;
    if( do_time_interpolation_    != propagator.do_time_interpolation_ )  return false;
    if( seed_                     != propagator.seed_ )                   return false;
    if( brems_                    != propagator.brems_ )                  return false;
    if( photo_                    != propagator.photo_ )                  return false;
    if( lpm_                      != propagator.lpm_ )                    return false;
    if( moliere_                  != propagator.moliere_ )                return false;
    if( do_exact_time_calulation_ != propagator.do_exact_time_calulation_ )return false;
    if( integrate_                != propagator.integrate_ )              return false;
    if( brems_multiplier_         != propagator.brems_multiplier_ )       return false;
    if( photo_multiplier_         != propagator.photo_multiplier_ )       return false;
    if( ioniz_multiplier_         != propagator.ioniz_multiplier_ )       return false;
    if( epair_multiplier_         != propagator.epair_multiplier_ )       return false;
    if( global_ecut_inside_       != propagator.global_ecut_inside_ )     return false;
    if( global_ecut_infront_      != propagator.global_ecut_infront_ )    return false;
    if( global_ecut_behind_       != propagator.global_ecut_behind_ )     return false;
    if( global_vcut_inside_       != propagator.global_vcut_inside_ )     return false;
    if( global_vcut_infront_      != propagator.global_vcut_infront_ )    return false;
    if( global_vcut_behind_       != propagator.global_vcut_behind_ )     return false;
    if( global_cont_inside_       != propagator.global_cont_inside_ )     return false;
    if( global_cont_infront_      != propagator.global_cont_infront_ )    return false;
    if( global_cont_behind_       != propagator.global_cont_behind_ )     return false;
    if( *collection_              != *propagator.collection_ )            return false;
    if( *time_particle_           != *propagator.time_particle_ )         return false;

    if( path_to_tables_.compare( propagator.path_to_tables_ )!=0 )        return false;

    if( interpol_time_particle_diff_ != NULL && propagator.interpol_time_particle_diff_ != NULL)
    {
        if( *interpol_time_particle_diff_   != *propagator.interpol_time_particle_diff_)        return false;
    }
    else if( interpol_time_particle_diff_ != propagator.interpol_time_particle_diff_)           return false;

    if( interpol_time_particle_ != NULL && propagator.interpol_time_particle_ != NULL)
    {
        if( *interpol_time_particle_   != *propagator.interpol_time_particle_)                  return false;
    }
    else if( interpol_time_particle_ != propagator.interpol_time_particle_)                     return false;

    //else
    return true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Propagator::operator!=(const Propagator &propagator) const
{
    return !(*this == propagator);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::swap(Propagator &propagator)
{
    using std::swap;

    swap( order_of_interpolation_   ,   propagator.order_of_interpolation_ );
    swap( debug_                    ,   propagator.debug_);
    swap( particle_interaction_     ,   propagator.particle_interaction_);
    swap( do_time_interpolation_    ,   propagator.do_time_interpolation_ );
    swap( seed_                     ,   propagator.seed_ );
    swap( brems_                    ,   propagator.brems_ );
    swap( photo_                    ,   propagator.photo_ );
    swap( lpm_                      ,   propagator.lpm_ );
    swap( moliere_                  ,   propagator.moliere_ );
    swap( do_exact_time_calulation_ ,   propagator.do_exact_time_calulation_ );
    swap( integrate_                ,   propagator.integrate_ );
    swap( brems_multiplier_         ,   propagator.brems_multiplier_ );
    swap( photo_multiplier_         ,   propagator.photo_multiplier_ );
    swap( ioniz_multiplier_         ,   propagator.ioniz_multiplier_ );
    swap( epair_multiplier_         ,   propagator.epair_multiplier_ );
    swap( global_ecut_inside_       ,   propagator.global_ecut_inside_ );
    swap( global_ecut_infront_      ,   propagator.global_ecut_infront_ );
    swap( global_ecut_behind_       ,   propagator.global_ecut_behind_ );
    swap( global_vcut_inside_       ,   propagator.global_vcut_inside_ );
    swap( global_vcut_infront_      ,   propagator.global_vcut_infront_ );
    swap( global_vcut_behind_       ,   propagator.global_vcut_behind_ );
    swap( global_cont_inside_       ,   propagator.global_cont_inside_ );
    swap( global_cont_infront_      ,   propagator.global_cont_infront_ );
    swap( global_cont_behind_       ,   propagator.global_cont_behind_ );

    path_to_tables_.swap( propagator.path_to_tables_ );

    particle_->swap( *propagator.particle_ );
    collection_->swap( *propagator.collection_ );
    time_particle_->swap(*propagator.time_particle_ );

    if( interpol_time_particle_ != NULL && propagator.interpol_time_particle_ != NULL)
    {
        interpol_time_particle_->swap(*propagator.interpol_time_particle_);
    }
    else if( interpol_time_particle_ == NULL && propagator.interpol_time_particle_ != NULL)
    {
        interpol_time_particle_ = new Interpolant(*propagator.interpol_time_particle_);
        propagator.interpol_time_particle_ = NULL;
    }
    else if( interpol_time_particle_ != NULL && propagator.interpol_time_particle_ == NULL)
    {
        propagator.interpol_time_particle_ = new Interpolant(*interpol_time_particle_);
        interpol_time_particle_ = NULL;
    }

    if( interpol_time_particle_diff_ != NULL && propagator.interpol_time_particle_diff_ != NULL)
    {
        interpol_time_particle_diff_->swap(*propagator.interpol_time_particle_diff_);
    }
    else if( interpol_time_particle_diff_ == NULL && propagator.interpol_time_particle_diff_ != NULL)
    {
        interpol_time_particle_diff_ = new Interpolant(*propagator.interpol_time_particle_diff_);
        propagator.interpol_time_particle_diff_ = NULL;
    }
    else if( interpol_time_particle_diff_ != NULL && propagator.interpol_time_particle_diff_ == NULL)
    {
        propagator.interpol_time_particle_diff_ = new Interpolant(*interpol_time_particle_diff_);
        interpol_time_particle_diff_ = NULL;
    }


}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------private member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::InitDefaultCollection()
{
    Medium* med             = new Medium("ice",1.);
    EnergyCutSettings* cuts = new EnergyCutSettings(500,-1);
    collection_             = new ProcessCollection(particle_ , med, cuts);

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::InitProcessCollections(ifstream &file)
{

    char buf[256];
    string str_buf;
    deque<string> *token = new deque<string>();
    string taux;

    cerr<<"Reading sector informations"<<endl;
    double ecut_inside  = -1;
    double ecut_infront = -1;
    double ecut_behind  = -1;

    double vcut_inside  = -1;
    double vcut_infront = -1;
    double vcut_behind  = -1;

    bool cont_inside    =   false;
    bool cont_infront   =   false;
    bool cont_behind    =   false;

    bool found_inside_cuts  =   false;
    bool found_behind_cuts  =   false;
    bool found_infront_cuts =   false;


    while(file.good())
    {
        file.getline(buf,256);
        str_buf =   string(buf);
        token   =   SplitString(str_buf," \t");
        if(!token->empty())
        {
            taux =   NextToken(token);
            if(taux.at(0)!='#')
            {
                break;
            }
        }
    }
    Geometry *geometry  = new Geometry();
    InitGeometry(geometry,token,taux);

    while(file.good())
    {
        file.getline(buf,256);
        str_buf =   string(buf);
        token   =   SplitString(str_buf," \t");
        if(token->empty())
        {
            continue;
        }
        taux =   NextToken(token);

        if(taux.at(0)=='#')
        {
            continue;
        }

        if(ToLowerCase(taux).compare("inside")==0)
        {
            if(token->size()==3)
            {
                taux    =   NextToken(token);
                try {
                    ecut_inside  = boost::lexical_cast<double>(taux);
                }
                catch(boost::bad_lexical_cast&) {
                    cout<<"Warning: ecut for inside of the detector must be a double! Set to 500."<<endl;
                    ecut_inside = 500;
                }

                taux    =   NextToken(token);
                try {
                    vcut_inside  = boost::lexical_cast<double>(taux);
                }
                catch(boost::bad_lexical_cast&) {
                    cout<<"Warning: vcut for inside of the detector must be a double! Set to -1."<<endl;
                    vcut_inside = -1;
                }

                taux    =   NextToken(token);
                try {
                    cont_inside  = boost::lexical_cast<bool>(taux);
                }
                catch(boost::bad_lexical_cast&) {
                    cout<<"Warning: cont for inside of the detector must be a bool! Set to false."<<endl;
                    cont_inside = false;
                }

                found_inside_cuts = true;

            }
            else
            {
                cerr<<"Warning: Expect 3 parameters afer keyword inside! Set inside cut settings to global cut settings"<<endl;
            }
        }
        else if(ToLowerCase(taux).compare("infront")==0)
        {
            if(token->size()==3)
            {
                taux    =   NextToken(token);
                try {
                    ecut_infront  = boost::lexical_cast<double>(taux);
                }
                catch(boost::bad_lexical_cast&) {
                    cout<<"Warning: ecut for infront of the detector must be a double! Set to -1."<<endl;
                    ecut_infront = -1;
                }

                taux    =   NextToken(token);
                try {
                    vcut_infront  = boost::lexical_cast<double>(taux);
                }
                catch(boost::bad_lexical_cast&) {
                    cout<<"Warning: vcut for infront of the detector must be a double! Set to 0.001."<<endl;
                    vcut_infront = 0.001;
                }

                taux    =   NextToken(token);
                try {
                    cont_infront  = boost::lexical_cast<bool>(taux);
                }
                catch(boost::bad_lexical_cast&) {
                    cout<<"Warning: cont for infront of the detector must be a bool! Set to true."<<endl;
                    cont_infront = true;
                }

                found_infront_cuts = true;

            }
            else
            {
                cerr<<"Warning: Expect 3 parameters afer keyword infront! Set inside cut settings to global cut settings"<<endl;
            }
        }
        else if(ToLowerCase(taux).compare("behind")==0)
        {
            if(token->size()==3)
            {
                taux    =   NextToken(token);
                try {
                    ecut_behind  = boost::lexical_cast<double>(taux);
                }
                catch(boost::bad_lexical_cast&) {
                    cout<<"Warning: ecut for behind of the detector must be a double! Set to -1."<<endl;
                    ecut_behind = -1;
                }

                taux    =   NextToken(token);
                try {
                    vcut_behind  = boost::lexical_cast<double>(taux);
                }
                catch(boost::bad_lexical_cast&) {
                    cout<<"Warning: vcut for behind of the detector must be a double! Set to -1."<<endl;
                    vcut_behind = -1;
                }

                taux    =   NextToken(token);
                try {
                    cont_behind  = boost::lexical_cast<bool>(taux);
                }
                catch(boost::bad_lexical_cast&) {
                    cout<<"Warning: cont for behind of the detector must be a bool! Set to false."<<endl;
                    cont_behind = false;
                }
                found_behind_cuts = true;

            }
            else
            {
                cerr<<"Warning: Expect 3 parameters afer keyword behind! Set inside cut settings to global cut settings"<<endl;
            }
        }
        else if(ToLowerCase(taux).compare("medium")==0)
        {
            string name     =   NextToken(token);

            double density_correction;

            taux    =   NextToken(token);
            try {
                density_correction  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: density correction factor must be a double! Set to 1."<<endl;
                density_correction = 1;
            }

            Medium *med     =   new Medium(name,density_correction);

            Particle *mu    =   new Particle("mu");
            Particle *tau   =   new Particle("tau");
            Particle *e     =   new Particle("e");

            EnergyCutSettings *inside;
            EnergyCutSettings *infront;
            EnergyCutSettings *behind;

            if(found_inside_cuts)
            {
                inside = new EnergyCutSettings(ecut_inside,vcut_inside);
            }
            else
            {
                inside = new EnergyCutSettings(global_ecut_inside_,global_vcut_inside_);
            }

            if(found_infront_cuts)
            {
                infront = new EnergyCutSettings(ecut_infront,vcut_infront);
            }
            else
            {
                infront = new EnergyCutSettings(global_ecut_infront_,global_vcut_infront_);
            }

            if(found_inside_cuts)
            {
                behind = new EnergyCutSettings(ecut_behind,vcut_behind);
            }
            else
            {
                behind = new EnergyCutSettings(global_ecut_behind_,global_vcut_behind_);
            }

            int former_size =collections_.size();

            collections_.push_back( new ProcessCollection(new Particle(*mu),new Medium(*med),new EnergyCutSettings(*infront)) );
            collections_.push_back( new ProcessCollection(new Particle(*mu),new Medium(*med),new EnergyCutSettings(*inside)) );
            collections_.push_back( new ProcessCollection(new Particle(*mu),new Medium(*med),new EnergyCutSettings(*behind)) );

            collections_.push_back( new ProcessCollection(new Particle(*tau),new Medium(*med),new EnergyCutSettings(*infront)) );
            collections_.push_back( new ProcessCollection(new Particle(*tau),new Medium(*med),new EnergyCutSettings(*inside)) );
            collections_.push_back( new ProcessCollection(new Particle(*tau),new Medium(*med),new EnergyCutSettings(*behind)) );

            collections_.push_back( new ProcessCollection(new Particle(*e),new Medium(*med),new EnergyCutSettings(*infront)) );
            collections_.push_back( new ProcessCollection(new Particle(*e),new Medium(*med),new EnergyCutSettings(*inside)) );
            collections_.push_back( new ProcessCollection(new Particle(*e),new Medium(*med),new EnergyCutSettings(*behind)) );

            for(unsigned int i = former_size ;i<collections_.size(); i++)
            {
                collections_.at(i)->SetGeometry(geometry);
                collections_.at(i)->SetDensityCorrection(density_correction);
            }
            delete med;

            delete mu;
            delete tau;
            delete e;

            delete inside;
            delete infront;
            delete behind;
            break;
        }
        else
        {
            cerr<<"Error! Last line in a sector segment must start with ’medium’!Exit"<<endl;
            exit(1);
        }
    }
    cout<<ecut_inside<<"\t"<<ecut_infront<<"\t"<<ecut_behind<<endl;
    cout<<vcut_inside<<"\t"<<vcut_infront<<"\t"<<vcut_behind<<endl;
    cout<<cont_inside<<"\t"<<cont_infront<<"\t"<<cont_behind<<endl;
    cout<<found_inside_cuts<<"\t"<<found_behind_cuts<<"\t"<<found_infront_cuts<<endl;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::InitGeometry(Geometry* geometry, std::deque<std::string>* token, string first_token)
{
    string taux = first_token;

    double origin_x =   0;
    double origin_y =   0;
    double origin_z =   0;
    if(ToLowerCase(taux).compare("cylinder")==0)
    {
        double radius;
        double inner_radius =   0;
        double height;

        // Only radius and height are specified
        if(token->size()==2)
        {
            taux    =   NextToken(token);
            try {
                radius  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Error: radius must be a double! Exit"<<endl;
                exit(1);
            }

            taux    =   NextToken(token);
            try {
                height  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Error: height must be a double! Exit"<<endl;
                exit(1);
            }
        }
        // Only radius, inner_radius and height are specified
        else if(token->size()==3)
        {
            taux    =   NextToken(token);
            try {
                radius  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Error: radius must be a double! Exit"<<endl;
                exit(1);
            }

            taux    =   NextToken(token);
            try {
                inner_radius  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: inner_radius must be a double! Set to 0"<<endl;
                inner_radius    =   0;
            }

            taux    =   NextToken(token);
            try {
                height  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Error: height must be a double! Exit"<<endl;
                exit(1);
            }
        }
        // origin_x, origin_y, origin_z, radius, inner_radius and height are specified
        else if(token->size()==6)
        {
            taux    =   NextToken(token);
            try {
                origin_x  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: origin_x must be a double! Set to 0"<<endl;
                origin_x    =   0;
            }

            taux    =   NextToken(token);
            try {
                origin_y  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: origin_y must be a double! Set to 0"<<endl;
                origin_y    =   0;
            }

            taux    =   NextToken(token);
            try {
                origin_z  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: origin_z must be a double! Set to 0"<<endl;
                origin_z    =   0;
            }

            taux    =   NextToken(token);
            try {
                radius  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Error: radius must be a double! Exit"<<endl;
                exit(1);
            }

            taux    =   NextToken(token);
            try {
                inner_radius  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: inner_radius must be a double! Set to 0"<<endl;
                inner_radius    =   0;
            }

            taux    =   NextToken(token);
            try {
                height  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Error: height must be a double! Exit"<<endl;
                exit(1);
            }
        }
        else
        {
            cerr<<"Error: Number of values after 'cylinder' must be 2,3 or 6. Exit!"<<endl;
            exit(1);
        }

        geometry->InitCylinder(origin_x,origin_y,origin_z,radius,inner_radius,height);
        cout<<origin_x<<"\t"<<origin_y<<"\t"<<origin_z<<"\t"<<radius<<"\t"<<inner_radius<<"\t"<<height<<endl;

    }
    if(ToLowerCase(taux).compare("sphere")==0)
    {
        double radius;
        double inner_radius =   0;

        // Only radius is specified
        if(token->size()==1)
        {
            taux    =   NextToken(token);
            try {
                radius  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Error: radius must be a double! Exit"<<endl;
                exit(1);
            }
        }
        // Only radius and inner_radius are specified
        else if(token->size()==2)
        {
            taux    =   NextToken(token);
            try {
                radius  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Error: radius must be a double! Exit"<<endl;
                exit(1);
            }

            taux    =   NextToken(token);
            try {
                inner_radius  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: inner_radius must be a double! Set to 0"<<endl;
                inner_radius    =   0;
            }

        }
        // origin_x, origin_y, origin_z, radius, inner_radius and height are specified
        else if(token->size()==5)
        {
            taux    =   NextToken(token);
            try {
                origin_x  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: origin_x must be a double! Set to 0"<<endl;
                origin_x    =   0;
            }

            taux    =   NextToken(token);
            try {
                origin_y  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: origin_y must be a double! Set to 0"<<endl;
                origin_y    =   0;
            }

            taux    =   NextToken(token);
            try {
                origin_z  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: origin_z must be a double! Set to 0"<<endl;
                origin_z    =   0;
            }

            taux    =   NextToken(token);
            try {
                radius  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Error: radius must be a double! Exit"<<endl;
                exit(1);
            }

            taux    =   NextToken(token);
            try {
                inner_radius  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: inner_radius must be a double! Set to 0"<<endl;
                inner_radius    =   0;
            }

        }
        else
        {
            cerr<<"Error: Number of values after 'sphere' must be 1,2 or 5. Exit!"<<endl;
            exit(1);
        }

        geometry->InitSphere(origin_x,origin_y,origin_z,radius,inner_radius);
        cout<<origin_x<<"\t"<<origin_y<<"\t"<<origin_z<<"\t"<<radius<<"\t"<<inner_radius<<endl;

    }
    else if(ToLowerCase(taux).compare("box")==0)
    {
        double width_x;
        double width_y;
        double width_z;

        // Only width_x,width_y and width_z are specified
        if(token->size()==3)
        {
            taux    =   NextToken(token);
            try {
                width_x  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Error: width_x must be a double! Exit"<<endl;
                exit(1);
            }

            taux    =   NextToken(token);
            try {
                width_y  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Error: width_y must be a double! Exit"<<endl;
                exit(1);
            }

            taux    =   NextToken(token);
            try {
                width_z  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Error: width_z must be a double! Exit"<<endl;
                exit(1);
            }
        }
        // origin_x, origin_y, origin_z, width_x,width_y and width_z are specified
        else if(token->size()==6)
        {
            taux    =   NextToken(token);
            try {
                origin_x  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: origin_x must be a double! Set to 0"<<endl;
                origin_x    =   0;
            }

            taux    =   NextToken(token);
            try {
                origin_y  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: origin_y must be a double! Set to 0"<<endl;
                origin_y    =   0;
            }

            taux    =   NextToken(token);
            try {
                origin_z  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Warning: origin_z must be a double! Set to 0"<<endl;
                origin_z    =   0;
            }

            taux    =   NextToken(token);
            try {
                width_x  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Error: width_x must be a double! Exit"<<endl;
                exit(1);
            }

            taux    =   NextToken(token);
            try {
                width_y  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Error: width_y must be a double! Exit"<<endl;
                exit(1);
            }

            taux    =   NextToken(token);
            try {
                width_z  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                cout<<"Error: width_z must be a double! Exit"<<endl;
                exit(1);
            }

        }
        else
        {
            cerr<<"Error: Number of values after 'box' must be 3 or 6. Exit!"<<endl;
            exit(1);
        }

        geometry->InitBox(origin_x,origin_y,origin_z,width_x,width_y,width_z);
        cout<<origin_x<<"\t"<<origin_y<<"\t"<<origin_z<<"\t"<<width_x<<"\t"<<width_y<<"\t"<<width_z<<endl;

    }
    else
    {
        cerr<<"Error: Unrecognized geometry: "<<taux<<" Must be cylinder, sphere or box! Exit"<<endl;
        exit(1);
    }

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------Functions to interpolate---------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Propagator::InterpolTimeParticle(double energy)
{
    return FunctionToTimeIntegral(energy);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Propagator::InterpolTimeParticleDiff(double energy)
{
    return time_particle_->Integrate(energy, particle_->GetLow(), boost::bind(&Propagator::FunctionToTimeIntegral, this, _1),4);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------Functions to integrate----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Propagator::FunctionToTimeIntegral(double energy)
{
    double aux;

    aux     =   collection_->FunctionToIntegral(energy);
    aux     *=  particle_->GetEnergy()/(particle_->GetMomentum()*SPEED);
    return aux;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Destructor---------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Propagator::~Propagator(){}





