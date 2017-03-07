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
#include "PROPOSAL/Output.h"

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
        rnddMin =   current_collection_->CalculateTrackingIntegal(initial_energy, rndd, false)/current_collection_->GetDensityCorrection();
    }

    rndiMin =   current_collection_->CalculateTrackingIntegal(initial_energy, rndi, true);
    //evaluating the energy loss
    if(rndd >= rnddMin || rnddMin<=0)
    {
        final.second =   particle_->GetLow();
    }
    else
    {
        final.second =   current_collection_->CalculateFinalEnergy( initial_energy, rndd*current_collection_->GetDensityCorrection(), false );
    }

    if(rndi >= rndiMin || rndiMin <= 0)
    {
        final.first =   particle_->GetLow();
    }
    else
    {
        final.first =   current_collection_->CalculateFinalEnergy( initial_energy, rndi, true );
    }

    return final;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


vector<PROPOSALParticle*> Propagator::Propagate( PROPOSALParticle *particle, double MaxDistance_cm )
{
    Output::getInstance().ClearSecondaryVector();

    Output::getInstance().GetSecondarys().reserve(1000);

    #if ROOT_SUPPORT
        Output::getInstance().StorePrimaryInTree(particle);
    #endif

    if(Output::store_in_ASCII_file_)Output::getInstance().StorePrimaryInASCII(particle);

    SetParticle(particle);

    double distance_to_collection_border    =   0;
    double distance_to_detector             =   0;
    double distance_to_closest_approach     =   0;
    double distance                         =   0;
    double result                           =   0;

    // These two variables are needed to calculate the energy loss inside the detector
    // energy_at_entry_point is initialized with the current energy because this is a
    // reasonable value for particle which starts inside the detector

    double energy_at_entry_point            =   particle_->GetEnergy();
    double energy_at_exit_point             =   0;


    bool starts_in_detector     =   detector_->IsParticleInside(particle_);
    bool is_in_detector     =   false;
    bool was_in_detector    =   false;

    while(1)
    {
        ChooseCurrentCollection(particle_);
        if(current_collection_ == NULL)
        {
            log_info("particle reached the border");
            break;
        }

        // Check if have have to propagate the particle through the whole collection
        // or only to the collection border

        distance_to_collection_border =
        current_collection_->GetGeometry()->DistanceToBorder(particle_).first;
        double tmp_distance_to_border;
        for(unsigned int i = 0 ; i < collections_.size() ; i++)
        {

            if (particle->GetType() != collections_.at(i)->GetParticle()->GetType())
                continue;

            if(detector_->IsParticleInfront(particle))
            {
                if(collections_.at(i)->GetLocation() != 0)
                    continue;
                else
                {
                    if(collections_.at(i)->GetGeometry()->GetHirarchy() >= current_collection_->GetGeometry()->GetHirarchy())
                    {
                        tmp_distance_to_border = collections_.at(i)->GetGeometry()->DistanceToBorder(particle_).first;
                        if(tmp_distance_to_border<=0)continue;
                        distance_to_collection_border = min(
                                      tmp_distance_to_border
                                    , distance_to_collection_border);
                    }
                }
            }

            else if(detector_->IsParticleInside(particle))
            {
                if(collections_.at(i)->GetLocation() != 1)
                    continue;
                else
                {
                    tmp_distance_to_border = collections_.at(i)->GetGeometry()->DistanceToBorder(particle_).first;
                    if(tmp_distance_to_border<=0)continue;
                    distance_to_collection_border = min(
                                  tmp_distance_to_border
                                , distance_to_collection_border);
                }

            }

            else if(detector_->IsParticleBehind(particle))
            {
                if(collections_.at(i)->GetLocation() != 2)
                    continue;
                else
                {
                    if(collections_.at(i)->GetGeometry()->GetHirarchy() >= current_collection_->GetGeometry()->GetHirarchy())
                    {
                        tmp_distance_to_border = collections_.at(i)->GetGeometry()->DistanceToBorder(particle_).first;
                        if(tmp_distance_to_border<=0)continue;
                        distance_to_collection_border = min(
                                      tmp_distance_to_border
                                    , distance_to_collection_border);
                    }
                    //The particle reached the border of all specified collections
                    else
                    {

                    }
                }
            }
        }

        distance_to_detector =
        detector_->DistanceToBorder(particle_).first;

        distance_to_closest_approach  =
        detector_->DistanceToClosestApproach(particle_);

        if(abs(distance_to_closest_approach) < GEOMETRY_PRECISION )
        {
            particle_->SetXc( particle_->GetX() );
            particle_->SetYc( particle_->GetY() );
            particle_->SetZc( particle_->GetZ() );
            particle_->SetEc( particle_->GetEnergy() );
            particle_->SetTc( particle_->GetT() );

            distance_to_closest_approach    =   0;

        }

        if(distance_to_detector > 0)
        {
            if(distance_to_closest_approach > 0)
            {
                if( distance_to_detector < distance_to_collection_border &&
                    distance_to_detector < distance_to_closest_approach )
                {
                    distance    =   distance_to_detector;
                }
                else if( distance_to_closest_approach < distance_to_collection_border)
                {
                    distance    =   distance_to_closest_approach;
                }
                else
                {
                    distance    =   distance_to_collection_border;
                }
            }
            else
            {
                if( distance_to_detector < distance_to_collection_border)
                {
                    distance    =   distance_to_detector;
                }
                else
                {
                    distance    =   distance_to_collection_border;
                }
            }

        }
        else
        {
            if(distance_to_closest_approach > 0)
            {
                if( distance_to_closest_approach < distance_to_collection_border)
                {
                    distance    =   distance_to_closest_approach;
                }
                else
                {
                    distance    =   distance_to_collection_border;
                }
            }
            else
            {
                distance    =   distance_to_collection_border;
            }
        }


        is_in_detector  =   detector_->IsParticleInside(particle_);
        // entry point of the detector
        if(!starts_in_detector && !was_in_detector && is_in_detector)
        {
            particle_->SetXi( particle_->GetX() );
            particle_->SetYi( particle_->GetY() );
            particle_->SetZi( particle_->GetZ() );
            particle_->SetEi( particle_->GetEnergy() );
            particle_->SetTi( particle_->GetT() );

            energy_at_entry_point   =   particle_->GetEnergy();

            was_in_detector =   true;
        }
        // exit point of the detector
        else if(was_in_detector && !is_in_detector)
        {
            particle_->SetXf( particle_->GetX() );
            particle_->SetYf( particle_->GetY() );
            particle_->SetZf( particle_->GetZ() );
            particle_->SetEf( particle_->GetEnergy() );
            particle_->SetTf( particle_->GetT() );

            energy_at_exit_point    =   particle_->GetEnergy();
            //we don't want to run in this case a second time so we set was_in_detector to false
            was_in_detector =   false;

        }
        // if particle starts inside the detector we only ant to fill the exit point
        else if(starts_in_detector && !is_in_detector)
        {
            particle_->SetXf( particle_->GetX() );
            particle_->SetYf( particle_->GetY() );
            particle_->SetZf( particle_->GetZ() );
            particle_->SetEf( particle_->GetEnergy() );
            particle_->SetTf( particle_->GetT() );

            energy_at_exit_point    =   particle_->GetEnergy();
            //we don't want to run in this case a second time so we set starts_in_detector to false
            starts_in_detector  =   false;

        }
        if(MaxDistance_cm <= particle_->GetPropagatedDistance() + distance)
        {
            distance = MaxDistance_cm - particle_->GetPropagatedDistance();
        }
        result  =   Propagate(distance);
        if(result<=0 || MaxDistance_cm <= particle_->GetPropagatedDistance()) break;

    }

    particle_->SetElost(energy_at_entry_point - energy_at_exit_point);

    #if ROOT_SUPPORT
        Output::getInstance().StorePropagatedPrimaryInTree(particle);
    #endif
        if(Output::store_in_ASCII_file_)Output::getInstance().StorePropagatedPrimaryInASCII(particle);

    return Output::getInstance().GetSecondarys();

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Propagator::Propagate( double distance )
{
    bool    flag;
    double  displacement;

    //cout << *particle_ << endl; //Tomasz
    double propagated_distance  =   0;

    double  initial_energy  =   particle_->GetEnergy();
    double  final_energy    =   particle_->GetEnergy();

    pair<double,string> decay;
    pair<double,PROPOSALParticle::ParticleType> energy_loss;


    int secondary_id    =   0;

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

    int NumInt = 0;
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

        //Calculate the displacement according to initial energy initial_energy and final_energy
        displacement  =   current_collection_->CalculateDisplacement(
                    initial_energy,
                    final_energy,
                    current_collection_->GetDensityCorrection()*(distance - propagated_distance)) / current_collection_->GetDensityCorrection();

        // The first interaction or decay happens behind the distance we want to propagate
        // So we calculate the final energy using only continuous losses
        if( displacement > distance - propagated_distance )
        {
            displacement  =   distance - propagated_distance;

            final_energy  =   current_collection_->CalculateFinalEnergy(initial_energy, current_collection_->GetDensityCorrection()*displacement);

        }
        //Advance the Particle according to the displacement
        //Initial energy and final energy are needed if Molier Scattering is enabled
        AdvanceParticle(displacement, initial_energy, final_energy);

        propagated_distance +=  displacement;

        if(abs(distance - propagated_distance) < abs(distance)*COMPUTER_PRECISION)
        {
            propagated_distance = distance;  // computer precision control
        }
        //Randomize the continuous energy loss if this option is enabled
        if( current_collection_->GetDoRandomization() )
        {
            if(final_energy != particle_->GetLow())
            {
                final_energy  = current_collection_->Randomize( initial_energy, final_energy );
            }

        }
        // Lower limit of particle energy is reached or
        // or complete particle is propagated the whole distance
        if( final_energy == particle_->GetLow() || propagated_distance == distance)
        {
            break;
        }

        //Set the particle energy to the current energy before making
        //stochatic losses or decay
        particle_->SetEnergy( final_energy );

        if(particle_interaction_)
        {
            NumInt++;//TOMASZ
            energy_loss     =   current_collection_->MakeStochasticLoss();
            final_energy    -=  energy_loss.first;
            log_debug("Energyloss: %d\t%d\t%d\t%d\t%d", energy_loss.first, energy_loss.second->GetName(), particle_->GetX(), particle_->GetY(), particle_->GetZ());
            secondary_id    =   particle_->GetParticleId() + 1;
            Output::getInstance().FillSecondaryVector(particle_, secondary_id, energy_loss, 0);
        }
        else
        {
            decay           =   current_collection_->MakeDecay();
            final_energy    =   0;
            log_debug("Decay of particle: %s", particle_->GetName().c_str());
            secondary_id    = particle_->GetParticleId()  +   1;
            Output::getInstance().FillSecondaryVector(particle_, secondary_id, decay ,0);

        }

        //break if the lower limit of particle energy is reached
        if(final_energy <= particle_->GetLow())
        {

            break;
        }

        //Next round: update the inital energy
        initial_energy  =   final_energy;

    }

    if(stopping_decay_)
    {
        if(propagated_distance!=distance && final_energy!=0 && particle_->GetLifetime()>=0)
        {
            particle_->SetEnergy(particle_->GetMass());

            double t    =   particle_->GetT() -particle_->GetLifetime()*log(RandomDouble());
            double product_energy   =   0;

            pair<double, string> decay_to_store;
            secondary_id    =   particle_->GetParticleId() + 1;

            particle_->SetT( t );

            if(particle_->GetType()==2)
            {
                product_energy  =   current_collection_->GetDecay()->CalculateProductEnergy(RandomDouble(), 0.5, RandomDouble());
            }
            else
            {
                product_energy  =   current_collection_->GetDecay()->CalculateProductEnergy(RandomDouble(), 0.5, 0.5);
            }

            decay_to_store.first    =   product_energy;
            decay_to_store.second   =   current_collection_->GetDecay()->GetOut();

            final_energy  =   0;

            Output::getInstance().FillSecondaryVector(particle_,secondary_id, decay_to_store, final_energy);
        }
    }


    //particle_->SetParticleId(NumInt); //TOMASZ: Hack to get the number of interactions
    particle_->SetEnergy(final_energy);

    //Particle reached the border, final energy is returned
    if(propagated_distance==distance)
    {
        return final_energy;
    }
    //The particle stopped/decayed, the propageted distance is return with a minus sign
    else
    {
        return -propagated_distance;
    }
    //Should never be here
    return 0;
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
        time   +=  current_collection_->CalculateParticleTime(ei, ef)/current_collection_->GetDensityCorrection();
    }
    else
    {
        time   +=  dr/SPEED;
    }


    if(scattering_model_!=-1)
    {
        switch(scattering_model_)
        {
            case 0:
                current_collection_->GetScattering()->Scatter(dr,ei,ef);
                break;

            case 1:
                scatteringFirstOrder_->Scatter(dr ,   particle_   ,   current_collection_->GetMedium());
                break;

            case 2:
                scatteringFirstOrderMoliere_->Scatter(dr ,   particle_   ,   current_collection_->GetMedium());
                break;
            default:
                log_error("Never should be here! scattering_model = %i !",scattering_model_);
        }

    }
    else
    {
        x   +=  particle_->GetSinTheta() * particle_->GetCosPhi() * dr;
        y   +=  particle_->GetSinTheta() * particle_->GetSinPhi() * dr;
        z   +=  particle_->GetCosTheta() * dr;
        particle_->SetX(x);
        particle_->SetY(y);
        particle_->SetZ(z);
    }
    particle_->SetPropagatedDistance(dist);
    particle_->SetT(time);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::ChooseCurrentCollection(PROPOSALParticle* particle)
{
    vector<unsigned int> crossed_collections;
    crossed_collections.resize(0);

    for(unsigned int i = 0 ; i < collections_.size() ; i++)
    {
        collections_.at(i)->RestoreBackup_particle();

        if(particle->GetType() != collections_.at(i)->GetParticle()->GetType())
            continue;

        if(detector_->IsParticleInfront(particle))
        {
            if(collections_.at(i)->GetLocation() != 0)
                continue;
            else
            {
                if(collections_.at(i)->GetGeometry()->IsParticleInside(particle))
                {
                    current_collection_ = collections_.at(i);
                    crossed_collections.push_back(i);
                }
            }
        }

        else if(detector_->IsParticleInside(particle))
        {
            if(collections_.at(i)->GetLocation() != 1)
                continue;
            else
            {
                if(collections_.at(i)->GetGeometry()->IsParticleInside(particle))
                {
                    current_collection_ = collections_.at(i);
                    crossed_collections.push_back(i);
                }
            }

        }

        else if(detector_->IsParticleBehind(particle))
        {
            if(collections_.at(i)->GetLocation() != 2)
                continue;
            else
            {
                if(collections_.at(i)->GetGeometry()->IsParticleInside(particle))
                {
                    current_collection_ = collections_.at(i);
                    crossed_collections.push_back(i);
                }
                //The particle reached the border of all specified collections
                else
                {

                }
            }
        }
    }

    //No process collection was found
    if(crossed_collections.size() == 0)
    {
        current_collection_ = NULL;
        log_fatal("No Cross Section was found!!!");
    }


    //Choose current collection when multiple collections are crossed!
    //
    //Choose by hirarchy of Geometry
    //If same hirarchys are available the denser one is choosen
    //If hirarchy and density are the same then the first found is taken.
    //
    for(unsigned int i = 0  ; i < crossed_collections.size() ;i++)
    {
        unsigned int ColNow     = crossed_collections.at(i);

        //Current Hirachy is bigger -> Nothing to do!
        //
        if(current_collection_->GetGeometry()->GetHirarchy() >
                collections_.at(ColNow)->GetGeometry()->GetHirarchy() )
        {
            continue;
        }
        //Current Hirachy is equal -> Look at the density!
        //
        else if( current_collection_->GetGeometry()->GetHirarchy() ==
                 collections_.at(ColNow)->GetGeometry()->GetHirarchy() )
        {
            //Current Density is bigger or same -> Nothing to do!
            //

            if( current_collection_->GetMedium()->GetMassDensity() >=
                    collections_.at(ColNow)->GetMedium()->GetMassDensity() )
            {
                continue;
            }

            //Current Density is smaller -> Set the new collection!
            //
            else
            {
                current_collection_ =  collections_.at(ColNow);
            }

        }

        //Current Hirachy is smaller -> Set the new collection!
        //
        else
        {
            current_collection_ =  collections_.at(ColNow);
        }
    }

    if(current_collection_ != NULL)
    {
        current_collection_->SetParticle(particle_);
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::ReadConfigFile(string config_file, bool DoApplyOptions)
{
    bool found_detector         =   false;

    //global

    if(!FileExist(config_file))
    {
        log_fatal("Error: config file %s does not exist!",config_file.c_str());
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
            log_info("Reading the global options");
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
                log_warn("The seed is set to %s but must be an integer! Set to 1", taux.c_str());
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
                log_warn("The bremsstrahlungs parametrization indentifier is set to %s but must be an integer! Set to 1", taux.c_str());
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
                log_warn("The photonuclear parametrization indentifier is set to %s but must be an integer! Set to 12", taux.c_str());
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
                log_warn("The bremsstrahlungs multiplier is set to %s but must be a double! Set to 1.", taux.c_str());
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
                log_warn("The photonuclear multiplier is set to %s but must be a double! Set to 1.", taux.c_str());
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
                log_warn("The epairproduction multiplier is set to %s but must be a double! Set to 1.", taux.c_str());
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
                log_warn("The ionization multiplier is set to %s but must be a double! Set to 1.", taux.c_str());
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
                log_warn("The ecut for inside the detector is set to %s but must be a double! Set to 500.", taux.c_str());
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
                log_warn("The ecut for behind the detector is set to %s but must be a double! Set to -1.", taux.c_str());
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
                log_warn("The ecut for infront of the detector is set to %s but must be a double! Set to -1.", taux.c_str());
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
                log_warn("The vcut for inside the detector is set to %s but must be a double! Set to -1.", taux.c_str());
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
                log_warn("The vcut for behind the detector is set to %s but must be a double! Set to -1.", taux.c_str());
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
                log_warn("The vcut for infront of the detector is set to %s but must be a double! Set to 0.001", taux.c_str());
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
                log_warn("cont for inside the detector is set to %s but must be a bool! Set to false", taux.c_str());
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
                log_warn("cont for behind the detector is set to %s but must be a double! Set to false.", taux.c_str());
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
                log_warn("cont for infront of the detector is set to %s but must be a double! Set to true", taux.c_str());
                global_vcut_infront_ = true;
            }
        }
        // lpm effect
        else if(ToLowerCase(taux).compare("lpm")==0)
        {
            lpm_ =   true;
        }
        // moliere scattering
        else if(ToLowerCase(taux).compare("moliere")==0 || ToLowerCase(taux).compare("scattering")==0)
        {
            if(ToLowerCase(taux).compare("moliere")==0)
            {
                moliere_ =   true;
                scattering_model_ = 0;
                continue;
            }
            taux = NextToken(token);
            if(ToLowerCase(taux).compare("moliere")==0)
            {
                moliere_ =   true;
                scattering_model_ = 0;
                continue;
            }
            if(ToLowerCase(taux).compare("firstorder")==0)
            {
                scatteringFirstOrder_ = new ScatteringFirstOrder();
                scattering_model_ = 1;
                continue;
            }
            if(ToLowerCase(taux).compare("firstordermoliere")==0)
            {
                scatteringFirstOrderMoliere_ = new ScatteringMoliere();
                scattering_model_ = 2;
                continue;
            }
            log_error("Scattering model not known. Defaulting to moliere!");
            moliere_ =   true;
            scattering_model_ = 0;
            continue;
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
        // save interpolation tables binary or not
        else if(ToLowerCase(taux).compare("raw")==0)
        {
            raw_ =   true;
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
                log_warn("Detector already specified. There can be only one. This one will be ignored");
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
            log_warn("Unrecognized option: %s",taux.c_str());
            continue;
        }
    }

    if(DoApplyOptions)
    {
        ApplyOptions();
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-----------------------Enable and disable interpolation---------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::EnableInterpolation(std::string path, bool raw)
{
    if(current_collection_ != NULL)
    {
        current_collection_->EnableInterpolation(path,raw);
    }
    for(unsigned int i = 0 ; i < collections_.size() ; i++)
    {
        collections_.at(i)->EnableInterpolation(path,raw);
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::DisableInterpolation()
{
    if(current_collection_ != NULL)
    {
        current_collection_->DisableInterpolation();
    }
    for(unsigned int i = 0 ; i < collections_.size() ; i++)
    {
        collections_.at(i)->DisableInterpolation();
    }
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
        ("propagator.exact_time",           po::value<bool>(&do_exact_time_calulation_)->implicit_value(false), "Do exact particle time calculations")
        ("propagator.interpol_order",       po::value<int>(&order_of_interpolation_)->default_value(5),         "number of interpolation points");


    po::options_description all("All options");
        all.add(general);
        all.add(propagator);

    for(unsigned int i =0 ; i < current_collection_->GetCrosssections().size(); i++)
    {
        all.add(current_collection_->GetCrosssections().at(i)->CreateOptions());
    }

    return all;

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
    ,seed_                      ( 1 )
    ,brems_                     ( 1 )
    ,photo_                     ( 12 )
    ,lpm_                       ( false )
    ,moliere_                   ( false )
    ,stopping_decay_            ( true )
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
    ,raw_                       ( false )
    ,scattering_model_          (-1)
    ,current_collection_        (NULL)
{
    particle_              = new PROPOSALParticle(PROPOSALParticle::ParticleType::MuMinus);
    detector_              = new Geometry();
    detector_->InitSphere(0,0,0,1e18,0);
    InitDefaultCollection(detector_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Propagator::Propagator(Medium* medium,
                       EnergyCutSettings* cuts,
                       PROPOSALParticle::ParticleType particle_type,
                       string path_to_tables,
                       bool moliere,
                       bool continuous_rand,
                       bool exact_time,
                       bool lpm,
                       int brems,
                       int photo,
                       double brems_multiplier,
                       double photo_multiplier,
                       double ioniz_multiplier,
                       double epair_multiplier,
                       bool integrate,
                       int scattering_model
                       )
    :order_of_interpolation_    ( 5 )
    ,debug_                     ( false )
    ,particle_interaction_      ( false )
    ,seed_                      ( 1 )
    ,brems_                     ( brems )
    ,photo_                     ( photo )
    ,lpm_                       ( lpm )
    ,moliere_                   ( moliere )
    ,stopping_decay_            ( true )
    ,do_exact_time_calulation_  ( exact_time )
    ,integrate_                 ( integrate )
    ,brems_multiplier_          ( brems_multiplier )
    ,photo_multiplier_          ( photo_multiplier )
    ,ioniz_multiplier_          ( ioniz_multiplier )
    ,epair_multiplier_          ( epair_multiplier )
    ,global_ecut_inside_        ( 500 )
    ,global_ecut_infront_       ( -1 )
    ,global_ecut_behind_        ( -1 )
    ,global_vcut_inside_        ( -1 )
    ,global_vcut_infront_       ( 0.001 )
    ,global_vcut_behind_        ( -1 )
    ,global_cont_inside_        ( false )
    ,global_cont_infront_       ( true )
    ,global_cont_behind_        ( false )
    ,path_to_tables_            ( path_to_tables )
    ,raw_                       ( true )
    ,scattering_model_          (scattering_model)
    ,current_collection_        (NULL)
{
    particle_              = new PROPOSALParticle(particle_type);
    current_collection_    = new ProcessCollection(particle_, medium, cuts);
    detector_              = new Geometry();
    detector_->InitSphere(0,0,0,1e18,0);
    current_collection_->SetGeometry(detector_);

    for(unsigned int i =0; i<current_collection_->GetCrosssections().size(); i++)
    {
        if(current_collection_->GetCrosssections().at(i)->GetName().compare("Bremsstrahlung")==0)
        {
            current_collection_->GetCrosssections().at(i)->SetParametrization(brems_);
            current_collection_->GetCrosssections().at(i)->SetMultiplier(brems_multiplier_);
            current_collection_->GetCrosssections().at(i)->EnableLpmEffect(lpm_);

        }
        else if(current_collection_->GetCrosssections().at(i)->GetName().compare("Ionization")==0)
        {
            current_collection_->GetCrosssections().at(i)->SetMultiplier(ioniz_multiplier_);
        }
        else if(current_collection_->GetCrosssections().at(i)->GetName().compare("Epairproduction")==0)
        {
            current_collection_->GetCrosssections().at(i)->SetMultiplier(epair_multiplier_);
            current_collection_->GetCrosssections().at(i)->EnableLpmEffect(lpm_);
        }
        else if(current_collection_->GetCrosssections().at(i)->GetName().compare("Photonuclear")==0)
        {
            current_collection_->GetCrosssections().at(i)->SetParametrization(photo_);
            current_collection_->GetCrosssections().at(i)->SetMultiplier(photo_multiplier_);
        }

    }

    if(continuous_rand)
    {
        current_collection_->EnableContinuousRandomization();
    }

    if(scattering_model_ != -1)
    {
        switch(scattering_model_)
        {
            case 0:
                moliere_ = true;
                scattering_model_ = 0;
                break;

            case 1:
                moliere_ = false;
                scatteringFirstOrder_ =   new ScatteringFirstOrder();
                break;
            case 2:
                moliere_ = false;
                scatteringFirstOrderMoliere_ =   new ScatteringMoliere();
                break;
            default:
                log_error("scattering model not known! defaulting to moliere!");
                moliere = true;
                scattering_model_ = 0;
        }
    }

    if(moliere_)
    {
        current_collection_->EnableScattering();
    }

    if(do_exact_time_calulation_)
    {
        current_collection_->EnableExactTimeCalculation();
    }

    if(!integrate_)
    {
        EnableInterpolation(path_to_tables_, raw_);
    }

}
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Propagator::Propagator(string config_file, bool DoApplyOptions)
    :order_of_interpolation_    ( 5 )
    ,debug_                     ( false )
    ,particle_interaction_      ( false )
    ,seed_                      ( 1 )
    ,brems_                     ( 1 )
    ,photo_                     ( 12 )
    ,lpm_                       ( false )
    ,moliere_                   ( false )
    ,stopping_decay_            ( true )
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
    ,raw_                       ( false )
    ,scattering_model_          (-1)
    ,current_collection_        (NULL)
{
    ReadConfigFile(config_file, DoApplyOptions);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Propagator::Propagator(const Propagator &propagator)
    :order_of_interpolation_    ( propagator.order_of_interpolation_ )
    ,debug_                     ( propagator.debug_ )
    ,particle_interaction_      ( propagator.particle_interaction_ )
    ,seed_                      ( propagator.seed_ )
    ,brems_                     ( propagator.brems_ )
    ,photo_                     ( propagator.photo_ )
    ,lpm_                       ( propagator.lpm_ )
    ,moliere_                   ( propagator.moliere_ )
    ,stopping_decay_            ( propagator.stopping_decay_ )
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
    ,raw_                       ( propagator.raw_ )
    ,particle_                  ( propagator.particle_ )
    //FirstOrderScattering
    ,scatteringFirstOrder_          ( propagator.scatteringFirstOrder_ )
    ,scatteringFirstOrderMoliere_   ( propagator.scatteringFirstOrderMoliere_ )
    ,scattering_model_          (propagator.scattering_model_)
    ,current_collection_        ( new ProcessCollection(*propagator.current_collection_) )

{

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
    //FirstOrderScattering
    if( scatteringFirstOrder_           != propagator.scatteringFirstOrder_ )           return false;
    if( scatteringFirstOrderMoliere_    != propagator.scatteringFirstOrderMoliere_ )    return false;
    if( scattering_model_               != propagator.scattering_model_)                return false;

    if( particle_interaction_     != propagator.particle_interaction_ )   return false;
    if( seed_                     != propagator.seed_ )                   return false;
    if( brems_                    != propagator.brems_ )                  return false;
    if( photo_                    != propagator.photo_ )                  return false;
    if( lpm_                      != propagator.lpm_ )                    return false;
    if( moliere_                  != propagator.moliere_ )                return false;
    if( stopping_decay_           != propagator.stopping_decay_ )         return false;
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
    if( *current_collection_      != *propagator.current_collection_ )    return false;
    if( raw_                      != propagator.raw_ )                    return false;

    if( path_to_tables_.compare( propagator.path_to_tables_ )!=0 )        return false;

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
    swap( seed_                     ,   propagator.seed_ );
    swap( brems_                    ,   propagator.brems_ );
    swap( photo_                    ,   propagator.photo_ );
    swap( lpm_                      ,   propagator.lpm_ );
    swap( moliere_                  ,   propagator.moliere_ );
    swap( stopping_decay_           ,   propagator.stopping_decay_ );
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
    swap( raw_                      ,   propagator.raw_ );

    path_to_tables_.swap( propagator.path_to_tables_ );

    particle_->swap( *propagator.particle_ );
    //FirstOrderScattering
    swap<ScatteringFirstOrder*> (scatteringFirstOrder_ ,propagator.scatteringFirstOrder_);
//    scatteringFirstOrderMoliere_->swap(*propagator.scatteringFirstOrderMoliere_);
    swap(scattering_model_ , propagator.scattering_model_);

    current_collection_->swap( *propagator.current_collection_ );
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------private member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::MoveParticle(double distance)
{

    double dist = particle_->GetPropagatedDistance();

    double x    = particle_->GetX();
    double y    = particle_->GetY();
    double z    = particle_->GetZ();

    dist   +=  distance;

    x   +=  particle_->GetSinTheta() * particle_->GetCosPhi() * distance;
    y   +=  particle_->GetSinTheta() * particle_->GetSinPhi() * distance;
    z   +=  particle_->GetCosTheta() * distance;
    particle_->SetX(x);
    particle_->SetY(y);
    particle_->SetZ(z);

    particle_->SetPropagatedDistance(dist);

}


void Propagator::InitDefaultCollection(Geometry* geom)
{
    Medium* med             = new Medium("ice",1.);
    EnergyCutSettings* cuts = new EnergyCutSettings(500,0.05);
    current_collection_     = new ProcessCollection(particle_ , med, cuts);
    current_collection_->SetGeometry(geom);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::InitProcessCollections(ifstream &file)
{

    char buf[256];
    string str_buf;
    deque<string> *token = new deque<string>();
    string taux;

    unsigned int hirarchy;
    bool found_hirarchy = false;

    log_info("Reading sector informations");
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
    geometry->SetHirarchy(0);

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

        if(ToLowerCase(taux).compare("hirarchy")==0)
        {
            found_hirarchy = true;
            taux = NextToken(token);
            try
            {
                hirarchy  = boost::lexical_cast<unsigned int>(taux);
            }
            catch(boost::bad_lexical_cast&)
            {
                log_warn("Chosen hirarchy %s but must be an unsigned int! Set to 0.", taux.c_str());
                hirarchy = 0;
            }
            geometry->SetHirarchy(hirarchy);
            continue;
        }
        else if(ToLowerCase(taux).compare("inside")==0)
        {
            if(token->size()==3)
            {
                taux    =   NextToken(token);
                try {
                    ecut_inside  = boost::lexical_cast<double>(taux);
                }
                catch(boost::bad_lexical_cast&) {
                    log_warn("ecut for inside of the detector is set to %s but must be a double! Set to 500.", taux.c_str());
                    ecut_inside = 500;
                }

                taux    =   NextToken(token);
                try {
                    vcut_inside  = boost::lexical_cast<double>(taux);
                }
                catch(boost::bad_lexical_cast&) {
                    log_warn("vcut for inside of the detector is set to %s but must be a double! Set to -1.", taux.c_str());
                    vcut_inside = -1;
                }

                taux    =   NextToken(token);
                try {
                    cont_inside  = boost::lexical_cast<bool>(taux);
                }
                catch(boost::bad_lexical_cast&) {
                    log_warn("cont for inside of the detector is set to %s but must be a bool! Set to false.", taux.c_str());
                    cont_inside = false;
                }

                found_inside_cuts = true;

            }
            else
            {
                log_warn("Expect 3 parameters afer keyword inside! Set inside cut settings to global cut settings");
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
                    log_warn("ecut for infront of the detector is set to %s but must be a double! Set to -1.", taux.c_str());
                    ecut_infront = -1;
                }

                taux    =   NextToken(token);
                try {
                    vcut_infront  = boost::lexical_cast<double>(taux);
                }
                catch(boost::bad_lexical_cast&) {
                    log_warn("vcut for infront of the detector is set %s to but must be a double! Set to 0.001.", taux.c_str());
                    vcut_infront = 0.001;
                }

                taux    =   NextToken(token);
                try {
                    cont_infront  = boost::lexical_cast<bool>(taux);
                }
                catch(boost::bad_lexical_cast&) {
                    log_warn("cont for infront of the detector is set to %s but must be a bool! Set to true.", taux.c_str());
                    cont_infront = true;
                }

                found_infront_cuts = true;

            }
            else
            {
                log_warn("Expect 3 parameters afer keyword infront! Set inside cut settings to global cut settings");
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
                    log_warn("ecut for behind of the detector is set to %s but must be a double! Set to -1.", taux.c_str());
                    ecut_behind = -1;
                }

                taux    =   NextToken(token);
                try {
                    vcut_behind  = boost::lexical_cast<double>(taux);
                }
                catch(boost::bad_lexical_cast&) {
                    log_warn("vcut for behind of the detector is set to %s but must be a double! Set to -1.", taux.c_str());
                    vcut_behind = -1;
                }

                taux    =   NextToken(token);
                try {
                    cont_behind  = boost::lexical_cast<bool>(taux);
                }
                catch(boost::bad_lexical_cast&) {
                    log_warn("cont for behind of the detector is set to %s but must be a bool! Set to false.", taux.c_str());
                    cont_behind = false;
                }
                found_behind_cuts = true;

            }
            else
            {
                log_warn("Expect 3 parameters afer keyword behind! Set inside cut settings to global cut settings");
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
                log_warn("density correction factor is set to %s but must be a double! Set to 1.",taux.c_str());
                density_correction = 1;
            }

            Medium *med     =   new Medium(name,density_correction);

            PROPOSALParticle *mu    =   new PROPOSALParticle(PROPOSALParticle::ParticleType::MuMinus);
            PROPOSALParticle *tau   =   new PROPOSALParticle(PROPOSALParticle::ParticleType::TauMinus);
            PROPOSALParticle *e     =   new PROPOSALParticle(PROPOSALParticle::ParticleType::EMinus);

            EnergyCutSettings *inside;
            EnergyCutSettings *infront;
            EnergyCutSettings *behind;

            ProcessCollection* mu_inside;
            ProcessCollection* tau_inside;
            ProcessCollection* e_inside;

            ProcessCollection* mu_infront;
            ProcessCollection* tau_infront;
            ProcessCollection* e_infront;

            ProcessCollection* mu_behind;
            ProcessCollection* tau_behind;
            ProcessCollection* e_behind;

            if(found_inside_cuts)
            {
                inside = new EnergyCutSettings(ecut_inside,vcut_inside);

                mu_inside   = new ProcessCollection(new PROPOSALParticle(*mu),new Medium(*med),new EnergyCutSettings(*inside));
                tau_inside  = new ProcessCollection(new PROPOSALParticle(*tau),new Medium(*med),new EnergyCutSettings(*inside));
                e_inside    = new ProcessCollection(new PROPOSALParticle(*e),new Medium(*med),new EnergyCutSettings(*inside));

                mu_inside->SetEnableRandomization(cont_inside);
                mu_inside->SetLocation(1);

                tau_inside->SetEnableRandomization(cont_inside);
                tau_inside->SetLocation(1);

                e_inside->SetEnableRandomization(cont_inside);
                e_inside->SetLocation(1);
            }
            else
            {
                inside = new EnergyCutSettings(global_ecut_inside_,global_vcut_inside_);

                mu_inside   = new ProcessCollection(new PROPOSALParticle(*mu),new Medium(*med),new EnergyCutSettings(*inside));
                tau_inside  = new ProcessCollection(new PROPOSALParticle(*tau),new Medium(*med),new EnergyCutSettings(*inside));
                e_inside    = new ProcessCollection(new PROPOSALParticle(*e),new Medium(*med),new EnergyCutSettings(*inside));

                mu_inside->SetEnableRandomization(global_cont_inside_);
                mu_inside->SetLocation(1);

                tau_inside->SetEnableRandomization(global_cont_inside_);
                tau_inside->SetLocation(1);

                e_inside->SetEnableRandomization(global_cont_inside_);
                e_inside->SetLocation(1);
            }

            if(found_infront_cuts)
            {
                infront = new EnergyCutSettings(ecut_infront,vcut_infront);

                mu_infront   = new ProcessCollection(new PROPOSALParticle(*mu),new Medium(*med),new EnergyCutSettings(*infront));
                tau_infront  = new ProcessCollection(new PROPOSALParticle(*tau),new Medium(*med),new EnergyCutSettings(*infront));
                e_infront    = new ProcessCollection(new PROPOSALParticle(*e),new Medium(*med),new EnergyCutSettings(*infront));

                mu_infront->SetEnableRandomization(cont_infront);
                mu_infront->SetLocation(0);

                tau_infront->SetEnableRandomization(cont_infront);
                tau_infront->SetLocation(0);

                e_infront->SetEnableRandomization(cont_infront);
                e_infront->SetLocation(0);
            }
            else
            {
                infront = new EnergyCutSettings(global_ecut_infront_,global_vcut_infront_);

                mu_infront   = new ProcessCollection(new PROPOSALParticle(*mu),new Medium(*med),new EnergyCutSettings(*infront));
                tau_infront  = new ProcessCollection(new PROPOSALParticle(*tau),new Medium(*med),new EnergyCutSettings(*infront));
                e_infront    = new ProcessCollection(new PROPOSALParticle(*e),new Medium(*med),new EnergyCutSettings(*infront));

                mu_infront->SetEnableRandomization(global_cont_infront_);
                mu_infront->SetLocation(0);

                tau_infront->SetEnableRandomization(global_cont_infront_);
                tau_infront->SetLocation(0);

                e_infront->SetEnableRandomization(global_cont_infront_);
                e_infront->SetLocation(0);
            }

            if(found_behind_cuts)
            {
                behind = new EnergyCutSettings(ecut_behind,vcut_behind);

                mu_behind   = new ProcessCollection(new PROPOSALParticle(*mu),new Medium(*med),new EnergyCutSettings(*behind));
                tau_behind  = new ProcessCollection(new PROPOSALParticle(*tau),new Medium(*med),new EnergyCutSettings(*behind));
                e_behind    = new ProcessCollection(new PROPOSALParticle(*e),new Medium(*med),new EnergyCutSettings(*behind));

                mu_behind->SetEnableRandomization(cont_behind);
                mu_behind->SetLocation(2);

                tau_behind->SetEnableRandomization(cont_behind);
                tau_behind->SetLocation(2);

                e_behind->SetEnableRandomization(cont_behind);
                e_behind->SetLocation(2);
            }
            else
            {
                behind = new EnergyCutSettings(global_ecut_behind_,global_vcut_behind_);

                mu_behind   = new ProcessCollection(new PROPOSALParticle(*mu),new Medium(*med),new EnergyCutSettings(*behind));
                tau_behind  = new ProcessCollection(new PROPOSALParticle(*tau),new Medium(*med),new EnergyCutSettings(*behind));
                e_behind    = new ProcessCollection(new PROPOSALParticle(*e),new Medium(*med),new EnergyCutSettings(*behind));

                mu_behind->SetEnableRandomization(global_cont_behind_);
                mu_behind->SetLocation(2);

                tau_behind->SetEnableRandomization(global_cont_behind_);
                tau_behind->SetLocation(2);

                e_behind->SetEnableRandomization(global_cont_behind_);
                e_behind->SetLocation(2);
            }


            int former_size =collections_.size();

            collections_.push_back( mu_infront );
            collections_.push_back( mu_inside );
            collections_.push_back( mu_behind );

            collections_.push_back( tau_infront );
            collections_.push_back( tau_inside );
            collections_.push_back( tau_behind );

            collections_.push_back( e_infront );
            collections_.push_back( e_inside );
            collections_.push_back( e_behind );

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
            log_fatal("Last line in a sector segment must start with ’medium’!");
            exit(1);
        }
    }

    if(!found_hirarchy)log_info("Hirarchy for geometry not set!");
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



int Propagator::GetSeed() const
{
    return seed_;
}

void Propagator::SetSeed(int seed)
{
    seed_ = seed;
}

int Propagator::GetBrems() const
{
    return brems_;
}

void Propagator::SetBrems(int brems)
{
    brems_ = brems;
}

int Propagator::GetPhoto() const
{
    return photo_;
}

void Propagator::SetPhoto(int photo)
{
    photo_ = photo;
}

std::string Propagator::GetPath_to_tables() const
{
    return path_to_tables_;
}

void Propagator::SetPath_to_tables(const std::string &path_to_tables)
{
    path_to_tables_ = path_to_tables;
}

Geometry *Propagator::GetDetector() const
{
    return detector_;
}

void Propagator::SetDetector(Geometry *detector)
{
    detector_ = detector;
}

bool Propagator::GetStopping_decay() const
{
    return stopping_decay_;
}

void Propagator::SetStopping_decay(bool stopping_decay)
{
    stopping_decay_ = stopping_decay;
}
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
                log_fatal("radius must be a double! Exit");
                exit(1);
            }

            taux    =   NextToken(token);
            try {
                height  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                log_fatal("height must be a double! Exit");
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
                log_fatal("adius must be a double! Exit");
                exit(1);
            }

            taux    =   NextToken(token);
            try {
                inner_radius  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                log_warn("inner_radius is set to %s but must be a double! Set to 0",taux.c_str());
                inner_radius    =   0;
            }

            taux    =   NextToken(token);
            try {
                height  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                log_fatal("Error: height must be a double! Exit");
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
                log_warn("origin_x is set to %s but must be a double! Set to 0",taux.c_str());
                origin_x    =   0;
            }

            taux    =   NextToken(token);
            try {
                origin_y  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                log_warn("origin_y is set to %s but must be a double! Set to 0", taux.c_str());
                origin_y    =   0;
            }

            taux    =   NextToken(token);
            try {
                origin_z  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                log_warn("origin_z is set to %s but must be a double! Set to 0",taux.c_str());
                origin_z    =   0;
            }

            taux    =   NextToken(token);
            try {
                radius  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                log_fatal("radius must be a double! Exit");
                exit(1);
            }

            taux    =   NextToken(token);
            try {
                inner_radius  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                log_warn("inner_radius is set to %s but must be a double! Set to 0", taux.c_str());
                inner_radius    =   0;
            }

            taux    =   NextToken(token);
            try {
                height  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                log_fatal("Error: height must be a double! Exit");
                exit(1);
            }
        }
        else
        {
            log_fatal("Number of values after 'cylinder' must be 2,3 or 6. Exit!");
            exit(1);
        }

        geometry->InitCylinder(origin_x,origin_y,origin_z,radius,inner_radius,height);

    }
    else if(ToLowerCase(taux).compare("sphere")==0)
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
                log_fatal("radius must be a double! Exit");
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
                log_fatal("radius must be a double! Exit");
                exit(1);
            }

            taux    =   NextToken(token);
            try {
                inner_radius  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                log_warn("inner_radius is set to %s but must be a double! Set to 0", taux.c_str());
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
                log_warn("origin_x is set to %s but must be a double! Set to 0", taux.c_str());
                origin_x    =   0;
            }

            taux    =   NextToken(token);
            try {
                origin_y  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                log_warn("origin_y is set to %s but must be a double! Set to 0", taux.c_str());
                origin_y    =   0;
            }

            taux    =   NextToken(token);
            try {
                origin_z  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                log_warn("origin_z is set to %s but must be a double! Set to 0", taux.c_str());
                origin_z    =   0;
            }

            taux    =   NextToken(token);
            try {
                radius  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                log_fatal("Error: radius must be a double! Exit");
                exit(1);
            }

            taux    =   NextToken(token);
            try {
                inner_radius  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                log_warn("inner_radius is set to %s but must be a double! Set to 0", taux.c_str());
                inner_radius    =   0;
            }

        }
        else
        {
            log_fatal("Number of values after 'sphere' must be 1,2 or 5. Exit!");
            exit(1);
        }

        geometry->InitSphere(origin_x,origin_y,origin_z,radius,inner_radius);

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
                log_fatal("width_x must be a double! Exit");
                exit(1);
            }

            taux    =   NextToken(token);
            try {
                width_y  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                log_fatal("width_y must be a double! Exit");
                exit(1);
            }

            taux    =   NextToken(token);
            try {
                width_z  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                log_fatal("width_z must be a double! Exit");
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
                log_warn("origin_x is set to %s but must be a double! Set to 0", taux.c_str());
                origin_x    =   0;
            }

            taux    =   NextToken(token);
            try {
                origin_y  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                log_warn("origin_y is set to %s but must be a double! Set to 0", taux.c_str());
                origin_y    =   0;
            }

            taux    =   NextToken(token);
            try {
                origin_z  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                log_warn("origin_z is set to %s but must be a double! Set to 0", taux.c_str());
                origin_z    =   0;
            }

            taux    =   NextToken(token);
            try {
                width_x  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                log_fatal("width_x must be a double! Exit");
                exit(1);
            }

            taux    =   NextToken(token);
            try {
                width_y  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                log_fatal("width_y must be a double! Exit");
                exit(1);
            }

            taux    =   NextToken(token);
            try {
                width_z  = boost::lexical_cast<double>(taux);
            }
            catch(boost::bad_lexical_cast&) {
                log_fatal("width_z must be a double! Exit");
                exit(1);
            }

        }
        else
        {
            log_fatal("Number of values after 'box' must be 3 or 6. Exit!");
            exit(1);
        }

        geometry->InitBox(origin_x,origin_y,origin_z,width_x,width_y,width_z);

    }
    else
    {
        log_fatal("Unrecognized geometry: %s Must be cylinder, sphere or box! Exit",taux.c_str());
        exit(1);
    }

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::ApplyOptions()
{
    for(unsigned int j = 0 ; j < collections_.size() ; j++)
    {
        for(unsigned int i =0; i<collections_.at(j)->GetCrosssections().size(); i++)
        {
            if(collections_.at(j)->GetCrosssections().at(i)->GetName().compare("Bremsstrahlung")==0)
            {
                collections_.at(j)->GetCrosssections().at(i)->SetParametrization(brems_);
                collections_.at(j)->GetCrosssections().at(i)->SetMultiplier(brems_multiplier_);
                collections_.at(j)->GetCrosssections().at(i)->EnableLpmEffect(lpm_);

            }
            else if(collections_.at(j)->GetCrosssections().at(i)->GetName().compare("Ionization")==0)
            {
                collections_.at(j)->GetCrosssections().at(i)->SetMultiplier(ioniz_multiplier_);
            }
            else if(collections_.at(j)->GetCrosssections().at(i)->GetName().compare("Epairproduction")==0)
            {
                collections_.at(j)->GetCrosssections().at(i)->SetMultiplier(epair_multiplier_);
                collections_.at(j)->GetCrosssections().at(i)->EnableLpmEffect(lpm_);
            }
            else if(collections_.at(j)->GetCrosssections().at(i)->GetName().compare("Photonuclear")==0)
            {
                collections_.at(j)->GetCrosssections().at(i)->SetParametrization(photo_);
                collections_.at(j)->GetCrosssections().at(i)->SetMultiplier(photo_multiplier_);
            }

        }

        if(collections_.at(j)->GetEnableRandomization())
        {
            collections_.at(j)->EnableContinuousRandomization();
        }

        if(moliere_)
        {
            collections_.at(j)->EnableScattering();
        }
        if(do_exact_time_calulation_)
        {
            collections_.at(j)->EnableExactTimeCalculation();
        }

    }
    if(!integrate_)
    {
        cout << "Starting Interpolation! This will take some time depending on the number of media you defined!\n";
        cout.flush();
        EnableInterpolation(path_to_tables_, raw_);
        cout << "Done!\n";
    }

}



//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------Functions to interpolate---------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------Functions to integrate----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

void Propagator::SetParticle(PROPOSALParticle* particle)
{
    particle_   =   particle;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Destructor---------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Propagator::~Propagator(){}





