
#include "PROPOSAL/Output.h"
#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

std::vector<DynamicData*> Output::secondarys_;
bool Output::store_in_root_trees_ = false;
bool Output::store_in_ASCII_file_ = false;


// ------------------------------------------------------------------------- //
void Output::FillSecondaryVector(const std::vector<Particle*>& particles)
{
    // Do not copy values. This method is used to store particles from the ouput of
    // the decay channel. The decay channel creates new particle. So the output becomes the owner
    secondarys_.insert(secondarys_.end(), particles.begin(), particles.end());
}

// ------------------------------------------------------------------------- //
void Output::FillSecondaryVector(const Particle& particle, const DynamicData::Type& secondary, double energyloss)
{
    DynamicData* data = NULL;

    if (secondary == DynamicData::Particle)
    {
        data = new Particle(particle);
        log_fatal("Returning a particle %s\n", particle.GetName().c_str());
    } else
    {
        data = new DynamicData(secondary);

        data->SetEnergy(energyloss);
        data->SetPosition(particle.GetPosition());
        data->SetDirection(particle.GetDirection());
        // TODO(mario): dedcide to have an id Mon 2017/09/11
        // data->SetParticleId(particle->GetParticleId() + 1);
        // data->SetParentParticleId(particle->GetParentParticleId());
        data->SetTime(particle.GetTime());
        data->SetParentParticleEnergy(particle.GetEnergy());
        data->SetPropagatedDistance(particle.GetPropagatedDistance());
    }

    secondarys_.push_back(data);
}

// ------------------------------------------------------------------------- //
void Output::FillSecondaryVector(DynamicData* continuous_loss)
{
    // same like the decay output but for the continuous energy losses
    secondarys_.push_back(continuous_loss);
}
// ------------------------------------------------------------------------- //
void Output::ClearSecondaryVector()
{
    for (unsigned int i = 0; i < secondarys_.size(); i++)
    {
        delete secondarys_[i];
    }
    secondarys_.clear();
}

// ------------------------------------------------------------------------- //
void Output::Close()
{
#if ROOT_SUPPORT
    if (store_in_root_trees_)
    {
        secondary_tree_->Write();
        primary_tree_->Write();
        propagated_primary_tree_->Write();
        rootfile_->Close();
    }
#endif
    if (store_in_ASCII_file_)
    {
        secondary_ascii_.close();
        primary_ascii_.close();
        propagated_primary_ascii_.close();
    }
}

// ------------------------------------------------------------------------- //
void Output::EnableROOTOutput(std::string rootfile_name)
{
#if ROOT_SUPPORT
    rootfile_ = new TFile(rootfile_name.c_str(), "RECREATE");

    secondary_tree_          = new TTree("secondarys", "secondarys");
    primary_tree_            = new TTree("primarys", "primarys");
    propagated_primary_tree_ = new TTree("propagated_primarys", "propagated_primarys");

    store_in_root_trees_ = true;

    // Set the Branch Addresses;
    secondary_tree_->Branch("position", &secondary_position_, "position/Vector3D");
    secondary_tree_->Branch("t", &secondary_t_, "t/D");
    secondary_tree_->Branch("direction", &secondary_direction_, "direction/Vector3D");
    secondary_tree_->Branch("energy", &secondary_energy_, "energy/D");
    secondary_tree_->Branch("parent_particle_id", &secondary_parent_particle_id_, "parent_particle_id/I");
    secondary_tree_->Branch("particle_id", &secondary_particle_id_, "particle_id/I");
    secondary_tree_->Branch("name", &secondary_name_, "name/C");
    secondary_tree_->Branch("current_primary_energy", &current_primary_energy_);

    primary_tree_->Branch("position", &primary_position_, "position/Vector3D");
    primary_tree_->Branch("t", &primary_t_, "t/D");
    primary_tree_->Branch("direction", &primary_direction_, "direction/Vector3D");
    primary_tree_->Branch("energy", &primary_energy_, "energy/D");
    primary_tree_->Branch("parent_particle_id", &primary_parent_particle_id_, "parent_particle_id/I");
    primary_tree_->Branch("particle_id", &primary_particle_id_, "particle_id/I");
    primary_tree_->Branch("name", &primary_name_, "name/C");

    propagated_primary_tree_->Branch("position", &prop_primary_position_, "position/Vector3D");
    propagated_primary_tree_->Branch("t", &prop_primary_t_, "t/D");
    propagated_primary_tree_->Branch("direction", &prop_primary_direction_, "direction/Vector3D");
    propagated_primary_tree_->Branch("energy", &prop_primary_energy_, "energy/D");
    propagated_primary_tree_->Branch("parent_particle_id", &prop_primary_parent_particle_id_, "parent_particle_id/I");
    propagated_primary_tree_->Branch("particle_id", &prop_primary_particle_id_, "particle_id/I");
    propagated_primary_tree_->Branch("name", &prop_primary_name_, "name/C");
    propagated_primary_tree_->Branch("entry_point", &prop_primary_entry_point_, "entry_point/Vector3D");
    propagated_primary_tree_->Branch("ti", &prop_primary_ti_, "ti/D");
    propagated_primary_tree_->Branch("ei", &prop_primary_ei_, "ei/D");
    propagated_primary_tree_->Branch("exit_point", &prop_primary_exit_point_, "exit_point/Vector3D");
    propagated_primary_tree_->Branch("tf", &prop_primary_tf_, "tf/D");
    propagated_primary_tree_->Branch("ef", &prop_primary_ef_, "ef/D");
    propagated_primary_tree_->Branch(
        "closest_approach_point", &prop_primary_closest_approach_point_, "closest_approach_point/Vector3D");
    propagated_primary_tree_->Branch("tc", &prop_primary_tc_, "tc/D");
    propagated_primary_tree_->Branch("ec", &prop_primary_ec_, "ec/D");
    propagated_primary_tree_->Branch(
        "propagated_distance", &prop_primary_propagated_distance_, "propagated_distance/D");
    propagated_primary_tree_->Branch("energy_lost", &prop_primary_elost_, "energy_lost/D");
#else
    (void) rootfile_name;
    log_error("NO ROOT SUPPORT! NOTHING WILL BE STORED IN TREES!");
#endif
}

// ------------------------------------------------------------------------- //
void Output::DisableROOTOutput()
{
#if ROOT_SUPPORT
    store_in_root_trees_ = false;
#else
    log_error("NO ROOT SUPPORT! NOTHING TO CLOSE!");
#endif
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

#if ROOT_SUPPORT
void Output::StorePrimaryInTree(Particle* primary)
{
    if (store_in_root_trees_)
    {
        primary_position_           = primary->GetPosition();
        primary_t_                  = primary->GetTime();
        primary_direction_          = primary->GetDirection();
        primary_energy_             = primary->GetEnergy();
        primary_parent_particle_id_ = primary->GetParentParticleId();
        primary_particle_id_        = primary->GetParticleId();
        primary_name_               = primary->GetName();

        primary_tree_->Fill();
    }
}

// ------------------------------------------------------------------------- //
void Output::StorePropagatedPrimaryInTree(Particle* prop_primary)
{
    if (store_in_root_trees_)
    {
        prop_primary_position_               = prop_primary->GetPosition();
        prop_primary_t_                      = prop_primary->GetTime();
        prop_primary_direction_              = prop_primary->GetDirection();
        prop_primary_energy_                 = prop_primary->GetEnergy();
        prop_primary_parent_particle_id_     = prop_primary->GetParentParticleId();
        prop_primary_particle_id_            = prop_primary->GetParticleId();
        prop_primary_name_                   = prop_primary->GetName();
        prop_primary_entry_point_            = prop_primary->GetEntryPoint();
        prop_primary_ti_                     = prop_primary->GetEntryTime();
        prop_primary_ei_                     = prop_primary->GetEntryEnergy();
        prop_primary_exit_point_             = prop_primary->GetExitPoint();
        prop_primary_tf_                     = prop_primary->GetExitTime();
        prop_primary_ef_                     = prop_primary->GetExitEnergy();
        prop_primary_closest_approach_point_ = prop_primary->GetClosestApproachPoint();
        prop_primary_tc_                     = prop_primary->GetClosestApproachTime();
        prop_primary_ec_                     = prop_primary->GetClosestApproachEnergy();
        prop_primary_propagated_distance_    = prop_primary->GetPropagatedDistance();
        prop_primary_elost_                  = prop_primary->GetElost();

        propagated_primary_tree_->Fill();
    }
}

#endif

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------------ASCII OUTPUT----------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

// ------------------------------------------------------------------------- //
void Output::EnableASCIIOutput(std::string ASCII_Prefix, bool append)
{
    if (store_in_ASCII_file_)
    {
        log_info("Storing in ascii file is already enabled!");
        return;
    }

    std::stringstream ss;

    ss.str("");
    ss << ASCII_Prefix << "_secondarys.txt";
    if (append)
    {
        secondary_ascii_.open(ss.str().c_str(), std::ios::app);
    } else
    {
        secondary_ascii_.open(ss.str().c_str(), std::ios::out);
    }

    ss.str("");
    ss << ASCII_Prefix << "_primarys.txt";
    if (append)
    {
        primary_ascii_.open(ss.str().c_str(), std::ios::app);
    } else
    {
        primary_ascii_.open(ss.str().c_str(), std::ios::out);
    }

    ss.str("");
    ss << ASCII_Prefix << "_propagated_primarys.txt";
    if (append)
    {
        propagated_primary_ascii_.open(ss.str().c_str(), std::ios::app);
    } else
    {
        propagated_primary_ascii_.open(ss.str().c_str(), std::ios::out);
    }

    if (!(propagated_primary_ascii_.is_open() && secondary_ascii_.is_open() && primary_ascii_.is_open()))
    {
        log_error("Could not open ASCII streams. ASCII output disabled.");
        store_in_ASCII_file_ = false;
    }

    store_in_ASCII_file_ = true;

    WriteDescriptionFile();
}

// ------------------------------------------------------------------------- //
void Output::DisableASCIIOutput()
{
    store_in_ASCII_file_ = false;
    secondary_ascii_.close();
    primary_ascii_.close();
    propagated_primary_ascii_.close();
}

// ------------------------------------------------------------------------- //
void Output::StorePrimaryInASCII(Particle* primary)
{
    if (store_in_ASCII_file_)
    {
        primary_ascii_ << primary->GetPosition() << "\t" << primary->GetTime() << "\t" << primary->GetDirection()
                       << "\t" << primary->GetEnergy() << "\t" << primary->GetParentParticleId() << "\t"
                       << primary->GetParticleId() << "\t" << primary->GetName() << std::endl;
    }
}

// ------------------------------------------------------------------------- //
void Output::StorePropagatedPrimaryInASCII(Particle* prop_primary)
{
    if (store_in_ASCII_file_)
    {
        propagated_primary_ascii_ << prop_primary->GetPosition() << "\t" << prop_primary->GetTime() << "\t"
                                  << prop_primary->GetDirection() << "\t" << prop_primary->GetEnergy() << "\t"
                                  << prop_primary->GetParentParticleId() << "\t" << prop_primary->GetParticleId()
                                  << "\t" << prop_primary->GetName() << "\t" << prop_primary->GetEntryPoint() << "\t"
                                  << prop_primary->GetEntryTime() << "\t" << prop_primary->GetEntryEnergy() << "\t"
                                  << prop_primary->GetExitPoint() << "\t" << prop_primary->GetExitTime() << "\t"
                                  << prop_primary->GetExitEnergy() << "\t" << prop_primary->GetClosestApproachPoint()
                                  << "\t" << prop_primary->GetClosestApproachTime() << "\t"
                                  << prop_primary->GetClosestApproachEnergy() << "\t"
                                  << prop_primary->GetPropagatedDistance() << "\t" << prop_primary->GetElost() << std::endl;
    }
}

// ------------------------------------------------------------------------- //
void Output::WriteDescriptionFile()
{
    std::ofstream description;
    description.open("ASCII_OUTPUT_DESCRIPTION.txt", std::ios_base::out);

    description << "Primary" << std::endl;
    description << "\tposition      position-coordinates of the particle in cm from detector center [cm]" << std::endl;
    description << "\tt				time-coordinate * [s]" << std::endl;
    description << "\tdirection     direction of the particle in cartesian [cm] and spherical [rad] coordinates"
                << std::endl;
    description << "\tenergy			particle energy [MeV]" << std::endl;
    description << "\tparentParticleId		Id of the parent particle (in real simulation e.g. HE4 Nucleus)"
                << std::endl;
    description << "\tparticleId			Id of the particle" << std::endl;
    description << "\tname			name of the particle" << std::endl;

    description << std::endl;

    description << "Propagated Primary" << std::endl;
    description << "\tposition      position-coordinates of the particle in cm from detector center [cm]" << std::endl;
    description << "\tt				time-coordinate * [s]" << std::endl;
    description << "\tdirection     direction of the particle in cartesian [cm] and spherical [rad] coordinates"
                << std::endl;
    description << "\tenergy			particle energy [MeV]" << std::endl;
    description << "\tparentParticleId		Id of the parent particle (in real simulation e.g. HE4 Nucleus)"
                << std::endl;
    description << "\tparticleId			Id of the particle" << std::endl;
    description << "\tname			name of the particle" << std::endl;
    description << "\tentry_point       coordinates of the particle entering the detector center [cm]" << std::endl;
    description << "\tti				time-coordinate * [s]" << std::endl;
    description << "\tenergyi			particle energy * [MeV]" << std::endl;
    description << "\texit_point        coordinate of the particle exiting the detector center [cm]" << std::endl;
    description << "\ttf				time-coordinate * [s]" << std::endl;
    description << "\tenergyf			particle energy * [MeV]" << std::endl;
    description << "\tclosest_approach_point    coordinate of the particles closes approach to the detector center [cm]"
                << std::endl;
    description << "\ttc				time-coordinate * [s]" << std::endl;
    description << "\tenergyc			particle energy * [MeV]" << std::endl;
    description << "\tpropagated distance	propagated distance of the particle [cm]" << std::endl;
    description << "\tElost				particle energy which was lost in the detector [MeV]" << std::endl;

    description << "Secondary" << std::endl;
    description << "\tposition      position-coordinates of the secondary in cm from detector center [cm]" << std::endl;
    description << "\tt				time-coordinate * [s]" << std::endl;
    description << "\tdirection     direction of the particle in cartesian [cm] and spherical [rad] coordinates"
                << std::endl;
    description << "\tenergy			secondary energy [MeV]" << std::endl;
    description << "\tparentParticleId		Id of the parent particle (in real simulation e.g. HE4 Nucleus)"
                << std::endl;
    description << "\tparticleId			Id of the particle (here propagated particle ID +1)" << std::endl;
    description << "\tname			name of the secondary" << std::endl;
    description << "\tParentEnergy		energy of the particle which created the secondary [MeV]" << std::endl;

    description.close();
}
