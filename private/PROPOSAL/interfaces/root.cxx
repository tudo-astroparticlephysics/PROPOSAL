#include "PROPOSAL/interfaces/root.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/Secondaries.h"

#include "TFile.h"
#include "TTree.h"

#include <string>

void Root::Root(std::string rootfile_name)
{
    rootfile_ = new TFile(rootfile_name.c_str(), "RECREATE");

    /* secondary_tree_          = new TTree("secondarys", "secondarys"); */
    primary_tree_            = new TTree("primarys", "primarys");
    /* propagated_primary_tree_ = new TTree("propagated_primarys", "propagated_primarys"); */

    store_in_root_trees_ = true;

    // Set the Branch Addresses;
    /* secondary_tree_->Branch("position", &secondary_position_, "position/Vector3D"); */
    /* secondary_tree_->Branch("t", &secondary_t_, "t/D"); */
    /* secondary_tree_->Branch("direction", &secondary_direction_, "direction/Vector3D"); */
    /* secondary_tree_->Branch("energy", &secondary_energy_, "energy/D"); */
    /* secondary_tree_->Branch("parent_particle_id", &secondary_parent_particle_id_, "parent_particle_id/I"); */
    /* secondary_tree_->Branch("particle_id", &secondary_particle_id_, "particle_id/I"); */
    /* secondary_tree_->Branch("name", &secondary_name_, "name/C"); */
    /* secondary_tree_->Branch("current_primary_energy", &current_primary_energy_); */

    primary_tree_->Branch("position", &primary_position_, "position/Vector3D");
    primary_tree_->Branch("t", &primary_t_, "t/D");
    primary_tree_->Branch("direction", &primary_direction_, "direction/Vector3D");
    primary_tree_->Branch("energy", &primary_energy_, "energy/D");
    primary_tree_->Branch("parent_particle_id", &primary_parent_particle_id_, "parent_particle_id/I");
    primary_tree_->Branch("particle_id", &primary_particle_id_, "particle_id/I");
    primary_tree_->Branch("name", &primary_name_, "name/C");

    /* propagated_primary_tree_->Branch("position", &prop_primary_position_, "position/Vector3D"); */
    /* propagated_primary_tree_->Branch("t", &prop_primary_t_, "t/D"); */
    /* propagated_primary_tree_->Branch("direction", &prop_primary_direction_, "direction/Vector3D"); */
    /* propagated_primary_tree_->Branch("energy", &prop_primary_energy_, "energy/D"); */
    /* propagated_primary_tree_->Branch("parent_particle_id", &prop_primary_parent_particle_id_, "parent_particle_id/I"); */
    /* propagated_primary_tree_->Branch("particle_id", &prop_primary_particle_id_, "particle_id/I"); */
    /* propagated_primary_tree_->Branch("name", &prop_primary_name_, "name/C"); */
    /* propagated_primary_tree_->Branch("entry_point", &prop_primary_entry_point_, "entry_point/Vector3D"); */
    /* propagated_primary_tree_->Branch("ti", &prop_primary_ti_, "ti/D"); */
    /* propagated_primary_tree_->Branch("ei", &prop_primary_ei_, "ei/D"); */
    /* propagated_primary_tree_->Branch("exit_point", &prop_primary_exit_point_, "exit_point/Vector3D"); */
    /* propagated_primary_tree_->Branch("tf", &prop_primary_tf_, "tf/D"); */
    /* propagated_primary_tree_->Branch("ef", &prop_primary_ef_, "ef/D"); */
    /* propagated_primary_tree_->Branch( */
    /*     "closest_approach_point", &prop_primary_closest_approach_point_, "closest_approach_point/Vector3D"); */
    /* propagated_primary_tree_->Branch("tc", &prop_primary_tc_, "tc/D"); */
    /* propagated_primary_tree_->Branch("ec", &prop_primary_ec_, "ec/D"); */
    /* propagated_primary_tree_->Branch( */
    /*     "propagated_distance", &prop_primary_propagated_distance_, "propagated_distance/D"); */
    /* propagated_primary_tree_->Branch("energy_lost", &prop_primary_elost_, "energy_lost/D"); */
}

void Root::StoreDynamicData(const DynamicData& primary)
{
    primary_position_           = primary.GetPosition();
    primary_t_                  = primary.GetTime();
    primary_direction_          = primary.GetDirection();
    primary_energy_             = primary.GetEnergy();
    primary_parent_particle_id_ = primary.GetParentParticleId();
    primary_particle_id_        = primary.GetParticleId();
    primary_name_               = primary.GetName();

    primary_tree_->Fill();
}

void Root::StoreSecondaries(const Secondaries& secondaries)
{
    for (auto dyn_data : secondaries.GetSecondaries()) {
        primary_position_           = dyn_data.GetPosition();
        primary_t_                  = dyn_data.GetTime();
        primary_direction_          = dyn_data.GetDirection();
        primary_energy_             = dyn_data.GetEnergy();
        primary_parent_particle_id_ = dyn_data.GetParentParticleId();
        primary_particle_id_        = dyn_data.GetParticleId();
        primary_name_               = dyn_data.GetName();

        primary_tree_->Fill();
    }
}

// ------------------------------------------------------------------------- //
/* void Root::StorePropagated(Particle* prop_primary) */
/* { */
/*     prop_primary_position_               = prop_primary->GetPosition(); */
/*     prop_primary_t_                      = prop_primary->GetTime(); */
/*     prop_primary_direction_              = prop_primary->GetDirection(); */
/*     prop_primary_energy_                 = prop_primary->GetEnergy(); */
/*     prop_primary_parent_particle_id_     = prop_primary->GetParentParticleId(); */
/*     prop_primary_particle_id_            = prop_primary->GetParticleId(); */
/*     prop_primary_name_                   = prop_primary->GetName(); */
/*     prop_primary_entry_point_            = prop_primary->GetEntryPoint(); */
/*     prop_primary_ti_                     = prop_primary->GetEntryTime(); */
/*     prop_primary_ei_                     = prop_primary->GetEntryEnergy(); */
/*     prop_primary_exit_point_             = prop_primary->GetExitPoint(); */
/*     prop_primary_tf_                     = prop_primary->GetExitTime(); */
/*     prop_primary_ef_                     = prop_primary->GetExitEnergy(); */
/*     prop_primary_closest_approach_point_ = prop_primary->GetClosestApproachPoint(); */
/*     prop_primary_tc_                     = prop_primary->GetClosestApproachTime(); */
/*     prop_primary_ec_                     = prop_primary->GetClosestApproachEnergy(); */
/*     prop_primary_propagated_distance_    = prop_primary->GetPropagatedDistance(); */
/*     prop_primary_elost_                  = prop_primary->GetElost(); */

/*     propagated_primary_tree_->Fill(); */
/* } */

void Root::Close()
{
    secondary_tree_->Write();
    primary_tree_->Write();
    propagated_primary_tree_->Write();
    rootfile_->Close();
}
