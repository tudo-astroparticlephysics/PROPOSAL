#include <string>
#include <iostream>
#include <cstring>
#include <memory>

#include "TFile.h"
#include "TTree.h"

#include "PROPOSAL/interfaces/root.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/Secondaries.h"


using namespace PROPOSAL;

Root::Root(std::string rootfile_name) 
{
    rootfile_ = new TFile(rootfile_name.c_str(), "RECREATE");
    secondaries_tree_ = new TTree("primarys", "primarys");


    // Set the Branch Addresses;
    secondaries_tree_->Branch("position", &position_, "x/D:y/D:z/D" );
    secondaries_tree_->Branch("t", &primary_t_, "t/D");
    secondaries_tree_->Branch("energy", &primary_energy_, "energy/D");
    secondaries_tree_->Branch("parent_particle_energy", &primary_parent_particle_energy, "parent_particle_energy/D");
    secondaries_tree_->Branch("propagated_distance", &primary_propagated_distance, "propagated_distance/D");
    secondaries_tree_->Branch("interaction", &primary_interaction_type, "interaction/C");
    
}

void Root::StoreDynamicData(const DynamicData& primary)
{

    position_                      = primary.GetPosition().GetCartesianCoordinates();
    primary_t_                     = primary.GetTime();
    primary_energy_                = primary.GetEnergy();
    primary_parent_particle_energy = primary.GetParentParticleEnergy();
    primary_propagated_distance    = primary.GetPropagatedDistance();

    std::cout << primary.GetTypeId() <<std::endl;
    std::cout << primary.GetNameFromType(primary.GetTypeId()) <<std::endl;

    std::strcpy(primary_interaction_type, primary.GetNameFromType(primary.GetTypeId()).c_str());

    secondaries_tree_->Fill();
}


void Root::Close()
{
    secondaries_tree_->Write();
    rootfile_->Close();
}
