#pragma once

#include <string>
#include <memory>

#include "TFile.h"
#include "TTree.h"

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/Secondaries.h"
#include "PROPOSAL/math/Vector3D.h"

namespace std {
template<typename T, typename... Args>
unique_ptr<T> make_unique(Args&&... args)
{
    return unique_ptr<T>(new T(forward<Args>(args)...));
}
}

namespace PROPOSAL {
class Root {
    public:
        Root(std::string);
        ~Root() {};
        // ~Root() {delete rootfile_; delete secondaries_tree_;};

        void StoreDynamicData(const DynamicData& primary);
        // void StoreVector(const Vector3D& vec);
        // void StoreSecondaries(const Secondaries& secondaries);
        void Close();

        TFile* rootfile_;
        TTree* secondaries_tree_;

        // Event* vector_tree_;
        // std::unique_ptr<TFile> rootfile_;
        // std::unique_ptr<TTree> secondaries_tree_;

    private:
        double primary_position_ ;
        double primary_t_;
        // Vector3D primary_direction_;
        double primary_energy_;
        double primary_parent_particle_energy;
        double primary_propagated_distance;
        char primary_interaction_type[32];

        Vector3D::CartesianCoordinates position_;
};
}
