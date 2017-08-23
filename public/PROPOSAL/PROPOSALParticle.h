/*! \file   PROPOSALParticle.h
*   \brief  Header file for the Particle routines.
*
*   For more details see the class documentation.
*
*   \date   2013.03.14
*   \author Jan-Hendrik Köhne
*/
#pragma once

#include <string>
#include <vector>

#include "PROPOSAL/ParticleDef.h"
#include "PROPOSAL/Vector3D.h"


namespace PROPOSAL{

// ----------------------------------------------------------------------------
/// @brief Particle enums
// ----------------------------------------------------------------------------

namespace ParticleType
{
    enum Enum
    {
        // The numbering is after the Monte Carlo Particle Numbering Scheme
        // of the Particle Data Group
        // Chin. Phys. C, 40, 100001 (2016)

        unknown = 0,

        // Leptons
        EMinus   =  11,
        EPlus    = -11,
        NuE      =  12,
        NuEBar   = -12,
        MuMinus  =  13,
        MuPlus   = -13,
        NuMu     =  14,
        NuMuBar  = -14,
        TauMinus =  15,
        TauPlus  = -15,
        NuTau    =  16,
        NuTauBar = -16,

        // Gauge Bosons
        Gamma  =  22,
        Z0     =  23,
        WPlus  =  24,
        WMinus = -24,

        // Mesons
        Pi0     =  111,
        PiPlus  =  211,
        PiMinus = -211,
        KPlus   =  321,
        KMinus  = -321,
        // K0_Long = 130,
        // K0_Short = 310,
        // K0 = 311
        // Eta = 221,
        // DPlus = 411,
        // DMinus = -411,
        // D0 = 421,
        // D0Bar = -421,
        // DsPlus = 431,
        // DsMinusBar = -431,

        // Baryons
        PPlus   =  2212,
        PMinus  = -2212,
        Neutron =  2112,
        // NeutronBar = -2112,
        // Lambda = 3122,
        // LambdaBar = -3122,
        // SigmaPlus = 3222,
        // SigmaPlusBar = -3112,
        // Sigma0 = 3212,
        // Sigma0Bar = -3212,
        // SigmaMinus = 3112,
        // SigmaMinusBar = -3222,
        // Xi0 = 3322,
        // Xi0Bar = -3322,
        // XiMinus = 3312,
        // XiPlusBar = -3312,
        // OmegaMinus = 3334,
        // OmegaPlusBar = -3334,
        // LambdacPlus = 4122,

        // Cross section types
        Brems = -1001,
        DeltaE = -1002,
        EPair = -1003,
        NuclInt = -1004,
        MuPair = -1005,
        Hadrons = -1006,
        ContinuousEnergyLoss = -1111,

        // Exotic particles
        Monopole  = -41,
        STauMinus =  1000015,
        STauPlus  = -1000015,
        StableMassiveParticle   =  1000016
    };
}


// ----------------------------------------------------------------------------
/// @brief This class provides the main particle properties and functions.
///
/// All coordinates, angles and physical values are stored in this class.
// ----------------------------------------------------------------------------
class PROPOSALParticle
{
    public:

    PROPOSALParticle();
    PROPOSALParticle(ParticleDef);

    // destructors
    virtual ~PROPOSALParticle() {}

    void swap(PROPOSALParticle& particle);

    // Operators
    PROPOSALParticle& operator=(const PROPOSALParticle&);
    bool operator==(const PROPOSALParticle& particle) const;
    bool operator!=(const PROPOSALParticle& particle) const;
    friend std::ostream& operator<<(std::ostream& os, PROPOSALParticle const& particle);

    // --------------------------------------------------------------------- //
    // Getter & Setter
    // --------------------------------------------------------------------- //

    // Setter
    void SetPropagatedDistance(double prop_dist);
    void SetMomentum(double momentum);
    void SetEnergy(double energy);
    void SetLow(double low);

    void SetParentParticleId(int parent_particle_id);
    void SetParentParticleEnergy(double parent_particle_energy);
    void SetParticleId(int particle_id);

    void SetPosition(const Vector3D& position);
    void SetT(double t);

    void SetDirection(const Vector3D& direction);

    void SetEntryPoint(Vector3D& entry_point);
    void SetTi(double ti);
    void SetEi(double ei);

    void SetExitPoint(Vector3D& exit_point);
    void SetTf(double tf);
    void SetEf(double ef);

    void SetClosestApproachPoint(Vector3D& closest_approach_point);
    void SetTc(double tc);
    void SetEc(double ec);

    void SetElost(double elost);

    // Getter
    const ParticleDef& GetParticleDef() const { return *particle_def_; }
    const DecayTable& GetDecayTable() const { return particle_def_->decay_table; }
    double GetPropagatedDistance() const { return propagated_distance_; }

    double GetT() const { return t_; }
    double GetMomentum() const { return momentum_; }
    double GetEnergy() const { return energy_; }
    double GetLow() const { return low_; }

    double GetMass() const { return particle_def_->mass; }
    double GetLifetime() const { return particle_def_->lifetime; }
    double GetCharge() const { return particle_def_->charge; }
    std::string GetName() const { return particle_def_->name; }

    Vector3D GetPosition() const { return position_; }
    Vector3D GetDirection() const { return direction_; }

    // ----------------------------------------------------------------------------
    /// @brief Return interpolation tables used for Photonuclear Crossection
    ///
    /// Only for muons and taus a non empty vector will be returned.
    ///
    /// @return 2d double vector
    // ----------------------------------------------------------------------------
    const HardBBTables::VecType* getHardBB() {return particle_def_->hardbb_table;}

    double GetParentParticleEnergy() const { return parent_particle_energy_; }
    int GetParentParticleId() const { return parent_particle_id_; }
    int GetParticleId() const { return particle_id_; }

    Vector3D GetEntryPoint() const { return entry_point_; }
    double GetTi() const { return ti_; }
    double GetEi() const { return ei_; }

    Vector3D GetExitPoint() const { return exit_point_; }
    double GetTf() const { return tf_; }
    double GetEf() const { return ef_; }

    Vector3D GetClosestApproachPoint() const { return closest_approach_point_; }
    double GetTc() const { return tc_; }
    double GetEc() const { return ec_; }

    double GetElost() const { return elost_; }

    protected:

    const ParticleDef* particle_def_; //!< static defenitions of the particle

    double propagated_distance_; //!< propagation distance [cm]

    double momentum_;        //!< momentum [MeV]
    double square_momentum_; //!< momentum square [MeV]
    double energy_;          //!< energy [MeV]
    double low_;       //!< energy below which the particle is lost [MeV]

    int parent_particle_id_;        //!< parent particle id
    double parent_particle_energy_; //!< energy of the parent particle
    int particle_id_;               //!< particle id

    Vector3D position_;          //!< position coordinates [cm]
    double t_;                   //!< age [sec]

    Vector3D direction_;         //!< direction vector, angles in [rad]

    Vector3D entry_point_; //!< entry point coordinates [cm]
    double ti_;            //!< t-coordinate entry Point [sec]
    double ei_;            //!< energy at entry point [MeV]

    Vector3D exit_point_; //!< exit point coordinates [cm]
    double tf_;           //!< t-coordinate exit Point [sec]
    double ef_;           //!< energy at exit point [MeV]

    Vector3D closest_approach_point_; // point of closest approach (to geometry center) [cm]
    double tc_;                       //!< t-coordinate at point of closest approach [sec]
    double ec_;                       //!< energy at at point of closest approach [MeV]

    double elost_; //!< energy lost in the detector volume [MeV]
};

} // namespace PROPOSAL
