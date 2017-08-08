/*! \file   PROPOSALParticle.h
*   \brief  Header file for the Particle routines.
*
*   For more details see the class documentation.
*
*   \date   2013.03.14
*   \author Jan-Hendrik Köhne
*/
#pragma once

#ifndef Particle_H
#define Particle_H

// #include <vector>
#include <string>

#include "PROPOSAL/Vector3D.h"

namespace PROPOSAL
{
    class PROPOSALParticle;
}

std::ostream& operator<<(std::ostream& os, PROPOSAL::PROPOSALParticle const& particle);

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


/**
  * \brief This class provides the main particle properties and functions.
  *
  * All coordinates, angles and physical values are stored in this class.
  */
class PROPOSALParticle
{
private:

    double propagated_distance_; //!< propagation distance [cm]
    Vector3D position_;         //!< position coordinates [cm]
    double t_;                   //!< age [sec]
    Vector3D direction_;        //!< direction vector, angles in [rad]

    long double costh_;          //!< cos(direction.theta)
    long double sinth_;          //!< sin(direction.theta)
    long double cosph_;          //!< cos(direction.phi)
    long double sinph_;          //!< sin(direction.phi)

    double momentum_;            //!< momentum [MeV]
    double square_momentum_;     //!< momentum square [MeV]
    double energy_;              //!< energy [MeV]
    double mass_;                //!< mass [MeV]
    double lifetime_;            //!< lifetime [sec]
    double charge_;              //!< charge

    std::string name_;           //!< name of the particle - Presetted to "mu"
    double low_;                 //!< energy below which the particle is lost [MeV]
    ParticleType::Enum type_;          //!< particle type: all particles, that can be propagated with PROPOSAL

    int parent_particle_id_;        //!< parent particle id
    double parent_particle_energy_; //!< energy of the parent particle
    int particle_id_;               //!< particle id

    Vector3D entry_point_;      //!< entry point coordinates [cm]
    double ti_;                  //!< t-coordinate entry Point [sec]
    double ei_;                  //!< energy at entry point [MeV]

    Vector3D exit_point_;       //!< exit point coordinates [cm]
    double tf_;                  //!< t-coordinate exit Point [sec]
    double ef_;                  //!< energy at exit point [MeV]

    Vector3D closest_approach_point_; // point of closest approach (to geometry center) [cm]
    double tc_;                  //!< t-coordinate at point of closest approach [sec]
    double ec_;                  //!< energy at at point of closest approach [MeV]

    double elost_;               //!< energy lost in the detector volume [MeV]


//----------------------------------------------------------------------------//

public:

    /**
     * \brief Default Constructor
     *
     * Constructor which sets "default" settings.
     */
    PROPOSALParticle();

//----------------------------------------------------------------------------//

    PROPOSALParticle(const PROPOSALParticle&);
    PROPOSALParticle& operator=(const PROPOSALParticle&);
    bool operator==(const PROPOSALParticle &particle) const;
    bool operator!=(const PROPOSALParticle &particle) const;
    friend std::ostream& operator<<(std::ostream& os, PROPOSALParticle const& particle);

//----------------------------------------------------------------------------//

    /*!
     * \brief Create particle from another particle.
     *
     * This function ist mostly used to store particle information.
     *
     * \param parent_particle_id      parent particle id
     * \param particle_id             particle id
     * \param name                    particle name
     * \param position                position-coordinates
     * \param direction               direction-coordinates (cartesian and spheric)
     * \param energy                  particle energy
     * \param t                       particle time
     * \param prop_dist               flight distance
     * \param *p                      source particle
     */

    PROPOSALParticle(int parent_particle_id,
             int particle_id,
             ParticleType::Enum type,
             Vector3D position,
             Vector3D direction,
             double energy,
             double t,
             double prop_dist,
             PROPOSALParticle *p);

//----------------------------------------------------------------------------//
    /*!
     * \brief Create particle with properties
     *
     * This function ist mostly used to store particle information.
     *
     * \param parent_particle_id      parent particle id
     * \param particle_id             particle id
     * \param name                    particle name
     * \param position                position-coordinates
     * \param direction               direction-coordinates (cartesian and spheric)
     * \param energy                  particle energy
     * \param t                       particle time
     * \param prop_dist               flight distance
     */
    PROPOSALParticle(int parent_particle_id,
             int particle_id,
             ParticleType::Enum type,
             Vector3D position,
             Vector3D direction,
             double energy,
             double t,
             double prop_dist,
             double prim_energy = 0);

//----------------------------------------------------------------------------//
    /*!
     * \brief Create particle with properties
     *
     * This function ist mostly used to store particle information.
     *
     * \param aname     particle name
     * \param position  position-coordinates
     * \param direction direction-coordinates (cartesian and spheric)
     * \param energy    particle energy
     * \param t         particle time
     */
    PROPOSALParticle(
            ParticleType::Enum type,
            Vector3D position,
            Vector3D direction,
            double energy,
            double t);

//----------------------------------------------------------------------------//
    /*!
     * \brief Create particle with properties
     *
     * This constructor is used to create a Particle only by name.
     * For example whne a ProcessCollection is created no other informations as
     * mass, lifetime etc. are needed and no infos about energy or location are necessary
     *
     * \param name     particle name
     */
    PROPOSALParticle(ParticleType::Enum type);

//----------------------------------------------------------------------------//
    /*!
     * initialize particle by its name
     *
     * \param name      particle type
     */

    // void InitByName(ParticleType::Enum type);
    void InitParticle(ParticleType::Enum type);


//----------------------------------------------------------------------------//
    // Memberfunctions

    /*!
     * initialize the location and direction of the particle,\n
     * time in sec, x, y, z in cm, theta and phi in deg
     *
     * \param time      particle time
     * \param position  position-coordinates
     * \param direction direction-coordinates (cartesian and spheric)
     */

    void Location(double time,
                  Vector3D& position,
                  Vector3D& direction);

//----------------------------------------------------------------------------//

    void swap(PROPOSALParticle &particle);

//----------------------------------------------------------------------------//

    //Setter
    // void SetProperties( int parent_particle_id = 0, int particle_id = 0, double energy = 0, double t=0,
    //                     Vector3D& position = Vector3D, Vector3D& direction = Vector3D,
    //                     Vector3D& entry_point = Vector3D, double ti = 0, double Ei = 0,
    //                     Vector3D& exit_point = Vector3D, double tf = 0, double Ef = 0,
    //                     Vector3D& closest_approach_point = Vector3D, double tc = 0, double Ec = 0);
    void SetEnergy(double energy);
    void SetPropagatedDistance(double prop_dist);
    void SetPosition(Vector3D& position);
    void SetT(double t);
    void SetDirection(Vector3D& direction);
    void SetMomentum(double momentum);
    void SetMass(double mass);
    void SetLifetime(double lifetime);
    void SetCharge(double charge);
    void SetName(std::string name);
    void SetLow(double low);
    void SetType(ParticleType::Enum type);
    void SetParentParticleId(int parent_particle_id);
    void SetParentParticleEnergy(double parent_particle_energy);
    void SetParticleId(int particle_id);

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


//----------------------------------------------------------------------------//
    // Getter

    double GetEnergy() const{return energy_;}
//----------------------------------------------------------------------------//
    double GetPropagatedDistance() const{return propagated_distance_;}
//----------------------------------------------------------------------------//
    Vector3D GetPosition() const{return position_;}
//----------------------------------------------------------------------------//
    double GetT() const{return t_;}
//----------------------------------------------------------------------------//
    Vector3D GetDirection() const{return direction_;}
//----------------------------------------------------------------------------//
    double GetSinTheta() const{return sinth_;}
//----------------------------------------------------------------------------//
    double GetSinPhi() const{return sinph_;}
//----------------------------------------------------------------------------//
    double GetCosTheta() const{return costh_;}
//----------------------------------------------------------------------------//
    double GetCosPhi() const{return cosph_;}
//----------------------------------------------------------------------------//
    double GetMomentum() const{return momentum_;}
//----------------------------------------------------------------------------//
    double GetMass() const{return mass_;}
//----------------------------------------------------------------------------//
    double GetLifetime() const{return lifetime_;}
//----------------------------------------------------------------------------//
    double GetCharge() const{return charge_;}
//----------------------------------------------------------------------------//
    std::string GetName() const{return name_;}
//----------------------------------------------------------------------------//
    static std::string GetName(ParticleType::Enum pt);
//----------------------------------------------------------------------------//
    double GetLow() const{return low_;}
//----------------------------------------------------------------------------//
    ParticleType::Enum GetType() const{return type_;}
//----------------------------------------------------------------------------//
    static ParticleType::Enum GetTypeFromName(std::string particle_name);
//----------------------------------------------------------------------------//
    int GetParentParticleId() const{return parent_particle_id_;}
//----------------------------------------------------------------------------//
    double GetParentParticleEnergy() const{return parent_particle_energy_;}
//----------------------------------------------------------------------------//
    int GetParticleId() const{return particle_id_;}
//----------------------------------------------------------------------------//
    Vector3D GetEntryPoint() const{return entry_point_;}
//----------------------------------------------------------------------------//
    double GetTi() const{return ti_;}
//----------------------------------------------------------------------------//
    double GetEi() const{return ei_;}
//----------------------------------------------------------------------------//
    Vector3D GetExitPoint() const{return exit_point_;}
//----------------------------------------------------------------------------//
    double GetTf() const{return tf_;}
//----------------------------------------------------------------------------//
    double GetEf() const{return ef_;}
//----------------------------------------------------------------------------//
    Vector3D GetClosestApproachPoint() const{return closest_approach_point_;}
//----------------------------------------------------------------------------//
    double GetTc() const{return tc_;}
//----------------------------------------------------------------------------//
    double GetEc() const{return ec_;}
//----------------------------------------------------------------------------//
    double GetElost() const{return elost_;}

//----------------------------------------------------------------------------//
    // destructors
    ~PROPOSALParticle() {}


};

}


#endif //PARTICLE_H
