/*
 * Propagator.h
 *
 *  Created on: 2013.05.06
 *      Author: Jan-Hendrik Köhne
 */
#pragma once

/*! \mainpage PROPOSAL:
 * <b>PR</b>opagator with <b>O</b>ptimal <b>P</b>recision and <b>O</b>ptimized <b>S</b>peed for <b>A</b>ll
 * <b>L</b>eptons.
 *
 * \section intro_sec Introduction
 *
 * This is the introduction.
 *
 * \section install_sec Installation
 *
 * \subsection requirements
 *
 * \section HowTo
 *
 */

// #include <utility>
#include <boost/property_tree/ptree.hpp>
#include <deque>
#include <vector>

#include "PROPOSAL/Output.h"
#include "PROPOSAL/sector/SectorFactory.h"

namespace PROPOSAL {

// class Geometry;
// class Sector;
// class Particle;
// class Vector3D;

class Propagator
{
    public:
    // Constructors
    Propagator();
    Propagator(const std::vector<Sector*>&, const Geometry&);
    Propagator(const ParticleDef&, const std::vector<SectorFactory::Definition>&, const Geometry&);
    Propagator(const ParticleDef&, const std::vector<SectorFactory::Definition>&, const Geometry&, const InterpolationDef&);
    Propagator(const ParticleDef&, const std::string&);

    Propagator(const Propagator&);
    ~Propagator();

    Propagator& operator=(const Propagator& propagator);
    bool operator==(const Propagator& propagator) const;
    bool operator!=(const Propagator& propagator) const;

    //----------------------------------------------------------------------------//
    // Memberfunctions
    /**
     * Propagates the particle of initial energy e to the distance r.
     * Returns the final energy if the
     * particle has survived or the track length to the
     * point of disappearance with a minus sign otherwise.
     *
     *  \param  distance   maximum track length
     *  \param  energy   initial energy
     *  \return energy at distance OR -(track length)
     */

    // double Propagate(double distance);

    // ----------------------------------------------------------------------------
    /// @brief Propagates the particle through the current set of Sectors
    ///
    /// @param MaxDistance_cm
    ///
    /// @return Secondary data
    // ----------------------------------------------------------------------------
    std::vector<DynamicData*> Propagate( double MaxDistance_cm = 1e20);

    //----------------------------------------------------------------------------//
    /*!
    * advances the particle by the given distance
    * Sets the x,y and z coordinates of particle_
    * and its time and propagation distance
    *
    * \param    dr  flight distance
    * \param    ei  initial energy
    * \param    ef  final energy
    */

    // void AdvanceParticle(double dr, double ei, double ef);

    //----------------------------------------------------------------------------//

    void swap(Propagator& propagator);

    // void ReadConfigFile(std::string config_file, bool DoApplyOptions=true);


    // --------------------------------------------------------------------- //
    // Getter
    // --------------------------------------------------------------------- //

    // --------------------------------------------------------------------- //
    const Sector* GetCurrentCollection() const { return current_sector_; }

    //----------------------------------------------------------------------------//
    std::vector<Sector*> GetSectors() const { return sectors_; }

    int GetSeed() const;
    void SetSeed(int seed);
    Geometry* GetDetector() const;
    Particle& GetParticle();
    // void SetDetector(Geometry *detector);
    // bool GetStopping_decay() const;
    // void SetStopping_decay(bool stopping_decay);

    private:

    // ----------------------------------------------------------------------------
    /// @brief Simple wrapper to initialize propagator from config file
    ///
    /// The value of var will be treated as a default value.
    ///
    /// @param var: the variable to initialize
    /// @param option: option in the property_tree
    /// @param property_tree
    // ----------------------------------------------------------------------------
    template<class T>
    void SetMember(T& var, std::string option, boost::property_tree::ptree& pt)
    {
        boost::optional<T> optional_param = pt.get_optional<T>(option);
        if (optional_param)
        {
            var = optional_param.get();
        }
        else
        {
            std::stringstream ss;
            ss << "Option " << option << " not set! Use default: " << var;
            log_warn("%s", ss.str().c_str());
        }
        // try
        // {
        //     var = pt.get<T>(option);
        // }
        // catch(std::exception& ex)
        // {
        //     std::stringstream ss;
        //     ss<<ex.what()<<"! Use default: "<<var;
        //     log_warn("%s", ss.str().c_str());
        // }
    }

    // ----------------------------------------------------------------------------
    /// @brief Choose the current collection the particle is in.
    ///
    /// @param particle_position
    /// @param particle_direction
    // ----------------------------------------------------------------------------
    void ChooseCurrentCollection(const Vector3D& particle_position, const Vector3D& particle_direction);

    // ----------------------------------------------------------------------------
    /// @brief Calculate the distance to propagate
    ///
    /// Calculate the distance to propagate and
    /// choose if the particle has to propagate through the whole sector
    /// or only to the collection border
    ///
    /// @param particle_position
    /// @param particle_direction
    ///
    /// @return distance
    // ----------------------------------------------------------------------------
    double CalculateEffectiveDistance(const Vector3D& particle_position, const Vector3D& particle_direction);

    int seed_; //!< seed of the random number generator
    // ParametrizationType::Enum  brems_;                     //!< Bremsstrahlungs parametrization
    // ParametrizationType::Enum  photo_;                     //!< Photonuclear parametrization
    static const double global_ecut_inside_; //!< ecut for inside the detector (it's used when not specified explicit for a sector in
    //!congiguration file)
    static const double global_ecut_infront_; //!< ecut for infront of the detector (it's used when not specified explicit for a
    //!sector in congiguration file)
    static const double global_ecut_behind_; //!< ecut for behind the detector (it's used when not specified explicit for a sector in
    //!congiguration file)
    static const double global_vcut_inside_; //!< vcut for inside the detector (it's used when not specified explicit for a sector in
    //!congiguration file)
    static const double global_vcut_infront_; //!< ecut for infront of the detector (it's used when not specified explicit for a
    //!sector in congiguration file)
    static const double global_vcut_behind_; //!< ecut for behind the detector (it's used when not specified explicit for a sector in
    //!congiguration file)
    static const double global_cont_inside_; //!< continuous randominzation flag for inside the detector (it's used when not
    //!specified explicit for a sector in congiguration file)
    static const double global_cont_infront_; //!< continuous randominzation flag for infront of the detector (it's used when not
    //!specified explicit for a sector in congiguration file)
    static const double global_cont_behind_;  //!< continuous randominzation flag for behind the detector (it's used when not
                                 //!specified explicit for a sector in congiguration file)
    static const bool interpolate_;  //!< Enable interpolation

    std::vector<Sector*> sectors_;
    Sector* current_sector_;

    Particle particle_;
    Geometry* detector_;
};

}
