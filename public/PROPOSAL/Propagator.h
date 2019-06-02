
/******************************************************************************
 *                                                                            *
 * This file is part of the simulation tool PROPOSAL.                         *
 *                                                                            *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,          *
 *                    Chair Experimental Physics 5b                           *
 *                                                                            *
 * This software may be modified and distributed under the terms of a         *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE".                                     *
 *                                                                            *
 * Modifcations to the LGPL License:                                          *
 *                                                                            *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the       *
 *         following reference:                                               *
 *                                                                            *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001                                          *
 *                                                                            *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *         GitHub webpage                                                     *
 *                                                                            *
 *         "https://github.com/tudo-astroparticlephysics/PROPOSAL"            *
 *                                                                            *
 ******************************************************************************/


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

// #include <deque>
#include <vector>

#include "PROPOSAL/Output.h"
#include "PROPOSAL/sector/Sector.h"

namespace PROPOSAL {

class Propagator
{
public:
    // Constructors
    Propagator(const std::vector<Sector*>&, const Geometry&);
    Propagator(const ParticleDef&, const std::vector<Sector::Definition>&, const Geometry&);
    Propagator(const ParticleDef&, const std::vector<Sector::Definition>&, const Geometry&, const InterpolationDef&);
    Propagator(const ParticleDef&, const std::string&);

    Propagator(const Propagator&);
    ~Propagator();

    bool operator==(const Propagator& propagator) const;
    bool operator!=(const Propagator& propagator) const;

    // ----------------------------------------------------------------------------
    /// @brief Propagates the particle through the current set of Sectors
    ///
    /// @param MaxDistance_cm
    ///
    /// @return Secondary data
    // ----------------------------------------------------------------------------
    std::vector<DynamicData*> Propagate(double MaxDistance_cm = 1e20);

    // --------------------------------------------------------------------- //
    // Getter
    // --------------------------------------------------------------------- //

    const Sector* GetCurrentSector() const { return current_sector_; }
    const std::vector<Sector*> GetSectors() const { return sectors_; }

    Geometry& GetDetector() const { return *detector_; };
    Particle& GetParticle() { return particle_; };

private:

    Propagator& operator=(const Propagator& propagator);

    // ----------------------------------------------------------------------------
    /// @brief Simple wrapper to initialize propagator from config file
    ///
    /// The value of var will be treated as a default value.
    ///
    /// @param var: the variable to initialize
    /// @param json_key: key in the json_object
    /// @param json_object
    // ----------------------------------------------------------------------------
    // template<class T>
    // void SetMember(T& var, const std::string& json_key, const nlohmann::json& json_object)
    // {
    //     if (json_object.find(json_key) != json_object.end())
    //     {
    //         var = json_object[json_key];
    //     } else
    //     {
    //         std::stringstream ss;
    //         ss << "Option " << json_key << " not set! Use default: " << var;
    //         log_debug("%s", ss.str().c_str());
    //     }
    // }

    InterpolationDef CreateInterpolationDef(const std::string& json_object_str);
    Sector::Definition CreateSectorDefinition(const std::string& json_object_str);
    std::string ParseCutSettings(const std::string& json_object_str,
                                 const std::string& json_key,
                                 double default_ecut,
                                 double default_vcut,
                                 bool default_contrand);

    // ----------------------------------------------------------------------------
    /// @brief Create geometry from json config file
    ///
    /// @param pt boost property tree
    ///
    /// @return new Geometry
    // ----------------------------------------------------------------------------
    Geometry* ParseGeometryConifg(const std::string& json_object_str);

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

    // --------------------------------------------------------------------- //
    // Global default values
    // --------------------------------------------------------------------- //

    static const int
        global_seed_; //!< Seed for the internal random number generator
    static const double
        global_ecut_inside_; //!< ecut for inside the detector (it's used when not specified explicit for a sector in
    //! congiguration file)
    static const double
        global_ecut_infront_; //!< ecut for infront of the detector (it's used when not specified explicit for a
    //! sector in congiguration file)
    static const double
        global_ecut_behind_; //!< ecut for behind the detector (it's used when not specified explicit for a sector in
    //! congiguration file)
    static const double
        global_vcut_inside_; //!< vcut for inside the detector (it's used when not specified explicit for a sector in
    //! congiguration file)
    static const double
        global_vcut_infront_; //!< ecut for infront of the detector (it's used when not specified explicit for a
    //! sector in congiguration file)
    static const double
        global_vcut_behind_; //!< ecut for behind the detector (it's used when not specified explicit for a sector in
    //! congiguration file)
    static const double
        global_cont_inside_; //!< continuous randominzation flag for inside the detector (it's used when not
    //! specified explicit for a sector in congiguration file)
    static const double
        global_cont_infront_; //!< continuous randominzation flag for infront of the detector (it's used when not
    //! specified explicit for a sector in congiguration file)
    static const double
        global_cont_behind_;        //!< continuous randominzation flag for behind the detector (it's used when not
                                    //! specified explicit for a sector in congiguration file)
    static const bool do_interpolation_; //!< Enable interpolation
    static const bool uniform_; //!< Enable uniform sampling of phase space points for decays

    // --------------------------------------------------------------------- //
    // Private Member
    // --------------------------------------------------------------------- //

    std::vector<Sector*> sectors_;
    Sector* current_sector_;

    Particle particle_;
    Geometry* detector_;
};

} // namespace PROPOSAL

