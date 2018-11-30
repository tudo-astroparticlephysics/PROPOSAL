
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

#ifndef ICECUBE

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "PROPOSAL/particle/Particle.h"

#if ROOT_SUPPORT
    #include "TFile.h"
    #include "TTree.h"
#endif

#if not PROPOSAL_STANDALONE
    #include <icetray/I3Logging.h>
#else // not PROPOSAL_STANDALONE
    #if LOG4CPLUS_SUPPORT
        #include <log4cplus/configurator.h>
        #include <log4cplus/logger.h>
        #include <log4cplus/loggingmacros.h>
    #endif // log4cplus
#endif // not PROPOSAL_STANDALONE

namespace PROPOSAL {

class Output
{
private:
    Output()
    {
#if LOG4CPLUS_SUPPORT
        log4cplus::PropertyConfigurator::doConfigure(LOG4CPLUS_TEXT(LOG4CPLUS_CONFIG));
        logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("PROPOSAL"));
#endif
    }

    Output(Output const&);         // Don't Implement.
    void operator=(Output const&); // Don't implement

    static std::vector<DynamicData*> secondarys_;

    static bool store_in_root_trees_;

    // ASCII Output
    std::ofstream secondary_ascii_;
    std::ofstream primary_ascii_;
    std::ofstream propagated_primary_ascii_;

#if ROOT_SUPPORT

    TTree* secondary_tree_;
    TTree* primary_tree_;
    TTree* propagated_primary_tree_;
    TFile* rootfile_;

    Vector3D secondary_position_;
    double secondary_t_;
    Vector3D secondary_direction_;
    double secondary_energy_;
    int secondary_parent_particle_id_;
    int secondary_particle_id_;
    std::string secondary_name_;
    double current_primary_energy_;

    Vector3D primary_position_;
    double primary_t_;
    Vector3D primary_direction_;
    double primary_energy_;
    int primary_parent_particle_id_;
    int primary_particle_id_;
    std::string primary_name_;

    Vector3D prop_primary_position_;
    double prop_primary_t_;
    Vector3D prop_primary_direction_;
    double prop_primary_energy_;
    int prop_primary_parent_particle_id_;
    int prop_primary_particle_id_;
    std::string prop_primary_name_;
    double prop_primary_propagated_distance_;
    Vector3D prop_primary_entry_point_;
    double prop_primary_ti_;
    double prop_primary_ei_;
    Vector3D prop_primary_exit_point_;
    double prop_primary_tf_;
    double prop_primary_ef_;
    Vector3D prop_primary_closest_approach_point_;
    double prop_primary_tc_;
    double prop_primary_ec_;
    double prop_primary_elost_;

#endif

public:
#if LOG4CPLUS_SUPPORT
    log4cplus::Logger logger;
#endif
    // ASCII
    static bool store_in_ASCII_file_;

    //----------------------------------------------------------------------------//

    static Output& getInstance()
    {
        static Output instance;
        return instance;
    }
    //----------------------------------------------------------------------------//

    void SetLoggingConfigurationFile(std::string file);

    //----------------------------------------------------------------------------//

    // ----------------------------------------------------------------------------
    /// @brief Fill secondary data
    ///
    /// This method is used to store decay products.
    ///
    /// @param std::vector
    // ----------------------------------------------------------------------------
    void FillSecondaryVector(const std::vector<Particle*>&);

    // ----------------------------------------------------------------------------
    /// @brief Fill secondary data
    ///
    /// New DynamicData is create with information of the given particle, type of
    /// DynamicData and the energyloss
    ///
    /// @param particle
    /// @param secondary
    /// @param energyloss
    // ----------------------------------------------------------------------------
    void FillSecondaryVector(const Particle& particle, const DynamicData::Type& secondary, double energyloss);

    // ----------------------------------------------------------------------------
    /// @brief Fill secondary data
    ///
    /// New DynamicData with Type ContinuousEnergyLoss is createed.
    /// The position is, where the continuous loss starts.
    /// The direction is the unnormed difference between the positions of the 2 stochastic losses.
    /// The time is the differnce between the 2 stochastic losses.
    /// The PropagatedDistance is the distance between the 2 stochastic losses.
    /// The Energy is the energy lost during the continuous loss (including the continuous randomization).
    /// The ParentParticleEnergy is the particle energy before the continuous loss.
    ///
    /// @param ContinuousEnergyLoss
    // ----------------------------------------------------------------------------
    void FillSecondaryVector(DynamicData* continuous_loss);

    //----------------------------------------------------------------------------//

    void ClearSecondaryVector();

    //----------------------------------------------------------------------------//

    void Close();

    //----------------------------------------------------------------------------//

    void EnableROOTOutput(std::string rootfile_name);

    //----------------------------------------------------------------------------//

    void DisableROOTOutput();

    //----------------------------------------------------------------------------//
#if ROOT_SUPPORT
    void StorePrimaryInTree(Particle* primary);

    //----------------------------------------------------------------------------//

    void StorePropagatedPrimaryInTree(Particle* prop_primary);

#endif

    // ASCII OUTPUT
    void EnableASCIIOutput(std::string ASCII_Prefix, bool append = false);

    //----------------------------------------------------------------------------//

    void DisableASCIIOutput();

    //----------------------------------------------------------------------------//

    void StorePrimaryInASCII(Particle* primary);

    //----------------------------------------------------------------------------//

    void StorePropagatedPrimaryInASCII(Particle* prop_primary);

    //----------------------------------------------------------------------------//

    void WriteDescriptionFile();

    // Getter
    std::vector<DynamicData*> GetSecondarys() const { return secondarys_; }
};
} // namespace PROPOSAL


#if PROPOSAL_STANDALONE

#if LOG4CPLUS_SUPPORT

template<typename... Args>
void log_error(Args ... args)
{
    LOG4CPLUS_ERROR_FMT(PROPOSAL::Output::getInstance().logger, args...);
}

template<typename... Args>
void log_fatal(Args ... args)
{
    LOG4CPLUS_FATAL_FMT(PROPOSAL::Output::getInstance().logger, args...);
    exit(1);
}

template<typename... Args>
void log_warn(Args ... args)
{
    LOG4CPLUS_WARN_FMT(PROPOSAL::Output::getInstance().logger, args...);
}

template<typename... Args>
void log_info(Args ... args)
{
    LOG4CPLUS_INFO_FMT(PROPOSAL::Output::getInstance().logger, args...);
}

template<typename... Args>
void log_trace(Args ... args)
{
    LOG4CPLUS_TRACE_FMT(PROPOSAL::Output::getInstance().logger, args...);
}

template<typename... Args>
void log_debug(Args ... args)
{
    LOG4CPLUS_DEBUG_FMT(PROPOSAL::Output::getInstance().logger, args...);
}

template<typename... Args>
void log_notice(Args ... args)
{
    LOG4CPLUS_NOTICE_FMT(PROPOSAL::Output::getInstance().logger, args...);
}

#else // log4cplus

template<typename... Args>
void log_error(Args ... args)
{
    void(args...);
}

template<typename... Args>
void log_fatal(Args ... args)
{
    printf("FATAL ERROR: ");
    printf(args...);
    exit(1);
}

template<typename... Args>
void log_warn(Args ... args)
{
    void(args...);
}

template<typename... Args>
void log_info(Args ... args)
{
    void(args...);
}

template<typename... Args>
void log_trace(Args ... args)
{
    void(args...);
}

template<typename... Args>
void log_debug(Args ... args)
{
    void(args...);
}

template<typename... Args>
void log_notice(Args ... args)
{
    void(args...);
}

#endif // log4cplus
#endif // PROPOSAL_STANDALONE

#endif // ICECUBE
