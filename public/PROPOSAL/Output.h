
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

#ifndef OUTPUT_H
#define OUTPUT_H
#ifndef ICECUBE


// #include <iostream>
#include <vector>
// #include <string>
// #include <utility>
// #include <iomanip>
// #include <sstream>
#include <fstream>

#include "PROPOSAL/PROPOSALParticle.h"

#if ROOT_SUPPORT
    #include "TTree.h"
    #include "TFile.h"
#endif

#if PROPOSAL_STANDALONE
    #if LOG4CPLUS_SUPPORT
        //Stuff for LOG4CPLUS
        #include <log4cplus/logger.h>
        #include <log4cplus/loggingmacros.h>
        #include <log4cplus/configurator.h>
        #include <boost/preprocessor/control/if.hpp>
        #include <boost/preprocessor/comparison/less.hpp>
        #include <boost/preprocessor/array/elem.hpp>
        #include <boost/preprocessor/comparison/equal.hpp>
        #include <boost/preprocessor/config/limits.hpp>

        #define VA_COUNT(...) VA_COUNT0(__VA_ARGS__, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1)
        #define VA_COUNT0(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, N,...) N
        #define VA_ARG(_i, ...) BOOST_PP_ARRAY_ELEM(_i, (VA_COUNT(__VA_ARGS__), (__VA_ARGS__)))

        #define log_error(...) BOOST_PP_IF( BOOST_PP_EQUAL( 1, VA_COUNT(__VA_ARGS__) ),\
            LOG4CPLUS_ERROR(Output::getInstance().logger, VA_ARG(0, __VA_ARGS__) ) , \
            LOG4CPLUS_ERROR_FMT(Output::getInstance().logger, __VA_ARGS__ ) )

        #define log_fatal(...) BOOST_PP_IF( BOOST_PP_EQUAL( 1, VA_COUNT(__VA_ARGS__) ),\
            LOG4CPLUS_FATAL(Output::getInstance().logger, VA_ARG(0, __VA_ARGS__) ) ; exit(1) , \
            LOG4CPLUS_FATAL_FMT(Output::getInstance().logger, __VA_ARGS__ ) ) ; exit(1)

        #define log_warn(...) BOOST_PP_IF( BOOST_PP_EQUAL( 1, VA_COUNT(__VA_ARGS__) ),\
            LOG4CPLUS_WARN(Output::getInstance().logger, VA_ARG(0, __VA_ARGS__) ) , \
            LOG4CPLUS_WARN_FMT(Output::getInstance().logger, __VA_ARGS__ ) )

        #define log_info(...) BOOST_PP_IF( BOOST_PP_EQUAL( 1, VA_COUNT(__VA_ARGS__) ),\
            LOG4CPLUS_INFO(Output::getInstance().logger, VA_ARG(0, __VA_ARGS__) ) , \
            LOG4CPLUS_INFO_FMT(Output::getInstance().logger, __VA_ARGS__ ) )

        #define log_trace(...) BOOST_PP_IF( BOOST_PP_EQUAL( 1, VA_COUNT(__VA_ARGS__) ),\
            LOG4CPLUS_TRACE(Output::getInstance().logger, VA_ARG(0, __VA_ARGS__) ) , \
            LOG4CPLUS_TRACE_FMT(Output::getInstance().logger, __VA_ARGS__ ) )

        #define log_debug(...) BOOST_PP_IF( BOOST_PP_EQUAL( 1, VA_COUNT(__VA_ARGS__) ),\
            LOG4CPLUS_DEBUG(Output::getInstance().logger, VA_ARG(0, __VA_ARGS__) ) , \
            LOG4CPLUS_DEBUG_FMT(Output::getInstance().logger, __VA_ARGS__ ) )

        #define log_notice(...) BOOST_PP_IF( BOOST_PP_EQUAL( 1, VA_COUNT(__VA_ARGS__) ),\
            LOG4CPLUS_NOTICE(Output::getInstance().logger, VA_ARG(0, __VA_ARGS__) ) , \
            LOG4CPLUS_NOTICE_FMT(Output::getInstance().logger, __VA_ARGS__ ) )

        using namespace log4cplus;
    #else
        #ifndef log_error
            #define log_error(fmt, ...) (void)0
            #define log_fatal(fmt, ...) (void)0
            #define log_warn(fmt, ...) (void)0
            #define log_info(fmt, ...) (void)0
            #define log_trace(fmt, ...) (void)0
            #define log_debug(fmt, ...) (void)0
            #define log_notice(fmt, ...) (void)0
        #endif
    #endif
#else
#include <icetray/I3Logging.h>
#endif




namespace PROPOSAL{

class Output
{
private:

    Output()
    {
        #if LOG4CPLUS_SUPPORT
        // PropertyConfigurator::doConfigure(LOG4CPLUS_TEXT("resources/log4cplus.conf"));
        PropertyConfigurator::doConfigure(LOG4CPLUS_TEXT(LOG4CPLUS_CONFIG));
        logger = Logger::getInstance(LOG4CPLUS_TEXT("PROPOSAL"));
        #endif
    }

    Output(Output const&);         // Don't Implement.
    void operator=(Output const&); // Don't implement

    static std::vector<PROPOSALParticle*> secondarys_;

    static bool store_in_root_trees_;

    //ASCII Output
    std::ofstream secondary_ascii_;
    std::ofstream primary_ascii_;
    std::ofstream propagated_primary_ascii_;

    #if ROOT_SUPPORT

        TTree *secondary_tree_;
        TTree *primary_tree_;
        TTree *propagated_primary_tree_;
        TFile *rootfile_;

        double secondary_x_;
        double secondary_y_;
        double secondary_z_;
        double secondary_t_;
        double secondary_theta_;
        double secondary_phi_;
        double secondary_energy_;
        int secondary_parent_particle_id_;
        int secondary_particle_id_;
        std::string secondary_name_;
        double current_primary_energy_;

        double primary_x_;
        double primary_y_;
        double primary_z_;
        double primary_t_;
        double primary_theta_;
        double primary_phi_;
        double primary_energy_;
        int primary_parent_particle_id_;
        int primary_particle_id_;
        std::string primary_name_;

        double prop_primary_x_;
        double prop_primary_y_;
        double prop_primary_z_;
        double prop_primary_t_;
        double prop_primary_theta_;
        double prop_primary_phi_;
        double prop_primary_energy_;
        int prop_primary_parent_particle_id_;
        int prop_primary_particle_id_;
        std::string prop_primary_name_;
        double prop_primary_propagated_distance_;
        double prop_primary_xi_;
        double prop_primary_yi_;
        double prop_primary_zi_;
        double prop_primary_ti_;
        double prop_primary_ei_;
        double prop_primary_xf_;
        double prop_primary_yf_;
        double prop_primary_zf_;
        double prop_primary_tf_;
        double prop_primary_ef_;
        double prop_primary_xc_;
        double prop_primary_yc_;
        double prop_primary_zc_;
        double prop_primary_tc_;
        double prop_primary_ec_;
        double prop_primary_elost_;

    #endif

public:

    #if LOG4CPLUS_SUPPORT
        Logger logger;
    #endif
    //ASCII
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

    void FillSecondaryVector(PROPOSALParticle *particle, int secondary_id, std::pair<double, ParticleType::Enum> energy_loss, double distance);

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
        void StorePrimaryInTree(PROPOSALParticle *primary);

//----------------------------------------------------------------------------//

        void StorePropagatedPrimaryInTree(PROPOSALParticle *prop_primary);

    #endif

// ASCII OUTPUT
        void EnableASCIIOutput(std::string ASCII_Prefix, bool append = false);

     //----------------------------------------------------------------------------//

         void DisableASCIIOutput();

     //----------------------------------------------------------------------------//

         void StorePrimaryInASCII(PROPOSALParticle *primary);

     //----------------------------------------------------------------------------//

         void StorePropagatedPrimaryInASCII(PROPOSALParticle *prop_primary);

     //----------------------------------------------------------------------------//

         void WriteDescriptionFile();


    //Getter
    std::vector<PROPOSALParticle*> GetSecondarys() const
    {
        return secondarys_;
    }

 };
}

#endif //ICECUBE

#endif //OUTPUT_H

