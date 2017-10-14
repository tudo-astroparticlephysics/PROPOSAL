#pragma once

#ifndef ICECUBE


#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "PROPOSAL/particle/Particle.h"

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

        // using namespace log4cplus;  // no using directives in header
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
        // log4cplus::PropertyConfigurator::doConfigure(LOG4CPLUS_TEXT("resources/log4cplus.conf"));
        log4cplus::PropertyConfigurator::doConfigure(LOG4CPLUS_TEXT(LOG4CPLUS_CONFIG));
        logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("PROPOSAL"));
        #endif
    }

    Output(Output const&);         // Don't Implement.
    void operator=(Output const&); // Don't implement

    static std::vector<DynamicData*> secondarys_;

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

    void FillSecondaryVector(Particle *particle, int secondary_id, std::pair<double, ParticleType::Enum> energy_loss, double distance);

    void FillSecondaryVector(const Particle& particle, const DynamicData::Type& secondary, double energyloss, double distance);
    void FillSecondaryVector(std::vector<Particle*>&);

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
        void StorePrimaryInTree(Particle *primary);

//----------------------------------------------------------------------------//

        void StorePropagatedPrimaryInTree(Particle *prop_primary);

    #endif

// ASCII OUTPUT
        void EnableASCIIOutput(std::string ASCII_Prefix, bool append = false);

     //----------------------------------------------------------------------------//

         void DisableASCIIOutput();

     //----------------------------------------------------------------------------//

         void StorePrimaryInASCII(Particle *primary);

     //----------------------------------------------------------------------------//

         void StorePropagatedPrimaryInASCII(Particle *prop_primary);

     //----------------------------------------------------------------------------//

         void WriteDescriptionFile();


    //Getter
    std::vector<DynamicData*> GetSecondarys() const
    {
        return secondarys_;
    }

 };
}

#endif //ICECUBE
