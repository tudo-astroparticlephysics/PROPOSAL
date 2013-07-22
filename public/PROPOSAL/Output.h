#ifndef OUTPUT_H
#define OUTPUT_H
#ifndef ICECUBE

//Stuff for LOG4CPLUS
#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
#include <log4cplus/configurator.h>
#include <iomanip>
#include <string>
#include <sstream>
#include <boost/preprocessor/control/if.hpp>
#include <boost/preprocessor/comparison/less.hpp>
#include <boost/preprocessor/array/elem.hpp>
#include <boost/preprocessor/comparison/equal.hpp>
#include <boost/preprocessor/config/limits.hpp>
#include <iostream>
#include <vector>
#include <string>
#include "PROPOSAL/Particle.h"
#include "utility"

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

class Output
{
private:
    Output()
    {
        PropertyConfigurator::doConfigure(LOG4CPLUS_TEXT("resources/log4cplus.conf"));
        logger = Logger::getInstance(LOG4CPLUS_TEXT("PORPOSAL"));
    }

    Output(Output const&);         // Don't Implement.
    void operator=(Output const&); // Don't implement

    static std::vector<Particle*> secondarys_;

public:

    Logger logger;

    static Output& getInstance()
    {
        static Output instance;
        return instance;
    }

    void SetLoggingConfigurationFile(std::string file);

    void FillSecondaryVector(Particle *particle, int secondary_id, std::pair<double,std::string> energy_loss, double distance);


    //Getter
    std::vector<Particle*> GetSecondarys() const
    {
        return secondarys_;
    }

 };

#endif //ICECUBE

#endif //OUTPUT_H























