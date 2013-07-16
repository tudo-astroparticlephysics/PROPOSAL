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
#include <boost/preprocessor/comparison/equal.hpp>

#define VA_COUNT(...) VA_COUNT0(__VA_ARGS__, 9, 8, 7, 6, 5, 4, 3, 2, 1)
#define VA_COUNT0(_1, _2, _3, _4, _5, _6, _7, _8, _9, N,...) N

#define log_error(para,...) LOG4CPLUS_ERROR_FMT(Output::getInstance().logger, para ,##__VA_ARGS__)
#define log_fatal(para,...) LOG4CPLUS_FATAL_FMT(Output::getInstance().logger, para ,##__VA_ARGS__)
//#define log_arsch(...) {if(VA_COUNT == 1) cout<<"ss.str()"<<endl; else cout<<"HALLO"<<endl}
#define log_arsch(...) BOOST_PP_IF( BOOST_PP_EQUAL(1, VA_COUNT(__VA_ARGS__)), \
    LOG4CPLUS_WARN(Output::getInstance().logger, LOG4CPLUS_TEXT(__VA_ARGS__) )   ,\
    LOG4CPLUS_WARN_FMT(Output::getInstance().logger, __VA_ARGS__ ))

#define log_info(para,...) LOG4CPLUS_INFO_FMT(Output::getInstance().logger, para ,##__VA_ARGS__)
#define log_trace(para,...) LOG4CPLUS_TRACE_FMT(Output::getInstance().logger, para ,##__VA_ARGS__)
#define log_debug(para,...) LOG4CPLUS_DEBUG_FMT(Output::getInstance().logger, para , ##__VA_ARGS__)
#define log_notice(para,...) LOG4CPLUS_NOTICE_FMT(Output::getInstance().logger, para , ##__VA_ARGS__)

//#define log_error(para,...) LOG4CPLUS_ERROR_FMT(Output::getInstance().logger, LOG4CPLUS_TEXT(para),##__VA_ARGS__)
//#define log_fatal(para,...) LOG4CPLUS_FATAL_FMT(Output::getInstance().logger, LOG4CPLUS_TEXT(para),##__VA_ARGS__)
#define log_warn(para,...) LOG4CPLUS_WARN_FMT(Output::getInstance().logger, LOG4CPLUS_TEXT(para),##__VA_ARGS__)
//#define log_info(para,...) LOG4CPLUS_INFO_FMT(Output::getInstance().logger, LOG4CPLUS_TEXT(para),##__VA_ARGS__)
//#define log_trace(para,...) LOG4CPLUS_TRACE_FMT(Output::getInstance().logger, LOG4CPLUS_TEXT(para),##__VA_ARGS__)
//#define log_debug(para,...) LOG4CPLUS_DEBUG_FMT(Output::getInstance().logger, LOG4CPLUS_TEXT(para),##__VA_ARGS__)
//#define log_notice(para,...) LOG4CPLUS_NOTICE_FMT(Output::getInstance().logger, LOG4CPLUS_TEXT(para),##__VA_ARGS__)


//    LOG4CPLUS_WARN(Output::getInstance().logger, __VA_ARGS__ )

using namespace log4cplus;

class Output
{
    public:
        static Output& getInstance()
        {
            static Output instance;
            return instance;
        }
        Logger logger;

        void SetLoggingConfigurationFile(std::string file);

    private:
        Output() {
            PropertyConfigurator::doConfigure(LOG4CPLUS_TEXT("resources/log4cplus.conf"));
            logger = Logger::getInstance(LOG4CPLUS_TEXT("PORPOSAL"));
        }
        Output(Output const&);              // Don't Implement.
        void operator=(Output const&); // Don't implement
 };

#endif //ICECUBE

#endif //OUTPUT_H
