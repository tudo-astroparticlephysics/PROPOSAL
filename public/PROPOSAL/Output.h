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
#include "boost/preprocessor/facilities/empty.hpp"

#define log_error(para,...) LOG4CPLUS_ERROR_FMT(Output::getInstance().logger, para ,##__VA_ARGS__)
#define log_fatal(para,...) LOG4CPLUS_FATAL_FMT(Output::getInstance().logger, para ,##__VA_ARGS__)
#define log_arsch(para,...) if(BOOST_PP_EMPTY) cout<<"ss.str()"<<endl;

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
