#ifndef OUTPUT_H
#define OUTPUT_H
#ifndef ICECUBE

//Stuff for LOG4CPLUS
#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
#include <log4cplus/configurator.h>
#include <iomanip>
#include <string>

#define log_error(para) LOG4CPLUS_ERROR(Output::getInstance().logger, LOG4CPLUS_TEXT(para))
#define log_fatal(para) LOG4CPLUS_FATAL(Output::getInstance().logger, LOG4CPLUS_TEXT(para))
#define log_warn(para) LOG4CPLUS_WARN(Output::getInstance().logger, LOG4CPLUS_TEXT(para))
#define log_info(para) LOG4CPLUS_INFO(Output::getInstance().logger, LOG4CPLUS_TEXT(para))
#define log_trace(para) LOG4CPLUS_TRACE(Output::getInstance().logger, LOG4CPLUS_TEXT(para))
#define log_debug(para) LOG4CPLUS_DEBUG(Output::getInstance().logger, LOG4CPLUS_TEXT(para))
#define log_notice(para) LOG4CPLUS_NOTICE(Output::getInstance().logger, LOG4CPLUS_TEXT(para))




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
