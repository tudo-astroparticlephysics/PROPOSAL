
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

#include <iostream>
#include <cstdio>
#include <string>

#ifdef ICECUBE_PROJECT
    #include <icetray/I3Logging.h>
#endif // ICECUBE_PROJECT

#if LOG4CPLUS_SUPPORT
    #include <log4cplus/configurator.h>
    #include <log4cplus/logger.h>
    #include <log4cplus/loggingmacros.h>
#endif // log4cplus

namespace PROPOSAL {

class Logging
{
private:
    Logging()
    {
#if LOG4CPLUS_SUPPORT
        log4cplus::PropertyConfigurator::doConfigure(LOG4CPLUS_TEXT(LOG4CPLUS_CONFIG));
        logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("PROPOSAL"));
#endif
    }

    Logging(Logging const&);         // Don't Implement.
    void operator=(Logging const&); // Don't implement

public:
#if LOG4CPLUS_SUPPORT
    log4cplus::Logger logger;
#endif

    static Logging& getInstance()
    {
        static Logging instance;
        return instance;
    }

    void SetLoggingConfigurationFile(std::string file)
    {
#if LOG4CPLUS_SUPPORT
        log4cplus::PropertyConfigurator::doConfigure(LOG4CPLUS_TEXT(file));
#else
        (void) file;
        std::cout << "Log4cplus not found! No log messages will be shown!" << std::endl;
#endif
    }

};
} // namespace PROPOSAL


#if LOG4CPLUS_SUPPORT

template<typename... Args>
void log_error(Args ... args)
{
    LOG4CPLUS_ERROR_FMT(PROPOSAL::Logging::getInstance().logger, args...);
}

template<typename... Args>
void log_fatal(Args ... args)
{
    LOG4CPLUS_FATAL_FMT(PROPOSAL::Logging::getInstance().logger, args...);
    exit(1);
}

template<typename... Args>
void log_warn(Args ... args)
{
    LOG4CPLUS_WARN_FMT(PROPOSAL::Logging::getInstance().logger, args...);
}

template<typename... Args>
void log_info(Args ... args)
{
    LOG4CPLUS_INFO_FMT(PROPOSAL::Logging::getInstance().logger, args...);
}

template<typename... Args>
void log_trace(Args ... args)
{
    LOG4CPLUS_TRACE_FMT(PROPOSAL::Logging::getInstance().logger, args...);
}

template<typename... Args>
void log_debug(Args ... args)
{
    LOG4CPLUS_DEBUG_FMT(PROPOSAL::Logging::getInstance().logger, args...);
}

template<typename... Args>
void log_notice(Args ... args)
{
    LOG4CPLUS_NOTICE_FMT(PROPOSAL::Logging::getInstance().logger, args...);
}

#else // log4cplus
#ifndef ICECUBE_PROJECT

template<typename... Args>
void log_error(Args ... args)
{
    printf(args...);
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
    printf(args...);
}

template<typename... Args>
void log_info(Args ... args)
{
    printf(args...);
}

template<typename... Args>
void log_trace(Args ... args)
{
    printf(args...);
}

template<typename... Args>
void log_debug(Args ... args)
{
    printf(args...);
}

template<typename... Args>
void log_notice(Args ... args)
{
    printf(args...);
}

#endif // not ICECUBE_PROJECT
#endif // log4cplus
