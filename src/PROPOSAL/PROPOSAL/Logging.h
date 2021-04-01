
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

#include <PROPOSAL/methods.h>
#include <spdlog/sinks/sink.h>
#include <spdlog/spdlog.h>

#include <memory>
#include <string>
#include <unordered_map>

namespace PROPOSAL {

using logger_ptr = std::shared_ptr<spdlog::logger>;

struct Logging {
    static std::unordered_map<std::string, logger_ptr> logger;
    static std::shared_ptr<spdlog::sinks::sink> sink;

    Logging() = delete;

    static logger_ptr Get(std::string const& name)
    {
        auto it = Logging::logger.find(name);
        if (it == logger.end())
            Logging::logger[name] = Logging::Create(name);
        return Logging::logger[name];
    }

    static void SetGlobalLoglevel(spdlog::level::level_enum loglevel)
    {
        for (auto& l : logger)
            l.second->set_level(loglevel);
        global_loglevel = loglevel;
    }

private:
    static logger_ptr Create(std::string const& name)
    {
        auto logger = std::make_shared<spdlog::logger>(name, sink);
        logger->set_level(global_loglevel);
        return logger;
    }

    static spdlog::level::level_enum global_loglevel;
};
} // namespace PROPOSAL
