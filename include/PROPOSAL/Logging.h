
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

#include <spdlog/spdlog.h>
#include <spdlog/sinks/sink.h>
#include <PROPOSAL/methods.h>

#include <unordered_map>
#include <memory>
#include <string>

using std::unordered_map;
using std::unique_ptr;
using std::shared_ptr;
using std::string;
using std::make_shared;

namespace PROPOSAL {
struct Logging {
    static unordered_map<std::string, unique_ptr<spdlog::logger>> logger;
    static shared_ptr<spdlog::sinks::sink> sink;

    Logging() = delete;

    static unique_ptr<spdlog::logger> Create(std::string const& name)
    {
        return PROPOSAL::make_unique<spdlog::logger>(name, sink);
    }

    static spdlog::logger* Get(std::string const& name)
    {
        auto it = Logging::logger.find(name);
        if (it == logger.end())
            Logging::logger[name] = Logging::Create(name);
        return Logging::logger[name].get();
    }
};
} // namespace PROPOSAL
