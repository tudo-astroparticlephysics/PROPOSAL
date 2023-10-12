/*! \file   methods.cxx
 *   \brief  Source file for the methods routines.
 *
 *   For more details see the class documentation.
 *
 *   \date   21.06.2010
 *   \author Jan-Hendrik Koehne
 */

// #include <stdlib.h>

#include "PROPOSAL/methods.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/Constants.h"
#include <algorithm>
#include <string>

namespace PROPOSAL {

std::string LogTableCreation::warn_for_path = "";

LogTableCreation::LogTableCreation(const std::string &path, const std::string &filename) {
    // TODO: use std::filesystem when we switch to c++17
    auto combined = path + "/" + filename;
    if (!Helper::file_exists(combined)) {
        if (warn_for_path != path) {
            warn_for_path = path;
            // we haven't logged a warning for this specific path yet
            if (InterpolationSettings::ALLOW_TABLE_CREATION) {
                if (!Helper::is_folder_writable(path)) {
                    Logging::Get("TableCreation")->warn("PROPOSAL is unable to write to the requested path '{}'", path);
                } else {
                    Logging::Get("TableCreation")->warn("Tables are not available and need to be created. "
                                                        "They will be written to '{}'. "
                                                        "This can take some minutes.", path);
                    return;
                }
            } else {
                Logging::Get("TableCreation")->warn("Creation of new tables not allowed due to "
                                                    "PROPOSAL::InterpolationSettings::ALLOW_TABLE_CREATION=False");
            }

            if (InterpolationSettings::ALLOW_TABLE_IN_MEMORY) {
                Logging::Get("TableCreation")->warn("PROPOSAL is unable to write to the requested path '{}'. "
                                                    "Therefore, tables will only be stored in memory.", path);
            } else {
                throw std::logic_error("Tables are not available, are not allowed to be created, and are "
                                       "note allows to be stored in memory due to "
                                       "PROPOSAL::InterpolationSettings::ALLOW_TABLE_IN_MEMORY=False. "
                                       "Unable to continue!");
            }
        }
    } else {
        Logging::Get("TableCreation")->debug("Tables are available and are read from file '{}'.", combined);
    }
}

namespace Helper {

    std::string Centered(int width, const std::string& str, char fill)
    {
        int len = str.length();
        if (width < len) {
            return str;
        }

        int diff = width - len;
        int pad1 = diff / 2;
        int pad2 = diff - pad1;
        return std::string(pad1, fill) + str + std::string(pad2, fill);
    }

    bool case_insensitive_comp::operator()(
        std::string const& lhs_ref, std::string const& rhs_ref) const
    {
        auto lhs = std::string(lhs_ref);
        auto rhs = std::string(rhs_ref);
        transform(lhs.begin(), lhs.end(), lhs.begin(), ::tolower);
        transform(rhs.begin(), rhs.end(), rhs.begin(), ::tolower);

        return lhs < rhs;
    }

} // namespace Helper

} // namespace PROPOSAL
