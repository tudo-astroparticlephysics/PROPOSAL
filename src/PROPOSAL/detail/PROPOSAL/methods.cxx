/*! \file   methods.cxx
 *   \brief  Source file for the methods routines.
 *
 *   For more details see the class documentation.
 *
 *   \date   21.06.2010
 *   \author Jan-Hendrik Koehne
 */

// #include <stdlib.h>

#include <climits> // for PATH_MAX
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <unistd.h>  // check for write permissions
#include <wordexp.h> // Used to expand path with environment variables

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Logging.h"
#include "PROPOSAL/methods.h"

namespace PROPOSAL {

namespace Helper {

    // -------------------------------------------------------------------------
    // //
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

    // -------------------------------------------------------------------------
    // //
    bool IsWritable(std::string table_dir)
    {
        bool writeable = false;

        if (access(table_dir.c_str(), F_OK) == 0) {
            if ((access(table_dir.c_str(), R_OK) == 0)
                && (access(table_dir.c_str(), W_OK) == 0)) {
                writeable = true;
                Logging::Get("proposal.methods")->debug("Table directory does exist and has read and write "
                          "permissions: %s",
                    table_dir.c_str());
            } else {
                if (access(table_dir.c_str(), R_OK) != 0)
                    Logging::Get("proposal.methods")->warn("Table directory is not readable: %s",
                        table_dir.c_str());
                else
                    Logging::Get("proposal.methods")->warn("Table directory is not writable: %s",
                        table_dir.c_str());
            }
        } else
            Logging::Get("proposal.methods")->warn("Table directory does not exist: %s", table_dir.c_str());

        return writeable;
    }

    // -------------------------------------------------------------------------
    // //
    bool IsReadable(std::string table_dir)
    {
        bool readable = false;

        if (access(table_dir.c_str(), F_OK) == 0) {
            if ((access(table_dir.c_str(), R_OK) == 0)
                && (access(table_dir.c_str(), W_OK) == 0)) {
                readable = true;
                Logging::Get("proposal.methods")->debug("Table directory does exist and has read and write "
                          "permissions: %s",
                    table_dir.c_str());
            } else {
                if (access(table_dir.c_str(), R_OK) == 0) {
                    readable = true;
                    Logging::Get("proposal.methods")->debug("Table directory does exist and has only read "
                              "permissions: "
                              "%s",
                        table_dir.c_str());
                } else
                    Logging::Get("proposal.methods")->warn("Table directory is not readable: %s",
                        table_dir.c_str());
            }
        } else
            Logging::Get("proposal.methods")->warn("Table directory does not exist: %s", table_dir.c_str());

        return readable;
    }

    // -------------------------------------------------------------------------
    // //
    std::string ResolvePath(const std::string& pathname, bool checkReadonly)
    {
        wordexp_t p;
        // Use WRDE_UNDEF to consider undefined shell variables as error
        int success = wordexp(pathname.c_str(), &p, WRDE_UNDEF);

        if (success != 0) {
            return "";
        }

        char full_path[PATH_MAX];
        char* resolved = realpath(*p.we_wordv, full_path);

        wordfree(&p);

        if (!resolved) {
            return "";
        }

        if (checkReadonly) {
            if (IsReadable(std::string(resolved))) {
                return std::string(resolved);
            } else {
                return "";
            }
        }
        if (IsWritable(std::string(resolved))) {
            return std::string(resolved);
        } else {
            return "";
        }
    }

    // -------------------------------------------------------------------------
    // //
    bool FileExist(const std::string path)
    {
        struct stat dummy_stat_return_val;
        if (stat(path.c_str(), &dummy_stat_return_val) != 0) {
            return false;
        } else {
            return true;
        }
    }

    bool case_insensitive_comp::operator()(const std::string &lhs, const std::string &rhs) const {
        return strcasecmp(lhs.c_str(), rhs.c_str()) < 0 ;
    }


} // namespace Helper

} // namespace PROPOSAL
