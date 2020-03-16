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

#include "PROPOSAL/crossection/parametrization/Parametrization.h"

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Logging.h"
#include "PROPOSAL/methods.h"

namespace PROPOSAL {

InterpolationDef::InterpolationDef(const nlohmann::json& config)
{
    nodes_propagate = config.value("nodes_propagate", 1000);
    nodes_continous_randomization
        = config.value("nodes_continous_randomization", 200);
    nodes_cross_section = config.value("nodes_cross_section", 100);
    max_node_energy = config.value("max_node_energy", 1e14);
    do_binary_tables = config.value("do_binary_tables", true);
    just_use_readonly_path = config.value("just_use_readonly_path", false);
    order_of_interpolation = config.value("order_of_interpolation", 5);

    if (!(nodes_propagate > 3))
        throw std::invalid_argument(
            "At least 3 nodes are required for qubic splines");
    if (!(nodes_continous_randomization > 4))
        throw std::invalid_argument(
            "At least 3 nodes are required for qubic splines");
    if (!(nodes_cross_section > 4))
        throw std::invalid_argument(
            "At least 3 nodes are required for qubic splines");
    if (!(max_node_energy > 0))
        throw std::invalid_argument("max_node_energy must be larger than "
                                    "highest primary particle energy.");
    if (!(order_of_interpolation > 1))
        throw std::invalid_argument(
            "Order of interpolation must be larger than one.");

    if (not config.contains("path_to_tables")) {
        log_warn("No valid writable path to interpolation tables found. Save "
                 "tables in memory, if readonly path is also not working!");
    } else {
        if (config.at("path_to_tables").is_string())
            path_to_tables = Helper::ResolvePath(config.at("path_to_tables"));

        if (config.at("path_to_tables").is_array()) {
            for (const auto& path : config.at("path_to_tables")) {
                if (path.is_string()) {
                    path_to_tables = Helper::ResolvePath(path);
                    if (path_to_tables != "")
                        break;
                }
            }
        }
        if (path_to_tables == "")
            throw std::invalid_argument(
                "Invalid input for option 'path_to_tables'. Expected a "
                "string or a list of strings.");
    };

    if (not config.contains("path_to_tables_readonly")) {
        log_warn("No valid writable path to interpolation tables readonly found."
                 "Save tables in memory, if readonly path is also not working!");
    } else {
        if (config.at("path_to_tables_readonly").is_string())
            path_to_tables_readonly = Helper::ResolvePath(config.at("path_to_tables_readonly"));

        if (config.at("path_to_tables_readonly").is_array()) {
            for (const auto& path : config.at("path_to_tables_readonly")) {
                if (path.is_string()) {
                    path_to_tables_readonly = Helper::ResolvePath(path);
                    if (path_to_tables_readonly != "")
                        break;
                }
            }
        }
        if (path_to_tables_readonly == "")
            throw std::invalid_argument(
                "Invalid input for option 'path_to_tables_readonly'. Expected a "
                "string or a list of strings.");
    };

}

// ------------------------------------------------------------------------- //
size_t InterpolationDef::GetHash() const
{
    size_t seed = 0;
    hash_combine(seed, order_of_interpolation, max_node_energy,
        nodes_cross_section, nodes_continous_randomization, nodes_propagate);

    return seed;
}

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
                log_debug("Table directory does exist and has read and write "
                          "permissions: %s",
                    table_dir.c_str());
            } else {
                if (access(table_dir.c_str(), R_OK) != 0)
                    log_warn("Table directory is not readable: %s",
                        table_dir.c_str());
                else
                    log_warn("Table directory is not writable: %s",
                        table_dir.c_str());
            }
        } else
            log_warn("Table directory does not exist: %s", table_dir.c_str());

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
                log_debug("Table directory does exist and has read and write "
                          "permissions: %s",
                    table_dir.c_str());
            } else {
                if (access(table_dir.c_str(), R_OK) == 0) {
                    readable = true;
                    log_debug("Table directory does exist and has only read "
                              "permissions: "
                              "%s",
                        table_dir.c_str());
                } else
                    log_warn("Table directory is not readable: %s",
                        table_dir.c_str());
            }
        } else
            log_warn("Table directory does not exist: %s", table_dir.c_str());

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

    // -------------------------------------------------------------------------
    // //
    void InitializeInterpolation(const std::string name,
        InterpolantBuilderContainer& builder_container,
        const std::vector<Parametrization*>& parametrizations,
        const InterpolationDef interpolation_def)
    {
        log_debug("Initialize %s interpolation.", name.c_str());

        // ---------------------------------------------------------------------
        // // Create hash for the file name
        // ---------------------------------------------------------------------
        // //

        size_t hash_digest = 0;
        if (parametrizations.size() == 1) {
            hash_digest = parametrizations[0]->GetHash();
        } else {
            for (std::vector<Parametrization*>::const_iterator it
                 = parametrizations.begin();
                 it != parametrizations.end(); ++it) {
                hash_combine(hash_digest, (*it)->GetHash(),
                    (*it)->GetMultiplier(), (*it)->GetParticleDef().low);
            }
            if (name.compare("decay") == 0) {
                hash_combine(hash_digest,
                    parametrizations[0]->GetParticleDef().lifetime);
            }
        }
        hash_combine(hash_digest, interpolation_def.GetHash());

        bool storing_failed = false;
        bool reading_worked = false;
        bool binary_tables = interpolation_def.do_binary_tables;
        bool just_use_readonly_path = interpolation_def.just_use_readonly_path;
        std::string pathname;
        std::stringstream filename;

        // ---------------------------------------------------------------------
        // // first check the reading paths if one of the reading paths already
        // has the required tables
        pathname = ResolvePath(interpolation_def.path_to_tables_readonly, true);
        if (!pathname.empty()) {
            filename << pathname << "/" << name << "_" << hash_digest;
            if (!binary_tables) {
                filename << ".txt";
            }
            if (FileExist(filename.str())) {
                std::ifstream input;
                if (binary_tables) {
                    input.open(filename.str().c_str(), std::ios::binary);
                } else {
                    input.open(filename.str().c_str());
                }

                // check if file is empty
                // this happens if multiple instances tries to load/create the
                // tables in parallel and another process already starts to
                // write this table now just hand over to writing process where
                // it might saves them in memory if the other instance is still
                // writing them down in the same path
                if (input.peek() == std::ifstream::traits_type::eof()) {
                    log_info("file %s is empty! Another process is presumably "
                             "writing. "
                             "Try another reading path or write in memory!",
                        filename.str().c_str());
                } else {
                    log_debug("%s tables will be read from file: %s",
                        name.c_str(), filename.str().c_str());

                    for (InterpolantBuilderContainer::iterator builder_it
                         = builder_container.begin();
                         builder_it != builder_container.end(); ++builder_it) {
                        // TODO(mario): read check Tue 2017/09/05
                        (*builder_it->second) = new Interpolant();
                        (*builder_it->second)->Load(input, binary_tables);
                    }
                    reading_worked = true;
                }

                input.close();

            } else {
                log_debug("In the readonly path to the interpolation tables, "
                          "the file %s "
                          "does not Exist",
                    filename.str().c_str());
            }
        } else {
            log_debug("No reading path was given, now the tables are read or "
                      "written to "
                      "the writing path.");
        }

        if (reading_worked) {
            log_debug("Initialize %s interpolation done.", name.c_str());
            return;
        }

        if (just_use_readonly_path) {
            log_fatal("The just_use_readonly_path option is enabled and the "
                      "table is not "
                      "in the readonly path.");
        }

        // ---------------------------------------------------------------------
        // // if none of the reading paths has the required interpolation table
        // the interpolation tables will be written in the path for writing
        pathname = ResolvePath(interpolation_def.path_to_tables);

        // clear the stringstream
        filename.str(std::string());
        filename.clear();
        filename << pathname << "/" << name << "_" << hash_digest;

        if (!binary_tables) {
            filename << ".txt";
        }

        if (!pathname.empty()) {
            if (FileExist(filename.str())) {
                std::ifstream input;

                if (binary_tables) {
                    input.open(filename.str().c_str(), std::ios::binary);
                } else {
                    input.open(filename.str().c_str());
                }

                // check if file is empty
                // this happens if multiple instances try to write the tables in
                // parallel now just one is writing them and the other just
                // saves them in memory
                if (input.peek() == std::ifstream::traits_type::eof()) {
                    log_info("file %s is empty! Another process is presumably "
                             "writing. "
                             "Save this table in memory!",
                        filename.str().c_str());
                    storing_failed = true;
                } else {
                    log_debug("%s tables will be read from file: %s",
                        name.c_str(), filename.str().c_str());

                    for (InterpolantBuilderContainer::iterator builder_it
                         = builder_container.begin();
                         builder_it != builder_container.end(); ++builder_it) {
                        // TODO(mario): read check Tue 2017/09/05
                        (*builder_it->second) = new Interpolant();
                        (*builder_it->second)->Load(input, binary_tables);
                    }
                }

                input.close();
            } else {
                log_debug("%s tables will be saved to file: %s", name.c_str(),
                    filename.str().c_str());

                std::ofstream output;

                if (binary_tables) {
                    output.open(filename.str().c_str(), std::ios::binary);
                } else {
                    output.open(filename.str().c_str());
                }

                if (output.good()) {
                    output.precision(16);

                    for (InterpolantBuilderContainer::iterator builder_it
                         = builder_container.begin();
                         builder_it != builder_container.end(); ++builder_it) {
                        (*builder_it->second) = builder_it->first->build();
                        (*builder_it->second)->Save(output, binary_tables);
                    }
                } else {
                    storing_failed = true;
                    log_warn(
                        "Can not open file %s for writing! Table will not be "
                        "stored!",
                        filename.str().c_str());
                }
                output.close();
            }
        }

        if (pathname.empty() || storing_failed) {
            log_debug("%s tables will be stored in memomy!", name.c_str());

            for (InterpolantBuilderContainer::iterator builder_it
                 = builder_container.begin();
                 builder_it != builder_container.end(); ++builder_it) {
                (*builder_it->second) = builder_it->first->build();
            }
        }

        log_debug("Initialize %s interpolation done.", name.c_str());
    }

} // namespace Helper

} // namespace PROPOSAL
