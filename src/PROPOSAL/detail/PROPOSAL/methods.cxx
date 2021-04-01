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

std::string InterpolationDef::path_to_tables = "";
std::string InterpolationDef::path_to_tables_readonly = "";
bool InterpolationDef::do_binary_tables = true;
bool InterpolationDef::just_use_readonly_path = false;

InterpolationDef::InterpolationDef(const nlohmann::json& config)
{
    do_binary_tables = config.value("do_binary_tables", true);
    just_use_readonly_path = config.value("just_use_readonly_path", false);

    if (not config.contains("path_to_tables")) {
        Logging::Get("proposal.methods")->warn("No valid writable path to interpolation tables found. Save "
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
        Logging::Get("proposal.methods")->warn(
            "No valid writable path to interpolation tables readonly found."
            "Save tables in memory, if readonly path is also not working!");
    } else {
        if (config.at("path_to_tables_readonly").is_string())
            path_to_tables
                = Helper::ResolvePath(config.at("path_to_tables_readonly"));

        if (config.at("path_to_tables_readonly").is_array()) {
            for (const auto& path : config.at("path_to_tables_readonly")) {
                if (path.is_string()) {
                    path_to_tables_readonly = Helper::ResolvePath(path);
                    if (path_to_tables_readonly != "")
                        break;
                }
            }
        }
        if (path_to_tables == "")
            std::invalid_argument("Invalid input for option "
                                  "'path_to_tables_readonly'. Expected a "
                                  "string or a list of strings.");
    };
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

    // -------------------------------------------------------------------------
    // //
    /* std::unique_ptr<Interpolant> InitializeInterpolation(std::string name, */
    /*     unique_ptr<InterpolantBuilder> builder, size_t hash_digest, */
    /*     const InterpolationDef& interpolation_def) */
    /* { */
    /*     // Simple wrapper for inizializing one Interpolant only */
    /*     Helper::InterpolantBuilderContainer builder_container; */
    /*     builder_container.push_back(std::move(builder)); */
    /*     auto return_vec = InitializeInterpolation( */
    /*         name, builder_container, hash_digest, interpolation_def); */
    /*     return std::move(return_vec.at(0)); */
    /* } */

    std::unique_ptr<Interpolant> InitializeInterpolation( std::string name,
            InterpolantBuilder const& builder, size_t hash_digest)
    {

        bool storing_failed = false;
        std::stringstream filename;

        // ---------------------------------------------------------------------
        // // first check the reading paths if one of the reading paths already
        // has the required tables

        auto pathname = ResolvePath(InterpolationDef::path_to_tables_readonly, true);
        if (!pathname.empty()) {
            filename << pathname << "/" << name << "_" << hash_digest;
            if (!InterpolationDef::do_binary_tables) {
                filename << ".txt";
            }
            if (FileExist(filename.str())) {
                std::ifstream input;
                if (InterpolationDef::do_binary_tables) {
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
                    Logging::Get("proposal.methods")->info("file %s is empty! Another process is presumably "
                             "writing. "
                             "Try another reading path or write in memory!",
                        filename.str().c_str());
                }
                Logging::Get("proposal.methods")->debug("%s tables will be read from file: %s", name.c_str(), filename.str().c_str());

                auto interpolant = std::make_unique<Interpolant>();
                if(interpolant->Load(input, InterpolationDef::do_binary_tables))
                    return interpolant;

            } else {
                Logging::Get("proposal.methods")->debug("In the readonly path to the interpolation tables, "
                          "the file %s "
                          "does not Exist",
                    filename.str().c_str());
            }
        } else {
            Logging::Get("proposal.methods")->debug("No reading path was given, now the tables are read or "
                      "written to "
                      "the writing path.");
        }

        if (InterpolationDef::just_use_readonly_path) {
            Logging::Get("proposal.methods")->critical("The just_use_readonly_path option is enabled and the "
                      "table is not "
                      "in the readonly path.");
        }

        // ---------------------------------------------------------------------
        // // if none of the reading paths has the required interpolation table
        // the interpolation tables will be written in the path for writing

        pathname = ResolvePath(InterpolationDef::path_to_tables);

        // clear the stringstream
        filename.str(std::string());
        filename.clear();
        filename << pathname << "/" << name << "_" << hash_digest;

        if (!InterpolationDef::do_binary_tables) {
            filename << ".txt";
        }

        if (!pathname.empty()) {
            if (FileExist(filename.str())) {
                std::ifstream input;

                if (InterpolationDef::do_binary_tables) {
                    input.open(filename.str().c_str(), std::ios::binary);
                } else {
                    input.open(filename.str().c_str());
                }

                // check if file is empty
                // this happens if multiple instances try to write the tables in
                // parallel now just one is writing them and the other just
                // saves them in memory
                if (input.peek() == std::ifstream::traits_type::eof()) {
                    Logging::Get("proposal.methods")->info("file %s is empty! Another process is presumably "
                             "writing. "
                             "Save this table in memory!",
                        filename.str().c_str());
                    storing_failed = true;
                } else {
                    Logging::Get("proposal.methods")->debug("%s tables will be read from file: %s",
                        name.c_str(), filename.str().c_str());

                    auto interpol = make_unique<Interpolant>();
                    if(interpol->Load(input, InterpolationDef::do_binary_tables))
                        return interpol;
                }

                input.close();
            } else {
                Logging::Get("proposal.methods")->debug("%s tables will be saved to file: %s", name.c_str(),
                    filename.str().c_str());

                std::ofstream output;

                if (InterpolationDef::do_binary_tables) {
                    output.open(filename.str().c_str(), std::ios::binary);
                } else {
                    output.open(filename.str().c_str());
                }

                if (output.good()) {
                    output.precision(16);

                    auto interpolant = builder.build();
                    if(interpolant->Save(output, InterpolationDef::do_binary_tables))
                        return interpolant;

                } else {
                    storing_failed = true;
                    Logging::Get("proposal.methods")->warn(
                        "Can not open file %s for writing! Table will not be "
                        "stored!",
                        filename.str().c_str());
                }
                output.close();
            }
        }

        if (pathname.empty() || storing_failed) {
            Logging::Get("proposal.methods")->debug("%s tables will be stored in memomy!", name.c_str());

            return builder.build();
        }

        throw std::logic_error("Never get here");
    }

    bool case_insensitive_comp::operator()(const std::string &lhs, const std::string &rhs) const {
        return strcasecmp(lhs.c_str(), rhs.c_str()) < 0 ;
    }

    /* std::vector<std::unique_ptr<Interpolant>> InitializeInterpolation( */
    /*     std::string name, const InterpolantBuilderContainer& builder_container, */
    /*     size_t hash_digest, const InterpolationDef& interpolation_def) */
    /* { */
    /*     std::vector<std::unique_ptr<Interpolant>> interpolants; */

    /*     bool storing_failed = false; */
    /*     bool reading_worked = false; */
    /*     bool binary_tables = interpolation_def.do_binary_tables; */
    /*     bool just_use_readonly_path = interpolation_def.just_use_readonly_path; */
    /*     std::string pathname; */
    /*     std::stringstream filename; */

    /*     // --------------------------------------------------------------------- */
    /*     // // first check the reading paths if one of the reading paths already */
    /*     // has the required tables */

    /*     pathname = ResolvePath(InterpolationDef::path_to_tables_readonly, true); */
    /*     if (!pathname.empty()) { */
    /*         filename << pathname << "/" << name << "_" << hash_digest; */
    /*         if (!binary_tables) { */
    /*             filename << ".txt"; */
    /*         } */
    /*         if (FileExist(filename.str())) { */
    /*             std::ifstream input; */
    /*             if (binary_tables) { */
    /*                 input.open(filename.str().c_str(), std::ios::binary); */
    /*             } else { */
    /*                 input.open(filename.str().c_str()); */
    /*             } */

    /*             // check if file is empty */
    /*             // this happens if multiple instances tries to load/create the */
    /*             // tables in parallel and another process already starts to */
    /*             // write this table now just hand over to writing process where */
    /*             // it might saves them in memory if the other instance is still */
    /*             // writing them down in the same path */
    /*             if (input.peek() == std::ifstream::traits_type::eof()) { */
    /*                 Logging::Get("proposal.methods")->info("file %s is empty! Another process is presumably " */
    /*                          "writing. " */
    /*                          "Try another reading path or write in memory!", */
    /*                     filename.str().c_str()); */
    /*             } else { */
    /*                 Logging::Get("proposal.methods")->debug("%s tables will be read from file: %s", */
    /*                     name.c_str(), filename.str().c_str()); */

    /*                 for (size_t i = 0; i < builder_container.size(); ++i) { */
    /*                     // TODO(mario): read check Tue 2017/09/05 */
    /*                     interpolants.emplace_back( */
    /*                         unique_ptr<Interpolant>(new Interpolant())); */
    /*                     interpolants.back()->Load(input, binary_tables); */
    /*                 } */
    /*                 reading_worked = true; */
    /*             } */

    /*             input.close(); */

    /*         } else { */
    /*             Logging::Get("proposal.methods")->debug("In the readonly path to the interpolation tables, " */
    /*                       "the file %s " */
    /*                       "does not Exist", */
    /*                 filename.str().c_str()); */
    /*         } */
    /*     } else { */
    /*         Logging::Get("proposal.methods")->debug("No reading path was given, now the tables are read or " */
    /*                   "written to " */
    /*                   "the writing path."); */
    /*     } */

    /*     if (reading_worked) { */
    /*         Logging::Get("proposal.methods")->debug("Initialize %s interpolation done.", name.c_str()); */
    /*         return interpolants; */
    /*     } */

    /*     if (just_use_readonly_path) { */
    /*         Logging::Get("proposal.methods")->critical("The just_use_readonly_path option is enabled and the " */
    /*                   "table is not " */
    /*                   "in the readonly path."); */
    /*     } */

    /*     // --------------------------------------------------------------------- */
    /*     // // if none of the reading paths has the required interpolation table */
    /*     // the interpolation tables will be written in the path for writing */

    /*     pathname = ResolvePath(InterpolationDef::path_to_tables); */

    /*     // clear the stringstream */
    /*     filename.str(std::string()); */
    /*     filename.clear(); */
    /*     filename << pathname << "/" << name << "_" << hash_digest; */

    /*     if (!binary_tables) { */
    /*         filename << ".txt"; */
    /*     } */

    /*     if (!pathname.empty()) { */
    /*         if (FileExist(filename.str())) { */
    /*             std::ifstream input; */

    /*             if (binary_tables) { */
    /*                 input.open(filename.str().c_str(), std::ios::binary); */
    /*             } else { */
    /*                 input.open(filename.str().c_str()); */
    /*             } */

    /*             // check if file is empty */
    /*             // this happens if multiple instances try to write the tables in */
    /*             // parallel now just one is writing them and the other just */
    /*             // saves them in memory */
    /*             if (input.peek() == std::ifstream::traits_type::eof()) { */
    /*                 Logging::Get("proposal.methods")->info("file %s is empty! Another process is presumably " */
    /*                          "writing. " */
    /*                          "Save this table in memory!", */
    /*                     filename.str().c_str()); */
    /*                 storing_failed = true; */
    /*             } else { */
    /*                 Logging::Get("proposal.methods")->debug("%s tables will be read from file: %s", */
    /*                     name.c_str(), filename.str().c_str()); */

    /*                 for (size_t i=0; i<builder_container.size(); ++i) { */
    /*                     // TODO(mario): read check Tue 2017/09/05 */
    /*                     interpolants.emplace_back(new Interpolant()); */
    /*                     interpolants.back()->Load(input, binary_tables); */
    /*                 } */
    /*             } */

    /*             input.close(); */
    /*         } else { */
    /*             Logging::Get("proposal.methods")->debug("%s tables will be saved to file: %s", name.c_str(), */
    /*                 filename.str().c_str()); */

    /*             std::ofstream output; */

    /*             if (binary_tables) { */
    /*                 output.open(filename.str().c_str(), std::ios::binary); */
    /*             } else { */
    /*                 output.open(filename.str().c_str()); */
    /*             } */

    /*             if (output.good()) { */
    /*                 output.precision(16); */

    /*                 for (const auto& builder : builder_container) { */
    /*                     interpolants.emplace_back(builder->build()); */
    /*                     interpolants.back()->Save(output, binary_tables); */
    /*                 } */
    /*             } else { */
    /*                 storing_failed = true; */
    /*                 Logging::Get("proposal.methods")->warn( */
    /*                     "Can not open file %s for writing! Table will not be " */
    /*                     "stored!", */
    /*                     filename.str().c_str()); */
    /*             } */
    /*             output.close(); */
    /*         } */
    /*     } */

    /*     if (pathname.empty() || storing_failed) { */
    /*         Logging::Get("proposal.methods")->debug("%s tables will be stored in memomy!", name.c_str()); */

    /*         for (const auto& builder : builder_container) */
    /*             interpolants.emplace_back(builder->build()); */
    /*     } */

    /*     Logging::Get("proposal.methods")->debug("Initialize %s interpolation done.", name.c_str()); */
    /*     return interpolants; */
    /* } */

} // namespace Helper

} // namespace PROPOSAL
