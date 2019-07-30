/*! \file   methods.cxx
 *   \brief  Source file for the methods routines.
 *
 *   For more details see the class documentation.
 *
 *   \date   21.06.2010
 *   \author Jan-Hendrik Koehne
 */

// #include <stdlib.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <sys/stat.h>
#include <climits> // for PATH_MAX
#include <wordexp.h> // Used to expand path with environment variables
#include <unistd.h> // check for write permissions

#include "PROPOSAL/crossection/parametrization/Parametrization.h"

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Constants.h"

#include "PROPOSAL/Logging.h"
#include "PROPOSAL/methods.h"


namespace PROPOSAL {

double inverseErrorFunction(double p)
{
    if (p <= 0 || p >= 1)
    {
        log_fatal("The inverse Error function can just handle values between 0 and 1.");
    }
    double a_arr[] = {-3.969683028665376e+01,
                      2.209460984245205e+02,
                      -2.759285104469687e+02,
                      1.383577518672690e+02,
                      -3.066479806614716e+01,
                      2.506628277459239e+00};
    double b_arr[] = {-5.447609879822406e+01,
                      1.615858368580409e+02,
                      -1.556989798598866e+02,
                      6.680131188771972e+01,
                      -1.328068155288572e+01};
    double c_arr[] = {-7.784894002430293e-03,
                      -3.223964580411365e-01,
                      -2.400758277161838e+00,
                      -2.549732539343734e+00,
                      4.374664141464968e+00,
                      2.938163982698783e+00};
    double d_arr[] = {7.784695709041462e-03,
                      3.224671290700398e-01,
                      2.445134137142996e+00,
                      3.754408661907416e+00};
    double p_low = 0.02425;
    double p_high = 1 - p_low;

    double q,r,x,e,u;
    if (p < p_low)
    {
        // Rational approximation for lower region.
        q = std::sqrt(-2 * std::log(p));
        x = (((((c_arr[0]*q+c_arr[1])*q+c_arr[2])*q+c_arr[3])*q+c_arr[4])*q+c_arr[5]) /
            ((((d_arr[0]*q+d_arr[1])*q+d_arr[2])*q+d_arr[3])*q+1);
    }
    else if (p_low <= p && p <= p_high)
    {
        // Rational approximation for central region.
        q = p - 0.5;
        r = q*q;
        x = (((((a_arr[0]*r+a_arr[1])*r+a_arr[2])*r+a_arr[3])*r+a_arr[4])*r+a_arr[5])*q /
            (((((b_arr[0]*r+b_arr[1])*r+b_arr[2])*r+b_arr[3])*r+b_arr[4])*r+1);
    }
    else
    {
        // Rational approximation for upper region.
        q = std::sqrt(-2 * std::log(1 - p));
        x = -(((((c_arr[0]*q+c_arr[1])*q+c_arr[2])*q+c_arr[3])*q+c_arr[4])*q+c_arr[5]) /
            ((((d_arr[0]*q+d_arr[1])*q+d_arr[2])*q+d_arr[3])*q+1);
    }

    // Refining the result:
    // One iteration of Halley’s rational method (third order)
    // gives full machine precision.
    e = 0.5 * std::erfc(-x/SQRT2) - p;
    u = e * std::sqrt(2 * PI) * std::exp(0.5 * x * x);
    x = x - u / (1 + x * u * 0.5);
    return x;
}

// ------------------------------------------------------------------------- //
double dilog(double x)
{
    double C_arr[] = { 0.42996693560813697,  0.40975987533077105, -0.01858843665014592, 0.00145751084062268,
                       -0.00014304184442340, 0.00001588415541880, -0.00000190784959387, 0.00000024195180854,
                       -0.00000003193341274, 0.00000000434545063, -0.00000000060578480, 0.00000000008612098,
                       -0.00000000001244332, 0.00000000000182256, -0.00000000000027007, 0.00000000000004042,
                       -0.00000000000000610, 0.00000000000000093, -0.00000000000000014, 0.00000000000000002 };

    double HF = 0.5;

    if (x == 1)
        return PI * PI / 6.0;
    else if (x == -1)
        return -PI * PI / 12.0;
    else
    {
        double _T = -x;
        double Y, S, A, B0, B1, B2, _H, _ALFA;

        if (_T <= -2.0)
        {
            Y  = -1.0 / (1.0 + _T);
            S  = 1.0;
            B1 = std::log(-_T);
            B2 = std::log(1.0 + 1.0 / _T);
            A  = -PI * PI / 3.0 + HF * (B1 * B1 - B2 * B2);
        } else if (_T < -1.0)
        {
            Y = -1.0 - _T;
            S = -1.0;
            A = std::log(-_T);
            A = -PI * PI / 6.0 + A * (A + std::log(1 + 1 / _T));
        } else if (_T <= -0.5)
        {
            Y = -(1 + _T) / _T;
            S = 1.0;
            A = std::log(-_T);
            A = -PI * PI / 6.0 + A * (-HF * A + std::log(1.0 + _T));
        } else if (_T < 0.0)
        {
            Y = -_T / (1 + _T);
            S = -1.0;
            A = std::log(1 + _T);
            A = HF * A * A;
        } else if (_T <= 1.0)
        {
            Y = _T;
            S = 1.0;
            A = 0.0;
        } else
        {
            Y = 1.0 / _T;
            S = -1.0;
            A = std::log(_T);
            A = PI * PI / 6 + HF * A * A;
        }

        _H    = Y + Y - 1;
        _ALFA = _H + _H;
        B1    = 0.0;
        B2    = 0.0;

        for (int i = 19; i >= 0.0; --i)
        {
            B0 = C_arr[i] + _ALFA * B1 - B2;
            B2 = B1;
            B1 = B0;
        }

        return -(S * (B0 - _H * B2) + A);
    }

    return x;
}

// ------------------------------------------------------------------------- //
size_t InterpolationDef::GetHash() const
{
    size_t seed = 0;
    hash_combine(seed,
                 order_of_interpolation,
                 max_node_energy,
                 nodes_cross_section,
                 nodes_continous_randomization,
                 nodes_propagate);

    return seed;
}

namespace Helper {

// ------------------------------------------------------------------------- //
std::string Centered(int width, const std::string& str, char fill)
{
    int len = str.length();
    if (width < len)
    {
        return str;
    }

    int diff = width - len;
    int pad1 = diff / 2;
    int pad2 = diff - pad1;
    return std::string(pad1, fill) + str + std::string(pad2, fill);
}

// ------------------------------------------------------------------------- //
bool IsWritable(std::string table_dir)
{
    bool writeable = false;

    if (access(table_dir.c_str(), F_OK) == 0)
    {
        if ((access(table_dir.c_str(), R_OK) == 0) && (access(table_dir.c_str(), W_OK) == 0))
        {
            writeable = true;
            log_debug("Table directory does exist and has read and write permissions: %s", table_dir.c_str());
        }
        else
        {
            if (access(table_dir.c_str(), R_OK) != 0)
                log_warn("Table directory is not readable: %s", table_dir.c_str());
            else
                log_warn("Table directory is not writable: %s", table_dir.c_str());
        }
    }
    else
        log_warn("Table directory does not exist: %s", table_dir.c_str());

    return writeable;
}


// ------------------------------------------------------------------------- //
bool IsReadable(std::string table_dir)
{
    bool readable = false;

    if (access(table_dir.c_str(), F_OK) == 0)
    {
        if ((access(table_dir.c_str(), R_OK) == 0) && (access(table_dir.c_str(), W_OK) == 0))
        {
            readable = true;
            log_debug("Table directory does exist and has read and write permissions: %s", table_dir.c_str());
        }
        else
        {
            if (access(table_dir.c_str(), R_OK) == 0)
            {
                readable = true;
                log_debug("Table directory does exist and has only read permissions: %s", table_dir.c_str());
            }
            else
                log_warn("Table directory is not readable: %s", table_dir.c_str());
        }
    }
    else
        log_warn("Table directory does not exist: %s", table_dir.c_str());

    return readable;
}


// ------------------------------------------------------------------------- //
std::string ResolvePath(const std::string& pathname, bool checkReadonly)
{
    wordexp_t p;
    // Use WRDE_UNDEF to consider undefined shell variables as error
    int success = wordexp(pathname.c_str(), &p, WRDE_UNDEF);

    if (success != 0)
    {
        log_warn("Invalid path given: \"%s\"", pathname.c_str());
        return "";
    }

    char full_path[PATH_MAX];
    char* resolved = realpath(*p.we_wordv, full_path);

    wordfree(&p);

    if (!resolved)
    {
        log_warn("Invalid path given: \"%s\"", pathname.c_str());
        return "";
    }

    if (checkReadonly)
    {
        if (IsReadable(std::string(resolved)))
        {
            return std::string(resolved);
        }
        else
        {
            return "";
        }
    }
    if (IsWritable(std::string(resolved)))
    {
        return std::string(resolved);
    }
    else
    {
        return "";
    }
}

// ------------------------------------------------------------------------- //
bool FileExist(const std::string path)
{
    struct stat dummy_stat_return_val;

    if (stat(path.c_str(), &dummy_stat_return_val) != 0)
    {
        return false;
    } else
    {
        return true;
    }
}

// ------------------------------------------------------------------------- //
void InitializeInterpolation(const std::string name,
                             InterpolantBuilderContainer& builder_container,
                             const std::vector<Parametrization*>& parametrizations,
                             const InterpolationDef interpolation_def)
{
    log_debug("Initialize %s interpolation.", name.c_str());

    // --------------------------------------------------------------------- //
    // Create hash for the file name
    // --------------------------------------------------------------------- //

    size_t hash_digest = 0;
    if (parametrizations.size() == 1)
    {
        hash_digest = parametrizations[0]->GetHash();
    }
    else
    {
        for (std::vector<Parametrization*>::const_iterator it = parametrizations.begin(); it != parametrizations.end();
             ++it)
        {
            hash_combine(hash_digest, (*it)->GetHash(), (*it)->GetMultiplier(), (*it)->GetParticleDef().low);
        }
        if (name.compare("decay") == 0)
        {
            hash_combine(hash_digest, parametrizations[0]->GetParticleDef().lifetime);
        }
    }
    hash_combine(hash_digest, interpolation_def.GetHash());

    bool storing_failed = false;
    bool reading_worked = false;
    bool binary_tables = interpolation_def.do_binary_tables;
    bool just_use_readonly_path = interpolation_def.just_use_readonly_path;
    std::string pathname;
    std::stringstream filename;

    // --------------------------------------------------------------------- //
    // first check the reading paths
    // if one of the reading paths already has the required tables
    pathname = ResolvePath(interpolation_def.path_to_tables_readonly, true);
    if (!pathname.empty())
    {
        filename << pathname << "/" << name << "_" << hash_digest;
        if (!binary_tables)
        {
            filename << ".txt";
        }
        if (FileExist(filename.str()))
        {
            std::ifstream input;
            if (binary_tables)
            {
                input.open(filename.str().c_str(), std::ios::binary);
            } else
            {
                input.open(filename.str().c_str());
            }

            // check if file is empty
            // this happens if multiple instances tries to load/create the tables in parallel
            // and another process already starts to write this table
            // now just hand over to writing process where it might saves them in memory
            // if the other instance is still writing them down in the same path
            if (input.peek() == std::ifstream::traits_type::eof())
            {
                log_info("file %s is empty! Another process is presumably writing. Try another reading path or write in memory!",
                        filename.str().c_str());
            }
            else
            {
                log_debug("%s tables will be read from file: %s", name.c_str(), filename.str().c_str());

                for (InterpolantBuilderContainer::iterator builder_it = builder_container.begin();
                     builder_it != builder_container.end();
                     ++builder_it)
                {
                    // TODO(mario): read check Tue 2017/09/05
                    (*builder_it->second) = new Interpolant();
                    (*builder_it->second)->Load(input, binary_tables);
                }
                reading_worked = true;
            }

            input.close();

        }
        else
        {
            log_debug("In the readonly path to the interpolation tables, the file %s does not Exist",
                      filename.str().c_str());
        }
    }
    else
    {
        log_debug("No reading path was given, now the tables are read or written to the writing path.");
    }

    if (reading_worked)
    {
        log_debug("Initialize %s interpolation done.", name.c_str());
        return;
    }

    if (just_use_readonly_path)
    {
        log_fatal("The just_use_readonly_path option is enabled and the table is not in the readonly path.");
    }

    // --------------------------------------------------------------------- //
    // if none of the reading paths has the required interpolation table
    // the interpolation tables will be written in the path for writing
    pathname = ResolvePath(interpolation_def.path_to_tables);

    // clear the stringstream
    filename.str(std::string());
    filename.clear();
    filename << pathname << "/" << name << "_" << hash_digest;

    if (!binary_tables)
    {
        filename << ".txt";
    }

    if (!pathname.empty())
    {
        if (FileExist(filename.str()))
        {
            std::ifstream input;

            if (binary_tables)
            {
                input.open(filename.str().c_str(), std::ios::binary);
            } else
            {
                input.open(filename.str().c_str());
            }

            // check if file is empty
            // this happens if multiple instances try to write the tables in parallel
            // now just one is writing them and the other just saves them in memory
            if (input.peek() == std::ifstream::traits_type::eof())
            {
                log_info("file %s is empty! Another process is presumably writing. Save this table in memory!",
                        filename.str().c_str());
                storing_failed = true;
            }
            else
            {
                log_debug("%s tables will be read from file: %s", name.c_str(), filename.str().c_str());

                for (InterpolantBuilderContainer::iterator builder_it = builder_container.begin();
                     builder_it != builder_container.end();
                     ++builder_it)
                {
                    // TODO(mario): read check Tue 2017/09/05
                    (*builder_it->second) = new Interpolant();
                    (*builder_it->second)->Load(input, binary_tables);
                }
            }

            input.close();
        }
        else
        {
            log_debug("%s tables will be saved to file: %s", name.c_str(), filename.str().c_str());

            std::ofstream output;

            if (binary_tables)
            {
                output.open(filename.str().c_str(), std::ios::binary);
            } else
            {
                output.open(filename.str().c_str());
            }

            if (output.good())
            {
                output.precision(16);

                for (InterpolantBuilderContainer::iterator builder_it = builder_container.begin();
                     builder_it != builder_container.end();
                     ++builder_it)
                {
                    (*builder_it->second) = builder_it->first->build();
                    (*builder_it->second)->Save(output, binary_tables);
                }
            } else
            {
                storing_failed = true;
                log_warn("Can not open file %s for writing! Table will not be stored!", filename.str().c_str());
            }
            output.close();
        }
    }

    if (pathname.empty() || storing_failed)
    {
        log_debug("%s tables will be stored in memomy!", name.c_str());

        for (InterpolantBuilderContainer::iterator builder_it = builder_container.begin();
             builder_it != builder_container.end();
             ++builder_it)
        {
            (*builder_it->second) = builder_it->first->build();
        }
    }

    log_debug("Initialize %s interpolation done.", name.c_str());
}

} // namespace Helper

} // namespace PROPOSAL
