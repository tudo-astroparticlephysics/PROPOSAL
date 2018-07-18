/*! \file   methods.cxx
 *   \brief  Source file for the methods routines.
 *
 *   For more details see the class documentation.
 *
 *   \date   21.06.2010
 *   \author Jan-Hendrik Koehne
 */

// #include <cmath>
// #include <stdlib.h>
#include <boost/functional/hash.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <wordexp.h> // Used to expand path with environment variables
#include <unistd.h> // check for write permissions

#include "PROPOSAL/crossection/parametrization/Parametrization.h"

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Constants.h"

#include "PROPOSAL/Output.h"
#include "PROPOSAL/methods.h"

// using namespace std;

namespace PROPOSAL {

float myErfInv2(float x)
{
    float tt1, tt2, lnx, sgn;
    sgn = (x < 0) ? -1.0f : 1.0f;

    x   = (1 - x) * (1 + x); // x = 1 - x*x;
    lnx = logf(x);

    tt1 = 2 / (3.141592653589793 * 0.147) + 0.5f * lnx;
    tt2 = 1 / (0.147) * lnx;

    return (sgn * sqrtf(-tt1 + sqrtf(tt1 * tt1 - tt2)));
}

// round a given double to the closest int
int RoundValue(double val)
{
    bool minus   = false;
    int valRound = 0;

    if (val < 0)
    {
        val *= -1;
        minus = true;
    }

    val += 0.5;
    valRound = (int)val;

    if (minus)
    {
        valRound *= -1;
    }

    return valRound;
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
            B1 = log(-_T);
            B2 = log(1.0 + 1.0 / _T);
            A  = -PI * PI / 3.0 + HF * (B1 * B1 - B2 * B2);
        } else if (_T < -1.0)
        {
            Y = -1.0 - _T;
            S = -1.0;
            A = log(-_T);
            A = -PI * PI / 6.0 + A * (A + log(1 + 1 / _T));
        } else if (_T <= -0.5)
        {
            Y = -(1 + _T) / _T;
            S = 1.0;
            A = log(-_T);
            A = -PI * PI / 6.0 + A * (-HF * A + log(1.0 + _T));
        } else if (_T < 0.0)
        {
            Y = -_T / (1 + _T);
            S = -1.0;
            A = log(1 + _T);
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
            A = log(_T);
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
            log_info("Table directory does exist and has read and write permissions: %s", table_dir.c_str());
        }
        else
        {
            if (access(table_dir.c_str(), R_OK) != 0)
                log_info("Table directory is not readable: %s", table_dir.c_str());
            else
                log_info("Table directory is not writable: %s", table_dir.c_str());
        }
    }
    else
        log_info("Table directory does not exist: %s", table_dir.c_str());

    return writeable;
}

// ------------------------------------------------------------------------- //
std::string ResolvePath(const std::string& pathname)
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
    using namespace std;

    log_info("Initialize %s interpolation.", name.c_str());

    bool storing_failed = false;
    bool reading_worked = true;

    std::string pathname = ResolvePath(interpolation_def.path_to_tables);
    bool raw             = interpolation_def.raw;

    // --------------------------------------------------------------------- //
    // Create filename out of hash
    // --------------------------------------------------------------------- //

    size_t hash_digest = 0;

    if (parametrizations.size() == 1)
    {
        hash_digest = parametrizations.at(0)->GetHash();
    } else
    {
        for (std::vector<Parametrization*>::const_iterator it = parametrizations.begin(); it != parametrizations.end();
             ++it)
        {
            boost::hash_combine(hash_digest, (*it)->GetHash());
        }
    }

    stringstream filename;
    filename << pathname << "/" << name << "_" << hash_digest;

    if (!raw)
    {
        filename << ".txt";
    }

    if (!pathname.empty())
    {
        if (FileExist(filename.str()))
        {
            log_info("%s tables will be read from file: %s", name.c_str(), filename.str().c_str());

            std::ifstream input;

            if (raw)
            {
                input.open(filename.str().c_str(), ios::binary);
            } else
            {
                input.open(filename.str().c_str());
            }

            for (InterpolantBuilderContainer::iterator builder_it = builder_container.begin();
                 builder_it != builder_container.end();
                 ++builder_it)
            {
                // TODO(mario): read check Tue 2017/09/05
                (*builder_it->second) = new Interpolant();
                (*builder_it->second)->Load(input, raw);
            }

            input.close();
        }

        if (!FileExist(filename.str()) || !reading_worked)
        {
            if (!reading_worked)
            {
                log_info("file %s is corrupted! Write it again!", filename.str().c_str());
            }

            log_info("%s tables will be saved to file: %s", name.c_str(), filename.str().c_str());

            ofstream output;

            if (raw)
            {
                output.open(filename.str().c_str(), ios::binary);
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
                    (*builder_it->second)->Save(output, raw);
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
        log_info("%s tables will be stored in memomy!", name.c_str());

        for (InterpolantBuilderContainer::iterator builder_it = builder_container.begin();
             builder_it != builder_container.end();
             ++builder_it)
        {
            (*builder_it->second) = builder_it->first->build();
        }
    }

    log_info("Initialize %s interpolation done.", name.c_str());
}

} // namespace Helper

} // namespace PROPOSAL
