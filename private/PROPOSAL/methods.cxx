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

#include <sys/stat.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <wordexp.h> // Used to expand path with environment variables

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"

#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/Output.h"

// using namespace std;

namespace PROPOSAL
{

float myErfInv2(float x)
{
   float tt1, tt2, lnx, sgn;
   sgn = (x < 0) ? -1.0f : 1.0f;

   x = (1 - x)*(1 + x);        // x = 1 - x*x;
   lnx = logf(x);

   tt1 = 2/(3.141592653589793*0.147) + 0.5f * lnx;
   tt2 = 1/(0.147) * lnx;

   return(sgn*sqrtf(-tt1 + sqrtf(tt1*tt1 - tt2)));
}

bool FileExist(std::string path)
{
    struct stat dummy_stat_return_val;

    if (stat(path.c_str(), &dummy_stat_return_val) != 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}

// round a given double to the closest int
int RoundValue(double val)
{
    bool minus      =   false;
    int valRound    =   0;

    if(val<0)
    {
        val     *=  -1;
        minus   =   true;
    }

    val         +=  0.5;
    valRound    =   (int)val;

    if(minus)
    {
        valRound    *=  -1;
    }

    return valRound;
}

InvalidPath::InvalidPath(const char* pathname, const char* err)
    : pathname_(pathname)
    , msg_()
{
    msg_ = std::string(err) + ": " + pathname_;
}

InvalidPath::InvalidPath(const std::string& pathname, const char* err)
    : pathname_(pathname)
    , msg_()
{
    msg_ = std::string(err) + ": " + pathname_;
}

const char* InvalidPath::what() const throw()
{
    return msg_.c_str();
}

namespace Helper
{

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
bool FileExist(std::string path)
{
    struct stat dummy_stat_return_val;

    if (stat(path.c_str(), &dummy_stat_return_val) != 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}

// ------------------------------------------------------------------------- //
void InitializeInterpolation(const std::string name,
                             InterpolantBuilderContainer& builder_container,
                             const std::vector<Parametrization*>& parametrizations,
                             InterpolationDef interpolation_def)
{
    using namespace std;

    log_info("Initialize %s interpolation.", name.c_str());

    bool storing_failed =   false;
    bool reading_worked =   true;


    // Expand environment variables in path and get full path
    wordexp_t p;
    int success = wordexp(interpolation_def.path_to_tables.c_str(), &p, 0 );

    char full_path[PATH_MAX];
    char* resolved = realpath(*p.we_wordv, full_path);

    wordfree( &p );

    if (!resolved)
    {
        throw InvalidPath(full_path, strerror(errno));
    }

    std::string pathname(resolved);

    bool raw = interpolation_def.raw;

    // --------------------------------------------------------------------- //
    // Create filename out of hash
    // --------------------------------------------------------------------- //

    size_t hash_digest = 0;

    if (parametrizations.size() == 1)
    {
        hash_digest = parametrizations.at(0)->GetHash();
    }
    else
    {
        for(std::vector<Parametrization*>::const_iterator it = parametrizations.begin(); it != parametrizations.end(); ++it)
        {
            boost::hash_combine(hash_digest, (*it)->GetHash());
        }
    }

    stringstream filename;
    filename << pathname << "/" << name << "_" << hash_digest;

    if(!raw)
    {
        filename << ".txt";
    }

    if(!pathname.empty())
    {
        if( FileExist(filename.str()) )
        {
            log_info("%s tables will be read from file: %s", name.c_str(), filename.str().c_str());

            std::ifstream input;

            if(raw)
            {
                input.open(filename.str().c_str(), ios::binary);
            }
            else
            {
                input.open(filename.str().c_str());
            }

            for (InterpolantBuilderContainer::iterator builder_it = builder_container.begin(); builder_it != builder_container.end(); ++builder_it)
            {
                //TODO(mario): read check Tue 2017/09/05
                (*builder_it->second) = new Interpolant();
                (*builder_it->second)->Load(input, raw);

                // for (std::vector<Interpolant*>::iterator interpol_it = builder_it->second.begin(); interpol_it != builder_it->second.end(); ++interpol_it)
                // {
                // }
            }

            input.close();
        }

        if(!FileExist(filename.str()) || !reading_worked )
        {
            if(!reading_worked)
            {
                log_info("file %s is corrupted! Write it again!",filename.str().c_str());
            }

            log_info("%s tables will be saved to file: %s",name.c_str(), filename.str().c_str());

            ofstream output;

            if(raw)
            {
                output.open(filename.str().c_str(), ios::binary);
            }
            else
            {
                output.open(filename.str().c_str());
            }

            if(output.good())
            {
                output.precision(16);

                for (InterpolantBuilderContainer::iterator builder_it = builder_container.begin(); builder_it != builder_container.end(); ++builder_it)
                {
                    (*builder_it->second) = builder_it->first->build();
                    (*builder_it->second)->Save(output, raw);
                }
            }
            else
            {
                storing_failed  =   true;
                log_warn("Can not open file %s for writing! Table will not be stored!",filename.str().c_str());
            }
            output.close();
        }
    }

    if(pathname.empty() || storing_failed)
    {
        log_info("%s tables will be stored in memomy: %s", name.c_str(), filename.str().c_str());

        for (InterpolantBuilderContainer::iterator builder_it = builder_container.begin(); builder_it != builder_container.end(); ++builder_it)
        {
            (*builder_it->second) = builder_it->first->build();
        }
    }

    log_info("Initialize %s interpolation done.", name.c_str());
}

} /* Helper */


}  // PROPOSAL
