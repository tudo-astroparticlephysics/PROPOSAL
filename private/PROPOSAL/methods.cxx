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

// #include <boost/random.hpp>
// #include <boost/generator_iterator.hpp>

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"

//TODO(mario): remove Thu 2017/09/07
// #include "PROPOSAL/crossection/parametrization/Ionization.h"
// #include "PROPOSAL/crossection/parametrization/Bremsstrahlung.h"
// #include "PROPOSAL/crossection/parametrization/EpairProduction.h"
// #include "PROPOSAL/crossection/parametrization/PhotoRealPhotonAssumption.h"

#include "PROPOSAL/math/Interpolant.h"
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

// double RationalApproximation(double t)
// {
//     // Abramowitz and Stegun formula 26.2.23.
//     // The absolute value of the error should be less than 4.5 e-4.
//     double c[] = {2.515517, 0.802853, 0.010328};
//     double d[] = {1.432788, 0.189269, 0.001308};
//     return t - ((c[2]*t + c[1])*t + c[0]) /
//                (((d[2]*t + d[1])*t + d[0])*t + 1.0);
// }
//
// double NormalCDFInverse(double p)
// {
//     if (p <= 0.0 || p >= 1.0)
//     {
//         std::stringstream os;
//         os << "Invalid input argument (" << p
//            << "); must be larger than 0 but less than 1.";
//         throw std::invalid_argument( os.str() );
//     }
//
//     // See article above for explanation of this section.
//     if (p < 0.5)
//     {
//         // F^-1(p) = - G^-1(p)
//         return -RationalApproximation( sqrt(-2.0*log(p)) );
//     }
//     else
//     {
//         // F^-1(p) = G^-1(1-p)
//         return RationalApproximation( sqrt(-2.0*log(1-p)) );
//     }
// }

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


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// returns true if text starts with token, false otherwise

// bool StartsWith(const std::string& text,const std::string& token)
// {

//     if(text.length() < token.length())
//     {
//         return false;
//     }
// 	return (text.compare(0, token.length(), token) == 0);
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// bool EndsWith(const std::string& text,const std::string& token)
// {
// 	return (text.compare(text.length()-1,token.length(),token) ==0);
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


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


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



// string ToLowerCase(string toConvert)
// {
//     string buffer;

//     for(int i=0; i<(int)toConvert.length(); i++)
//     {
//         buffer  +=  (char)tolower(toConvert.at(i));
//     }


//     return buffer;
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// string ReplaceAll(string toConvert, const char oldChar, const char newChar)
// {
//     string buffer="";

//     for(int i =0; i<(int)toConvert.size(); i++)
//     {
//         if(toConvert.at(i)==oldChar)
//         {
//             buffer  +=  newChar;
//         }
//         else
//         {
//             buffer  +=  toConvert.at(i);
//         }
//     }

//     return buffer;
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// returns a random Double

// double Old_RandomDouble()
// {
// 	double result;

//     result  =   rand() + 0.0;
//     result  =   result / RAND_MAX;

// 	return result;
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// deque<string>* SplitString(string args, string Delimiters)
// {

//     deque<string> *Tokens   =   new deque<string>;
//     string::size_type Next, EndNext =   0;

//     while (EndNext != string::npos)
//     {
//         // Find next token
//         Next    =   args.find_first_not_of(Delimiters, EndNext);
//         EndNext =   args.find_first_of(Delimiters, Next);

//         // Stop if end of string reached
//         if (Next == string::npos)
//         {
//             break;
//         }

//         // Add token to vector.
//         Tokens->push_back(args.substr(Next, EndNext - Next));
//     }

//     return Tokens;
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// string NextToken(deque<string> *Tokens)
// {
//     string next;
//     next    =   Tokens->front();
//     Tokens->pop_front();

//     return next;

// }


// ----------------------------------------------------------------------------
/// @brief Helper to initialize Interpolants with interpolation tables
// ----------------------------------------------------------------------------
// class InitInterpolationHelper
// {
//     public:
//
//     bool FileExist(std::string path);
//
//     static InitInterpolationHelper& Get()
//     {
//         static InitInterpolationHelper instance;
//         return instance;
//     }
//
//     private:
//     InitInterpolationHelper(){};
//     ~InitInterpolationHelper(){};
// };
//
// bool InitInterpolationHelper::FileExist(std::string path)
// {
//     struct stat dummy_stat_return_val;
//
//     if (stat(path.c_str(), &dummy_stat_return_val) != 0)
//     {
//         return false;
//     }
//     else
//     {
//         return true;
//     }
// }

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

// void ParseParametrizations(const std::vector<Parametrization*> parametrizations)
// {
//     if (parametrizations.empty())
//     {
//         log_fatal("At least one parametrization is needed for interpolation initialization!");
//     }
//
//     if (parametrizations.size() == 1)
//     {
//         Parametrization* param = parametrizations.at(0);
//         const ParticleDef& particle_def = param->GetParticleDef();
//         const Medium& medium = param->GetMedium();
//         const EnergyCutSettings& cut_settings = param->GetEnergyCuts();
//
//         std::string param_name;
//
//         BremsPetrukhinShestakov* brems = dynamic_cast<BremsPetrukhinShestakov*>(param);
//         if (brems)
//         {
//             param_name = "B_PS";
//         }
//         BremsKelnerKokoulinPetrukhin* brems;
//         BremsCompleteScreening* brems;
//         BremsAndreevBezrukovBugaev* brems;
//     }
//
//     const ParticleDef& particle_def = parametrizations.at(0)->GetParticleDef();
//     const Medium& medium = parametrizations.at(0)->GetMedium();
//     const EnergyCutSettings& cut_settings = parametrizations.at(0)->GetEnergyCuts();
//
//     std::string param_name;
//     if (dynamic_cast<Ionization*>(parametrizations.at(0)))
//     {
//         param_name = "ioniz";
//     }
//     else if (dynamic_cast<Bremsstrahlung*>(parametrizations.at(0)))
//     {
//         param_name = "brems";
//     }
//     else if (dynamic_cast<EpairProduction*>(parametrizations.at(0)))
//     {
//         param_name = "epair";
//     }
//     else if (dynamic_cast<Photonuclear*>(parametrizations.at(0)))
//     {
//         param_name = "photo";
//     }
//     // std::string param_name =
//     //     parametrizations.size() == 1 ? std::string(typeid(parametrizations[0]).name()) : "";
//
//     // clang-format off
//     filename << pathname
//              << "/"
//              << name
//              << "_"
//              << param_name
//              << "_particle_" << particle_def.name
//              << "_mass_" << particle_def.mass
//              << "_charge_" << particle_def.charge
//              << "_lifetime_" << particle_def.lifetime
//              //  <<"_para_"<<parametrization_
//              << "_med_" << medium.GetName()
//              << "_" << medium.GetMassDensity()
//              << "_ecut_" << cut_settings.GetEcut()
//              << "_vcut_" << cut_settings.GetVcut();
//     // clang-format on
//
//     if (parametrizations.size() == 1)
//     {
//         // clang-format off
//         filename << "_param_" << typeid(parametrizations[0]).name()
//                  << "_lpm_"<<parametrizations[0]->GetDefinition().lpm_effect_enabled
//                  <<"_multiplier_"<< parametrizations[0]->GetDefinition().multiplier;
//         // clang-format on
//     }
//     else
//     {
//         for(std::vector<Parametrization*>::const_iterator it = parametrizations.begin(); it != parametrizations.end(); ++it)
//         {
//             // clang-format off
//             std::string param_name = std::string(typeid(*it).name()).substr(0, 1);
//             filename << "_" << param_name
//                      << "_lpm_" << (*it)->GetDefinition().lpm_effect_enabled
//                      << "_multiplier_" << (*it)->GetDefinition().multiplier;
//             // clang-format on
//
//         }
//     }
//
//     if(!raw)
//     {
//         filename << ".txt";
//     }
// }
// }

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

    std::string pathname = interpolation_def.path_to_tables;
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
        // dndx_interpolant_1d_.resize( medium.GetNumComponents() );
        // dndx_interpolant_2d_.resize( medium.GetNumComponents() );

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
                    // for (std::vector<Interpolant*>::iterator interpol_it = builder_it->second.begin(); interpol_it != builder_it->second.end(); ++interpol_it)
                    // {
                    //     *interpol_it = builder_it->first->build();
                    //     (*interpol_it)->Save(output, raw);
                    // }
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
        for (InterpolantBuilderContainer::iterator builder_it = builder_container.begin(); builder_it != builder_container.end(); ++builder_it)
        {
            (*builder_it->second) = builder_it->first->build();

            // for (std::vector<Interpolant*>::iterator interpol_it = builder_it->second.begin(); interpol_it != builder_it->second.end(); ++interpol_it)
            // {
            //     *interpol_it = builder_it->first->build();
            // }
        }
    }

    log_info("Initialize %s interpolation done.", name.c_str());
}
// ------------------------------------------------------------------------- //
// void InitializeInterpolation(const std::string name,
//                              InterpolantBuilderContainer& builder_container,
//                              size_t hash_digest)
// {
//     using namespace std;
//
//     log_info("Initialize %s interpolation.", name.c_str());
//
//     bool storing_failed =   false;
//     bool reading_worked =   true;
//
//     std::string pathname = "";
//     bool raw = true;
//
//     stringstream filename;
//
//     // --------------------------------------------------------------------- //
//     // Read in data from parametrization and create the filename
//     // --------------------------------------------------------------------- //
//
//     try
//     {
//         //TODO(mario): Maybe check if all data are equal Wed 2017/09/06
//         const ParticleDef& particle_def = parametrizations.at(0)->GetParticleDef();
//         const Medium& medium = parametrizations.at(0)->GetMedium();
//         const EnergyCutSettings& cut_settings = parametrizations.at(0)->GetEnergyCuts();
//
//         pathname = parametrizations.at(0)->GetDefinition().path_to_tables;
//         raw = parametrizations.at(0)->GetDefinition().raw;
//
//         //TODO(mario): remove when enums are introduced Wed 2017/09/06
//         std::string param_name;
//         if (dynamic_cast<Ionization*>(parametrizations.at(0)))
//         {
//             param_name = "ioniz";
//         }
//         else if (dynamic_cast<Bremsstrahlung*>(parametrizations.at(0)))
//         {
//             param_name = "brems";
//         }
//         else if (dynamic_cast<EpairProduction*>(parametrizations.at(0)))
//         {
//             param_name = "epair";
//         }
//         else if (dynamic_cast<Photonuclear*>(parametrizations.at(0)))
//         {
//             param_name = "photo";
//         }
//         // std::string param_name =
//         //     parametrizations.size() == 1 ? std::string(typeid(parametrizations[0]).name()) : "";
//
//         // clang-format off
//         filename << pathname
//                  << "/"
//                  << name
//                  << "_"
//                  << param_name
//                  << "_particle_" << particle_def.name
//                  << "_mass_" << particle_def.mass
//                  << "_charge_" << particle_def.charge
//                  << "_lifetime_" << particle_def.lifetime
//                  //  <<"_para_"<<parametrization_
//                  << "_med_" << medium.GetName()
//                  << "_" << medium.GetMassDensity()
//                  << "_ecut_" << cut_settings.GetEcut()
//                  << "_vcut_" << cut_settings.GetVcut();
//         // clang-format on
//
//         if (parametrizations.size() == 1)
//         {
//             // clang-format off
//             filename << "_param_" << typeid(parametrizations[0]).name()
//                      << "_lpm_"<<parametrizations[0]->GetDefinition().lpm_effect_enabled
//                      <<"_multiplier_"<< parametrizations[0]->GetDefinition().multiplier;
//             // clang-format on
//         }
//         else
//         {
//             for(std::vector<Parametrization*>::const_iterator it = parametrizations.begin(); it != parametrizations.end(); ++it)
//             {
//                 // clang-format off
//                 std::string param_name = std::string(typeid(*it).name()).substr(0, 1);
//                 filename << "_" << param_name
//                          << "_lpm_" << (*it)->GetDefinition().lpm_effect_enabled
//                          << "_multiplier_" << (*it)->GetDefinition().multiplier;
//                 // clang-format on
//
//             }
//         }
//
//         if(!raw)
//         {
//             filename << ".txt";
//         }
//     }
//     catch (const std::out_of_range& e)
//     {
//         log_fatal("At least one parametrization is needed for initializing interpolations.");
//     }
//
//     if(!pathname.empty())
//     {
//         // dndx_interpolant_1d_.resize( medium.GetNumComponents() );
//         // dndx_interpolant_2d_.resize( medium.GetNumComponents() );
//
//         if( FileExist(filename.str()) )
//         {
//             log_debug("Bremsstrahlungs parametrisation tables (dNdx) will be read from file:\t%s",filename.str().c_str());
//
//             std::ifstream input;
//
//             if(raw)
//             {
//                 input.open(filename.str().c_str(), ios::binary);
//             }
//             else
//             {
//                 input.open(filename.str().c_str());
//             }
//
//             for (InterpolantBuilderContainer::iterator builder_it = builder_container.begin(); builder_it != builder_container.end(); ++builder_it)
//             {
//                 //TODO(mario): read check Tue 2017/09/05
//                 (*builder_it->second) = new Interpolant();
//                 (*builder_it->second)->Load(input, raw);
//
//                 // for (std::vector<Interpolant*>::iterator interpol_it = builder_it->second.begin(); interpol_it != builder_it->second.end(); ++interpol_it)
//                 // {
//                 // }
//             }
//
//             input.close();
//         }
//
//         if(!FileExist(filename.str()) || !reading_worked )
//         {
//             if(!reading_worked)
//             {
//                 log_info("file %s is corrupted! Write it again!",filename.str().c_str());
//             }
//
//             log_info("Info: Bremsstrahlungs parametrisation tables (dNdx) will be saved to file:\t%s",filename.str().c_str());
//
//             ofstream output;
//
//             if(raw)
//             {
//                 output.open(filename.str().c_str(), ios::binary);
//             }
//             else
//             {
//                 output.open(filename.str().c_str());
//             }
//
//             if(output.good())
//             {
//                 output.precision(16);
//
//                 for (InterpolantBuilderContainer::iterator builder_it = builder_container.begin(); builder_it != builder_container.end(); ++builder_it)
//                 {
//                     (*builder_it->second) = builder_it->first->build();
//                     (*builder_it->second)->Save(output, raw);
//                     // for (std::vector<Interpolant*>::iterator interpol_it = builder_it->second.begin(); interpol_it != builder_it->second.end(); ++interpol_it)
//                     // {
//                     //     *interpol_it = builder_it->first->build();
//                     //     (*interpol_it)->Save(output, raw);
//                     // }
//                 }
//             }
//             else
//             {
//                 storing_failed  =   true;
//                 log_warn("Can not open file %s for writing! Table will not be stored!",filename.str().c_str());
//             }
//             output.close();
//         }
//     }
//
//     if(pathname.empty() || storing_failed)
//     {
//         for (InterpolantBuilderContainer::iterator builder_it = builder_container.begin(); builder_it != builder_container.end(); ++builder_it)
//         {
//             (*builder_it->second) = builder_it->first->build();
//
//             // for (std::vector<Interpolant*>::iterator interpol_it = builder_it->second.begin(); interpol_it != builder_it->second.end(); ++interpol_it)
//             // {
//             //     *interpol_it = builder_it->first->build();
//             // }
//         }
//     }
// }

} /* Helper */


}  // PROPOSAL
