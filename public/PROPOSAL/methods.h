
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

#include <deque>

#include <boost/math/special_functions/erf.hpp>
#include <vector>

// necessary since the boost function does not work for [-1,1] but for (-1,1)
#define FIXPREC 0.9999999999999999
// #define erfInv(x) myErfInv2(FIXPREC*x)
#define erfInv(x) boost::math::erf_inv(FIXPREC* x)

namespace PROPOSAL {

// ----------------------------------------------------------------------------
/// @brief Implementation of inverse error function
///
/// It is done for performance test with valgrind --tool=callgrind.
/// This tool has problems with boost::math::erf_inv.
///
/// @param x
///
/// @return inverse error function of at x
// ----------------------------------------------------------------------------
float myErfInv2(float x);

//----------------------------------------------------------------------------//
int RoundValue(double val);

// ----------------------------------------------------------------------------
/// @brief Calculate the dilogarithums
///
/// @param x
///
/// @return dilog(x)
// ----------------------------------------------------------------------------
double dilog(double x);

// ----------------------------------------------------------------------------
/// @brief Definition needed to initialize interpolation
// ----------------------------------------------------------------------------
struct InterpolationDef
{
    InterpolationDef()
        : order_of_interpolation(5)
        , path_to_tables("")
        , raw(true)
    {
    }

    int order_of_interpolation;
    std::string path_to_tables;
    bool raw;
};

class Parametrization;
class Interpolant;
class InterpolantBuilder;

namespace Helper {


// ----------------------------------------------------------------------------
/// @brief Check if a given directory has write permissions
///
/// @param table_dir directory to check write permissions
///
/// @return bool
// ----------------------------------------------------------------------------
bool IsWritable(std::string table_dir);

// ----------------------------------------------------------------------------
/// @brief Resolve given path
//
/// Environment variables are tried to expand and relative path will be
/// converted to absolute paths.
///
/// @param std::string path
///
/// @return resolved path or empty path if errors occured.
// ----------------------------------------------------------------------------
std::string ResolvePath(const std::string&);

// ----------------------------------------------------------------------------
/// @brief Uses stat() from sys/stat.h to determine if a file exists
///        in the file system.
///
/// @param path: path to file
///
/// @return bool
// ----------------------------------------------------------------------------
bool FileExist(const std::string path);

// ----------------------------------------------------------------------------
/// @brief Center string
///
/// @param width: length of the resulting string
/// @param str: the body of the centered string
/// @param fill: fill the white space with this character
///
/// @return the centered formatted string
// ----------------------------------------------------------------------------
std::string Centered(int width, const std::string& str, char fill = '=');

typedef std::vector<std::pair<InterpolantBuilder*, Interpolant**> > InterpolantBuilderContainer;

// ----------------------------------------------------------------------------
/// @brief Helper for interpolation initialization
///
/// @param name: subject of resulting file name
/// @param InterpolantBuilderContainer:
///        vector of builder, pointer to Interplant pairs
/// @param std::vector: vector of parametrizations used to create
///        the interpolation tables with
// ----------------------------------------------------------------------------
void InitializeInterpolation(const std::string name,
                             InterpolantBuilderContainer&,
                             const std::vector<Parametrization*>&,
                             const InterpolationDef);

} // namespace Helper

} // namespace PROPOSAL
