/*! \file   methods.h
 *   \brief  Header file for the methods routines.
 *
 *   Some methods which were implemented for often used algorithms.
 *
 *   \date   21.06.2010
 *   \author Jan-Hendrik Koehne
 */
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
