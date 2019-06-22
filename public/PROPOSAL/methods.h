
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

#include <vector>
#include <functional>

#define PROPOSAL_MAKE_HASHABLE(type, ...) \
    namespace std {\
        template<> struct hash<type> {\
            std::size_t operator()(const type &t) const {\
                std::size_t ret = 0;\
                PROPOSAL::hash_combine(ret, __VA_ARGS__);\
                return ret;\
            }\
        };\
    }

namespace PROPOSAL {

inline void hash_combine(std::size_t& seed) { (void) seed; }

// ----------------------------------------------------------------------------
/// @brief Function to combine hash values
///
/// This is the implementation of boost::hash_combine as
/// variadic template to combine multiple values at once
///
/// @param seed
/// @param v
/// @param rest
// ----------------------------------------------------------------------------
template <typename T, typename... Rest>
inline void hash_combine(std::size_t& seed, const T& v, Rest... rest) {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    hash_combine(seed, rest...);
}

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
double inverseErrorFunction(double x);

// ----------------------------------------------------------------------------
/// @brief Calculate the dilogarithm
///
/// @param x real argument
///
/// Originally translated by R.Brun from CERNLIB DILOG function C332
///
/// Implemented as a truncated series expansion in terms of Chebyshev
/// polynomials, see [Yudell L. Luke: Mathematical functions and their
/// approximations, Academic Press Inc., New York 1975, p.67].
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
        , path_to_tables(std::string())
        , path_to_tables_readonly(std::string())
        , max_node_energy(1e14)
        , nodes_cross_section(100)
        , nodes_continous_randomization(200)
        , nodes_propagate(1000)
        , do_binary_tables(true)
        , just_use_readonly_path(false)
    {
    }

    int order_of_interpolation;
    std::string path_to_tables;
    std::string path_to_tables_readonly;
    double max_node_energy;
    int nodes_cross_section;
    int nodes_continous_randomization;
    int nodes_propagate;
    bool do_binary_tables;
    bool just_use_readonly_path;

    size_t GetHash() const;
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
/// @brief Check if a given directory has read permissions
///
/// @param table_dir directory to check read permissions
///
/// @return bool
// ----------------------------------------------------------------------------
bool IsReadable(std::string table_dir);

// ----------------------------------------------------------------------------
/// @brief Resolve given path
//
/// Environment variables are tried to expand and relative path will be
/// converted to absolute paths.
///
/// @param std::string path
/// @param bool only check if the path is readable, not writable
///
/// @return resolved path or empty path if errors occured.
// ----------------------------------------------------------------------------
std::string ResolvePath(const std::string&, bool checkReadonly=false);

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
