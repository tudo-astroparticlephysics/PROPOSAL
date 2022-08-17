
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

#include <string>
#include <functional>
#include <map>
#include <memory>
#include <vector>
#include <sys/stat.h>
// use unistd.h for access on POSIX os, use io.h for windows systems
// todo: use filesystem instead when switching to c++17
#ifdef _WIN32
#include <io.h>
   #define access    _access_s
#else
    #include <unistd.h>
#endif
#define PROPOSAL_MAKE_HASHABLE(type, ...)                                      \
    namespace std {                                                            \
    template <> struct hash<type> {                                            \
        std::size_t operator()(const type& t) const                            \
        {                                                                      \
            std::size_t ret = 0;                                               \
            PROPOSAL::hash_combine(ret, __VA_ARGS__);                          \
            return ret;                                                        \
        }                                                                      \
    };                                                                         \
    }

namespace PROPOSAL {

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

inline void hash_combine(std::size_t& seed) { (void)seed; }

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
inline void hash_combine(std::size_t& seed, const T& v, Rest... rest)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    hash_combine(seed, rest...);
}

class Interpolant;
struct InterpolantBuilder;

// This object checks whether the table under path/filename is already existing. Needs to be created before the table
// creation process has started.
struct LogTableCreation {
    LogTableCreation(const std::string& path, const std::string& filename);

private:
    static std::string warn_for_path;
};

namespace Helper {

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

    // Comparison function to be used for case-insensitive strings in a map
    struct case_insensitive_comp {
        bool operator()(std::string const& lhs, std::string const& rhs) const;
    };

    // TODO: use std::filesystem when we switch to c++17
    inline bool file_exists (const std::string& path_to_file) {
        return access( path_to_file.c_str(), 0 ) == 0;
    }

    inline bool is_folder_writable(const std::string& path_to_file) {
        return access( path_to_file.c_str(), 2 ) == 0;
    }

} // namespace Helper


} // namespace PROPOSAL
