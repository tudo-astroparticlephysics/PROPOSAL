
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
#include <map>
#include "PROPOSAL/json.hpp"

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
/// @brief Definition needed to initialize interpolation
// ----------------------------------------------------------------------------
struct InterpolationDef
{
    InterpolationDef()
        : order_of_interpolation(5)
        , path_to_tables(std::string())
        , path_to_tables_readonly(std::string())
        , max_node_energy(1e14) // upper energy bound for Interpolation (MeV)
        , nodes_cross_section(100) // number of interpolation in cross section
        , nodes_continous_randomization(200) // number of interpolation in continuous randomization
        , nodes_propagate(1000) // number of interpolation in propagate
        , do_binary_tables(true)
        , just_use_readonly_path(false)
    {
    }

    InterpolationDef(const nlohmann::json&);

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

// ----------------------------------------------------------------------------
/// @brief Simple map structure where keys and values can be used for indexing
// ----------------------------------------------------------------------------
template <typename K, typename V>
class Bimap
{
    typedef std::map<K, V> left_map_type;
    typedef std::map<V, K> right_map_type;

public:
    // ------------------------------------------------------------------------
    /// @brief Insert a key value pair into the Bimap
    ///
    /// @param key: Key value.
    /// @param value: Value.
    // ------------------------------------------------------------------------
    void insert(K key, V value)
    {
        if (_left.count(key) == 0 && _right.count(value) == 0)
        {
            _left[key] = value;
            _right[value] = key;
        }
    };

    // ------------------------------------------------------------------------
    /// @brief Remove the key value pair with given key from the Bimap
    ///
    /// @param key: Key value to be removed.
    // ------------------------------------------------------------------------
    void remove(K key)
    {
        auto it = _left.find(key);
        if (it != _left.end())
        {
            _left.erase(it);
            _right.erase(it->second);
        }
    }

    // ------------------------------------------------------------------------
    /// @brief Remove the key value pair with the given value from the Bimap
    ///
    /// @param value: Value to be removed acting as a key.
    // ------------------------------------------------------------------------
    void remove(V value)
    {
        auto it = _right.find(value);
        if (it != _right.end())
        {
            _right.erase(it);
            _left.erase(it->second);
        }
    }

    // ------------------------------------------------------------------------
    /// @brief Get a 'view' of the (key, value) map
    // ------------------------------------------------------------------------
    const left_map_type& GetLeft() const
    {
        return _left;
    }

    // ------------------------------------------------------------------------
    /// @brief Get a 'view' of the (value, key) map
    // ------------------------------------------------------------------------
    const right_map_type& GetRight() const
    {
        return _right;
    }

    // ------------------------------------------------------------------------
    /// @brief Clear the whole map, size is zero after that
    // ------------------------------------------------------------------------
    void clear()
    {
        _left.clear();
        _right.clear();
    }

    // ------------------------------------------------------------------------
    /// @brief Get the current size (number of keys) of the map
    // ------------------------------------------------------------------------
    size_t size() const
    {
        return _left.size();
    }

private:
    left_map_type _left;
    right_map_type _right;
};


} // namespace Helper

} // namespace PROPOSAL
