/*! \file   methods.h
*   \brief  Header file for the methods routines.
*
*   Some methods which were implemented for often used algorithms.
*
*   \date   21.06.2010
*   \author Jan-Hendrik Koehne
*/
#pragma once

// #include <string>
// #include <iostream>
// #include <fstream>
// #include <vector>
#include <deque>

#include <boost/math/special_functions/erf.hpp>
// #include <utility>
#include <vector>

//necessary since the boost function does not work for [-1,1] but for (-1,1)
#define FIXPREC 0.9999999999999999
// #define erfInv(x) myErfInv2(FIXPREC*x)
#define erfInv(x) boost::math::erf_inv(FIXPREC*x)

namespace PROPOSAL
{

// double RationalApproximation(double t);
// double NormalCDFInverse(double p);
float myErfInv2(float x);

//----------------------------------------------------------------------------//

int RoundValue(double val);


// // ----------------------------------------------------------------------------
// /// @brief CTRP to add name to a class
// // ----------------------------------------------------------------------------
// template <class Base, class Derivied>
// class NamedClass: virtual public Base
// {
//     public:
//     const std::string GetName() { return name_; }
//
//     private:
//     static const std::string name_;
// };
//
// // ----------------------------------------------------------------------------
// /// @brief Clone CTRP
// // ----------------------------------------------------------------------------
// template <class Base, class Derivied>
// class Clonable: virtual public Base
// {
//     public:
//     Base* clone() { return Derivied(static_cast<Derivied>(*this)); }
// };



// ----------------------------------------------------------------------------
/// @brief Error class for system path checks
// ----------------------------------------------------------------------------
// class InvalidPath : public std::exception
// {
//     public:
//     explicit InvalidPath(const char* pathname, const char* err);
//     explicit InvalidPath(const std::string& pathname, const char* err);
//
//     virtual ~InvalidPath() throw() {}
//     virtual const char* what() const throw();
//
//     private:
//     std::string pathname_;
//     std::string msg_;
// };


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

// class ParticleDef;
// class Medium;
// class EnergyCutSettings;
class Parametrization;
class Interpolant;
class InterpolantBuilder;

namespace  Helper
{

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

bool FileExist(const std::string path);

// ----------------------------------------------------------------------------
/// @brief Center string
///
/// @param width
/// @param str
/// @param fill
///
/// @return
// ----------------------------------------------------------------------------
std::string Centered(int width, const std::string& str, char fill = '=');

typedef std::vector<std::pair<InterpolantBuilder*, Interpolant** > > InterpolantBuilderContainer;

// ------------------------------------------------------------------------- //
// ----------------------------------------------------------------------------
/// @brief Helper for interpolation initialization
///
/// @param name
/// @param InterpolantBuilderContainer
/// @param std::vector
// ----------------------------------------------------------------------------
void InitializeInterpolation(const std::string name,
                             InterpolantBuilderContainer&,
                             const std::vector<Parametrization*>&,
                             const InterpolationDef);

} /*  Helper */


} /* PROPOSAL */
