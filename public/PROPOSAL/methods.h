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
#define erfInv(x) myErfInv2(x)
// #define erfInv(x) boost::math::erf_inv(FIXPREC*x)

namespace PROPOSAL
{

// double RationalApproximation(double t);
// double NormalCDFInverse(double p);
float myErfInv2(float x);

bool FileExist(std::string path);

//----------------------------------------------------------------------------//

// bool StartsWith(const std::string& text,const std::string& token);

//----------------------------------------------------------------------------//

// bool EndsWith(const std::string& text,const std::string& token);

//----------------------------------------------------------------------------//

int RoundValue(double val);

//----------------------------------------------------------------------------//

// std::deque<std::string>* SplitString(std::string args, std::string Delimiters);

//----------------------------------------------------------------------------//

// std::string ToLowerCase(std::string toConvert);

//----------------------------------------------------------------------------//

// std::string ReplaceAll(std::string toConvert, char oldChar, char newChar);

//----------------------------------------------------------------------------//

// double Old_RandomDouble();

//----------------------------------------------------------------------------//

// std::string NextToken(std::deque<std::string> *Tokens);

//----------------------------------------------------------------------------//


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
                             InterpolationDef);

} /*  Helper */


} /* PROPOSAL */
