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
#define erfInv(x) boost::math::erf_inv(FIXPREC*x)

namespace PROPOSAL
{

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

class ParticleDef;
class Medium;
class EnergyCutSettings;
class Parametrization;
class Interpolant;
class InterpolantBuilder;

namespace  Helper
{

// typedef std::vector<std::pair<InterpolantBuilder*, std::vector<Interpolant*> > > InterpolantBuilderContainer;
typedef std::vector<std::pair<InterpolantBuilder*, Interpolant** > > InterpolantBuilderContainer;

// ------------------------------------------------------------------------- //
void InitializeInterpolation(const std::string name,
                             InterpolantBuilderContainer&,
                             const std::vector<Parametrization*>&);

} /*  Helper */


} /* PROPOSAL */

#define SWAP(a, b,T) {T t; t = a; a = b; b = t;}
