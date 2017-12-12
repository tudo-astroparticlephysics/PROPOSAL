
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

#ifndef METHODS_H_
#define METHODS_H_

// #include <string>
// #include <iostream>
// #include <fstream>
// #include <vector>
#include <deque>

#include <boost/math/special_functions/erf.hpp>

//necessary since the boost function does not work for [-1,1] but for (-1,1)
#define FIXPREC 0.9999999999999999
#define erfInv(x)   boost::math::erf_inv(FIXPREC*x)

namespace PROPOSAL
{

bool FileExist(std::string path);

//----------------------------------------------------------------------------//

bool StartsWith(const std::string& text,const std::string& token);

//----------------------------------------------------------------------------//

bool EndsWith(const std::string& text,const std::string& token);

//----------------------------------------------------------------------------//

int RoundValue(double val);

//----------------------------------------------------------------------------//

std::deque<std::string>* SplitString(std::string args, std::string Delimiters);

//----------------------------------------------------------------------------//

std::string ToLowerCase(std::string toConvert);

//----------------------------------------------------------------------------//

std::string ReplaceAll(std::string toConvert, char oldChar, char newChar);

//----------------------------------------------------------------------------//

double Old_RandomDouble();

//----------------------------------------------------------------------------//

std::string NextToken(std::deque<std::string> *Tokens);

//----------------------------------------------------------------------------//

}

#define SWAP(a, b,T) {T t; t = a; a = b; b = t;}

#endif /* METHODS_H_ */
