/*! \file   methods.h
*   \brief  Header file for the methods routines.
*
*   Some methods which were implemented for often used algorithms.
*
*   \date   21.06.2010
*   \author Jan-Hendrik Koehne
*/

#ifndef METHODS_H_
#define METHODS_H_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "deque"

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

bool EqualsIgnoreCase(std::string s1, std::string s2);

//----------------------------------------------------------------------------//

std::string ToLowerCase(std::string toConvert);

//----------------------------------------------------------------------------//

std::string ReplaceAll(std::string toConvert, char oldChar, char newChar);

//----------------------------------------------------------------------------//

double Old_RandomDouble();

//----------------------------------------------------------------------------//

bool Compare(std::string name1, std::string name2 , bool output);

//----------------------------------------------------------------------------//

std::string NextToken(std::deque<std::string> *Tokens);

//----------------------------------------------------------------------------//

#define SWAP(a, b,T) {T t; t = a; a = b; b = t;}

#endif /* METHODS_H_ */
