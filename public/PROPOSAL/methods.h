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



bool StartsWith(const std::string& text,const std::string& token);

//----------------------------------------------------------------------------//

bool endsWith(const std::string& text,const std::string& token);

//----------------------------------------------------------------------------//

int roundValue(double val);

//----------------------------------------------------------------------------//

std::deque<std::string>* splitString(std::string args, std::string Delimiters);

//----------------------------------------------------------------------------//

bool equalsIgnoreCase(std::string s1, std::string s2);

//----------------------------------------------------------------------------//

std::string toLowerCase(std::string toConvert);

//----------------------------------------------------------------------------//

std::string replaceAll(std::string toConvert, char oldChar, char newChar);

//----------------------------------------------------------------------------//

double old_RandomDouble();

//----------------------------------------------------------------------------//

bool compare(std::string name1, std::string name2 , bool output);

//----------------------------------------------------------------------------//

std::string nextToken(std::deque<std::string> *Tokens);

//----------------------------------------------------------------------------//

#define SWAP(a, b,T) {T t; t = a; a = b; b = t;}

#endif /* METHODS_H_ */
