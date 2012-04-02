/*
 * methods.h
 *
 *  Created on: 21.06.2010
 *      Author: koehne
 */

#ifndef METHODS_H_
#define METHODS_H_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "deque"



bool StartsWith(const std::string& text,const std::string& token);

bool endsWith(const std::string& text,const std::string& token);

int roundValue(double val);

std::deque<std::string>* splitString(std::string args, std::string Delimiters);

bool equalsIgnoreCase(std::string s1, std::string s2);

std::string toLowerCase(std::string toConvert);

std::string replaceAll(std::string toConvert, char oldChar, char newChar);

double old_RandomDouble();

bool compare(std::string name1, std::string name2 , bool output);

std::string nextToken(std::deque<std::string> *Tokens);

//double randomdouble();

#endif /* METHODS_H_ */
