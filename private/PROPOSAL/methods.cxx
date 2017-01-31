/*! \file   methods.cxx
*   \brief  Source file for the methods routines.
*
*   For more details see the class documentation.
*
*   \date   21.06.2010
*   \author Jan-Hendrik Koehne
*/


#include <string>
#include "PROPOSAL/methods.h"
#include <cmath>
#include <stdlib.h>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <sys/stat.h>


using namespace std;



bool FileExist(string path)
{
    struct stat dummy_stat_return_val;

    if (stat(path.c_str(), &dummy_stat_return_val) != 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// returns true if text starts with token, false otherwise

bool StartsWith(const std::string& text,const std::string& token)
{

    if(text.length() < token.length())
    {
        return false;
    }
	return (text.compare(0, token.length(), token) == 0);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool EndsWith(const std::string& text,const std::string& token)
{
	return (text.compare(text.length()-1,token.length(),token) ==0);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// round a given double to the closest int

int RoundValue(double val)
{
    bool minus      =   false;
    int valRound    =   0;

    if(val<0)
    {
        val     *=  -1;
        minus   =   true;
    }

    val         +=  0.5;
    valRound    =   (int)val;

    if(minus)
    {
        valRound    *=  -1;
    }

    return valRound;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// Compares strings case-insensitive

bool EqualsIgnoreCase(std::string s1, std::string s2)
{
    if(s1.length()==s2.length())
    {
        for(int i=0; i<(int)s1.length(); i++)
        {
            char buffer;

            if(s1.at(i)==s2.at(i))
            {
                continue;
            }
            else
            {
                buffer  =   (char)toupper(s1.at(i));

                if(s2.at(i)==buffer)
                {
                    continue;
                }
                else
                {
                    buffer  =   (char)tolower(s1.at(i));

                    if(s2.at(i)==buffer)
                    {
                        continue;
                    }
                    else
                    {
                        return false;
                    }
                }
            }
        }

        return true;
    }

    else
    {
        return false;
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



string ToLowerCase(string toConvert)
{
    string buffer;

    for(int i=0; i<(int)toConvert.length(); i++)
    {
        buffer  +=  (char)tolower(toConvert.at(i));
    }


    return buffer;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


string ReplaceAll(string toConvert, const char oldChar, const char newChar)
{
    string buffer="";

    for(int i =0; i<(int)toConvert.size(); i++)
    {
        if(toConvert.at(i)==oldChar)
        {
            buffer  +=  newChar;
        }
        else
        {
            buffer  +=  toConvert.at(i);
        }
    }

    return buffer;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// returns a random Double

double Old_RandomDouble()
{
	double result;

    result  =   rand() + 0.0;
    result  =   result / RAND_MAX;

	return result;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


deque<string>* SplitString(string args, string Delimiters)
{

    deque<string> *Tokens   =   new deque<string>;
    string::size_type Next, EndNext =   0;

    while (EndNext != string::npos)
    {
        // Find next token
        Next    =   args.find_first_not_of(Delimiters, EndNext);
        EndNext =   args.find_first_of(Delimiters, Next);

        // Stop if end of string reached
        if (Next == string::npos)
        {
            break;
        }

        // Add token to vector.
        Tokens->push_back(args.substr(Next, EndNext - Next));
    }

    return Tokens;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


string NextToken(deque<string> *Tokens)
{
    string next;
    next    =   Tokens->front();
    Tokens->pop_front();

    return next;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Compare(string name1, string name2 , bool output)
{
    ifstream file;
    ifstream file2;
    double zeile;
    double zeile2;

    int linenumber  =   1;
    bool fail       =   false;

    // open file
    file.open(name1.c_str(), ios::in);
    file2.open (name2.c_str(), ios::in);

    while (file.good() && file2.good())
    {
        // read files line by line

        file >> zeile;
        file2 >> zeile2;

        if(!file.eof() && !file2.eof())
        {
            if( zeile - zeile2 != 0)
            {
                if(output)
                cout << "Problem in line " << linenumber << "\t Value in: "<<name1<<":" << zeile << "\t Value in "<<name2<<": " << zeile2 << endl;
                fail=true;
            }
            ++linenumber;
        }
    } // end of while

    if(fail==true)
    {
        return false;
    }
    else
    {
        return true;
    }
}

