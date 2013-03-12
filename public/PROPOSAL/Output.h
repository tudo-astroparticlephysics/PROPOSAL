/*! \file   Output.h
*   \brief  Headerfile for the output routines.
*
*   For more details see the class documentation.
*
*   \date   2013.03.12
*   \author Jan-Hednrik Köhne
*/

#ifndef OUTPUT_H_
#define OUTPUT_H_


#include <string>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include "deque"



/**
 * \brief Class containing the Output Methods.
 *
 * This Class will provide just the fundamentel Methods so far.
 * You will need to implement more methods later
 * Look in the Java source file to do so.
 */


class Output
{

private:


public:

//----------------------------------------------------------------------------//
    static std::fstream     ReadWrite;

    static bool         outf;
    static bool         inf;
    static bool         raw;
    static bool         texi;
    static int          OUTNUM;
 //----------------------------------------------------------------------------//
    // Constructors

    Output();

//----------------------------------------------------------------------------//
    //Memberfunctions


    static double read();
//----------------------------------------------------------------------------//

    static void write(double g);

//----------------------------------------------------------------------------//

    /**
     *format the double
    */

    static std::string f(double d);
//----------------------------------------------------------------------------//


    static void open(std::string name);

//----------------------------------------------------------------------------//

    static void close();

//----------------------------------------------------------------------------//

    static bool FileExist(std::string path);

//----------------------------------------------------------------------------//
    static void Delete(std::string name);


//----------------------------------------------------------------------------//

    static std::deque<std::string> splitString(std::string args, std::string Delimiters);


};

#endif /* OUTPUT_H_ */
