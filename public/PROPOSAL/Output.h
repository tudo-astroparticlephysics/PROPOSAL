/*! \file   Output.h
*   \brief  Headerfile for the output routines.
*
*   For more details see the class documentation.
*
*   \date   29.09.2010
*   \author Martin Schmitz
*/


#include "PROPOSAL/CrossSections.h"


#ifndef OUTPUT_H_
#define OUTPUT_H_



#include <string>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include "MMCexception.h"
#include "deque"



/**
 * \brief Class containing the Output Methods.
 *
 * This Class will provide just the fundamentel Methods so far.
 * You will need to implement more methods later
 * Look in the Java source file to do so.
 */

class PROPOSALParticle;

class Output
{

private:

    PROPOSALParticle    *p;        // PROPOSALParticle which needs to be printed somehow
    static int          OUTNUM;

    static std::vector<char>    c;        // needed to format the double.
    const static char           z0='0';   // needed to format the double.
    const static char           zm='-';   // needed to format the double.
    const static char           zd='.';   // needed to format the double.
    const static char           ze='e';   // needed to format the double.
    static long                 POWOUT;   // needed to format the double.


public:

//----------------------------------------------------------------------------//
    static std::ofstream    output_P;
    std::stringstream       history;
    static std::fstream     ReadWrite;

    int     igen;
    int     gens;

    static std::string  version;

    const static int    HISTSIZE   =   1000;
    int                 HIST;                   //!< output "styler", preset is -1
    static bool         AMASIM;

    std::streambuf      *stdout_backup;

    std::fstream        err;

    static bool         outf;
    static bool         inf;
    static bool         I3flag;
    static bool         RecDec;

    bool            DEBUG;

    static bool     raw;
    static bool     texi;

    std::vector<PROPOSALParticle *> I3hist;

    std::stringstream out;

 //----------------------------------------------------------------------------//
    // Constructors


    Output(PROPOSALParticle *part);

//----------------------------------------------------------------------------//
    Output();

//----------------------------------------------------------------------------//
    //Memberfunctions

//----------------------------------------------------------------------------//
    /**
     * enables redirection of the stderr output into a file
     */

    static bool setStderr(std::string name) throw (MMCexception);
//----------------------------------------------------------------------------//

    /**
     * enables redirection of the stdout output into a string
     */

    void setStdout() ;

//----------------------------------------------------------------------------//


    static void particleToFileAscii(PROPOSALParticle *p,std::string filename, bool isPrimary);

//----------------------------------------------------------------------------//


    static void particleToFileRaw(PROPOSALParticle *p, std::string filename, bool isPrimary);

//----------------------------------------------------------------------------//

    /**
     * gets the stdout output string, resets the stdout output
     */

    std::string getStdout();// throw (MMCexception);

//----------------------------------------------------------------------------//

    static double read();
//----------------------------------------------------------------------------//

    static std::string readStr();
//----------------------------------------------------------------------------//
    static void write(double g);

//----------------------------------------------------------------------------//
    //Tomasz
    static void write(std::string g);

//----------------------------------------------------------------------------//
    /**
     *format the double
    */

    static std::string f(double d);
//----------------------------------------------------------------------------//

    void output(int wint, std::string comp, double de, double ef);

//----------------------------------------------------------------------------//

    void init(std::string name);

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

//----------------------------------------------------------------------------//
    /**
     *  initialize for regular history output
     */

    void initDefault(std::string type);

//----------------------------------------------------------------------------//
    /**
     *  initialize for regular history output
     */

    void initDefault(int igen, int gens, std::string type, double time, double x, double y, double z, double theta, double phi);

//----------------------------------------------------------------------------//
    /**
     * initialize for F2000 history output
     */

    void initF2000(int igen, int gens, std::string type, double time, double x, double y, double z, double theta, double phi);

//----------------------------------------------------------------------------//
    // Getter

    bool get_DEBUG()
    {
        return DEBUG;
    }

    bool get_AMASIM()
    {
        return AMASIM;
    }
//----------------------------------------------------------------------------//
    // Setter

    static void set_I3flag(bool new_flag)
    {
        I3flag = new_flag;
    }

    void set_DEBUG(bool newDebug)
    {
        DEBUG = newDebug;
    }


};

#endif /* OUTPUT_H_ */
