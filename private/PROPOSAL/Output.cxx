/*! \file   Output.cxx
*   \brief  Source file for the output routines.
*
*   For more details see the class documentation.
*
*   \date   2013.03.12
*   \author Jan-Hendrik Köhne
*/


#include "PROPOSAL/methods.h"
#include <cmath>
#include <sstream>
#include "PROPOSAL/Output.h"
#include <sys/stat.h>
#include <iomanip>

bool Output::inf            =   false;
bool Output::outf           =   false;
bool Output::raw            =   false;
int Output::OUTNUM          =   13;
std::fstream Output::ReadWrite;

using namespace std;

// Constructors

Output::Output()
{
    inf     =   false;
    outf    =   false;
}


//----------------------------------------------------------------------------//


double Output::read()
{

    double data;

    if(raw)
    {
        if(ReadWrite.good())
        {
            ReadWrite.read( reinterpret_cast<char *>( &data ), sizeof data );
        }
        else
        {
            throw 0;
            cerr<< "***WARNING***: read failed: "<<endl;
            cerr<< "parametrization table read error"<<endl;
        }

        return data;
    }
    else
    {
        if(ReadWrite.good())
        {
            ReadWrite>>data;
        }
        else
        {
            throw 0;
            cerr<< "***WARNING***: read failed: "<<endl;
            cerr<< "parametrization table read error"<<endl;
        }

        return data;
    }
}


//----------------------------------------------------------------------------//

void Output::write(double g)
{
    if(ReadWrite.good())
    {
        if(raw)
        {
            ReadWrite.write(reinterpret_cast <char *> (&g), sizeof( g ));
        }
        else
        {
            ReadWrite<<g<<endl;
        }
    }
    else
    {
        cerr<<"***WARNING***: write failed: "<<endl;
    }

}


//----------------------------------------------------------------------------//

void Output::open(std::string name)
{
    inf     =   false;
    outf    =   false;
    ReadWrite.precision(16);

    if(FileExist(name))
    {
        if(raw)
        {
            ReadWrite.open(name.c_str(), ios::in|ios::binary);
        }
        else
        {
            ReadWrite.open(name.c_str(), ios::in);
        }

        inf =   true;
        cerr<<"Parametrization tables will be read in from the file "<<name<<endl;
    }
    else
    {
        if(raw)
        {
            ReadWrite.open(name.c_str(), ios::out |ios::binary);
        }
        else
        {
            ReadWrite.open(name.c_str(), ios::out);
        }

        outf    =   true;
        cerr<<"Parametrization tables will be saved to the file "<<name<<endl;
    }


}

//----------------------------------------------------------------------------//


string Output::f(double d)
{

    stringstream ss;
    ss.precision(OUTNUM);
    ss<<scientific<<d;
    return ss.str();

}

//----------------------------------------------------------------------------//

bool Output::FileExist(string path)
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

void Output::close()
{
    inf     =   false;
    outf    =   false;

    ReadWrite.close();
}

//----------------------------------------------------------------------------//

void Output::Delete(string name)
{

    Output::close();
    cout<<"... Deleting corrupt file "<<endl;

    if( remove( name.c_str() ) != 0 )
    {
        cout<< "Error deleting file" <<endl;
    }
    else
    {
        cout<< "File successfully deleted" <<endl;
    }
}
