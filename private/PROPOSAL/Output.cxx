/*! \file   Output.cxx
*   \brief  Source file for the output routines.
*
*   For more details see the class documentation.
*
*   \date   29.09.2010
*   \author Martin Schmitz
*/


#include "PROPOSAL/CrossSections.h"
#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/methods.h"
#include <cmath>
#include <sstream>
#include "PROPOSAL/MMCexception.h"
#include "PROPOSAL/Output.h"
#include <sys/stat.h>
#include <iomanip>

bool Output::I3flag         =   false;
std::string Output::version =   "Lepton Propagation Code in C++: PROPOSAL";
long Output::POWOUT         =   0;
bool Output::AMASIM         =   false;
int Output::OUTNUM          =   8;
bool Output::raw            =   false;
bool Output::texi           =   false;
bool Output::inf            =   false;
bool Output::outf           =   false;
bool Output::RecDec         =   false;

std::fstream Output::ReadWrite;
std::vector<char> Output::c;
std::ofstream Output::output_P;

using namespace std;

// Constructors

Output::Output(PROPOSALParticle *part)
{

    version="Lepton Propagation Code in C++: PROPOSAL";

    AMASIM  =   false;
    this->p =   part;
    DEBUG   =   false;
    HIST    =   -1;
    OUTNUM  =   8;

    c.resize(16);

    POWOUT  =   (long)roundValue(pow((double) 10, 8));
    texi    =   false;

    I3hist.resize(HISTSIZE);

    inf     =   false;
    outf    =   false;


}

//----------------------------------------------------------------------------//

bool Output::setStderr(std::string name) throw(MMCexception)
{

    try
    {
        if(freopen(name.c_str(),"a",stderr)==0)
        {
            throw MMCexception(0,"output", name);
        }

        return true;
    }
    catch(MMCexception &error)
    {
        cout<<error.what()<<endl;
        return false;
    }

}

//----------------------------------------------------------------------------//

void Output::setStdout()
{
    out.str("");
    out.clear();
    stdout_backup   =   cout.rdbuf();
    cout.rdbuf( out.rdbuf() );

}

//----------------------------------------------------------------------------//

string Output::getStdout()
{
    cout.flush();
    cout.rdbuf( stdout_backup );
    return out.str();

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

string Output::readStr()
{

    string data;

    if(raw)
    {
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
            cerr<< "***WARNING***: readStr failed: "<<endl;
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

void Output::write(string g)
{
    if(ReadWrite.good())
    {
        if(raw)
        {
            //ReadWrite.write(reinterpret_cast <char*> (&g), sizeof( g ));
            //muss implementiert werden
        }
        else
        {
            ReadWrite<<g<<endl;
        }
    }
    else
    {
        cerr<<"***WARNING***: writeStr failed: "<<endl;
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

void Output::output(int wint, string comp, double de, double ef)
{
    string str;

    if(DEBUG || (HIST==1))
    {
        str =   " . ";

        if(wint==0)
        {
            if(comp=="conti")
            {
                str +=  "lost continuously " + f(de) + " MeV over " + f(ef) + " cm up to " + f(p->r) + " cm\n";
            }
            else
            {
                str +=  "muon with energy "+f(de)+" MeV traveled  "+f(ef)+" cm\n";
            }
        }
        else
        {
            if(wint==1)
            {
                str +=  "decayed into " + comp;
            }
            else if(wint==2)
            {
                str +=  "ionized by " + comp;
            }
            else if(wint==3)
            {
                str +=  "bremsed by " + comp;
            }
            else if(wint==4)
            {
                str +=  "photoed by " + comp;
            }
            else if(wint==5)
            {
                str +=  "epaired by " + comp;
            }

            str +=  " ... lost " + f(de) + " MeV, ef = " + f(ef) + " MeV, xf = " + f(p->r) + " cm\n";
        }

        if(HIST==1)
        {
            history<<str<<endl;
        }

        if(DEBUG)
        {
            cerr<<str<<endl;
        }
    }

    if(HIST==2 || (RecDec && wint==1))
    {
        if(wint==0 || wint==1)
        {
            str =   comp;
        }
        else if(wint==2)
        {
            str =   "delta";
        }
        else if(wint==3)
        {
            str =   "brems";
        }
        else if(wint==4)
        {
            str =   "munu";
        }
        else if(wint==5)
        {
            str =   "epair";
        }
        else
        {
            str =   "";
        }

        if(str!="")
        {
            gens++;

            if(I3flag)
            {
                PROPOSALParticle *PROPOSALParticleToPrint   =   new PROPOSALParticle(igen, gens, str, p->x, p->y, p->z, p->theta, p->phi, de, p->t, wint==0?ef:0);
                I3hist.push_back(PROPOSALParticleToPrint);

            }
            else
            {
                std::ostringstream Print_igen;
                std::ostringstream Print_gens;
                Print_gens <<gens;
                Print_igen <<igen;
                str =   "TR " + Print_gens.str() + " " + Print_igen.str() + " " + str + " " + f(p->x*1.e-2) + " " + f(p->y*1.e-2) + " " + f(p->z*1.e-2);
                str+=" " + f(180-p->theta) + " " + f(p->phi<180?p->phi+180:p->phi-180) + " " + f(wint==0?ef*1.e-2:0) + " " + f(de*1.e-3) + " " + f(p->t*1.e9) + "\n";
                history<<str<<endl;
            }
        }
    }
}


// split the string into substrings using given delimeter list
deque<string> Output::splitString(string args, string Delimiters)
{

    deque<string> Tokens;
    string::size_type Next, EndNext=0;

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
        Tokens.push_back(args.substr(Next, EndNext - Next));
    }

    return Tokens;
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

void Output::init(string name)
{

    HIST    =   0;
    p->location(name, 0, 0, 0, 0, 0, 0);

}

//----------------------------------------------------------------------------//

void Output::initDefault(string type)
{
    HIST    =   1;
    p->location(type, 0, 0, 0, 0, 0, 0);
}

//----------------------------------------------------------------------------//

void Output::initDefault(int igen, int gens, string type, double time, double x, double y, double z, double theta, double phi)
{
    HIST        =   3;
    this->igen  =   igen;
    this->gens  =   gens;
    p->location(type, time, x, y, z, theta, phi);

    if(RecDec)
    {
        if(I3flag)
        {
            I3hist.clear();
        }
    }

}

//----------------------------------------------------------------------------//

void Output::initF2000(int igen, int gens, string type, double time, double x, double y, double z, double theta, double phi)
{
    HIST        =   2;
    this->igen  =   igen;
    this->gens  =   gens;
    p->location(type, time, x, y, z, theta, phi);

    if(I3flag)
    {
        I3hist.clear();
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

//----------------------------------------------------------------------------//

void Output::particleToFileAscii(PROPOSALParticle *p, string filename, bool isPrimary)
{


    output_P.precision(15);

    if(FileExist(filename))
    {
        output_P.open(filename.c_str(),ios_base::app);
    }
    else
    {
        output_P.open(filename.c_str());
        //output<<"# The following line is the branchDiscriptor which is needed by the TTree::ReadFile method"<<endl;
        //output<<"x/D:y:z:t:l:theta:phi:e:type/I:xi/D:yi:zi:Ei:xc:yc:zc:Ec:xf:yf:zf:Ef:Elost:name/C"<<endl;
    }

    if(isPrimary)
    {

        output_P<<"#\t";
        output_P<<p->r<<"\t";
        output_P<<p->x<<"\t";
        output_P<<p->y<<"\t";
        output_P<<p->z<<"\t";
        output_P<<p->t<<"\t";
        output_P<<p->l<<"\t";
        output_P<<p->theta<<"\t";
        output_P<<p->phi<<"\t";
        output_P<<p->e<<"\t";
        output_P<<p->type<<"\t";
        output_P<<p->xi<<"\t";
        output_P<<p->yi<<"\t";
        output_P<<p->zi<<"\t";
        output_P<<p->Ei<<"\t";
        output_P<<p->xc<<"\t";
        output_P<<p->yc<<"\t";
        output_P<<p->zc<<"\t";
        output_P<<p->Ec<<"\t";
        output_P<<p->xf<<"\t";
        output_P<<p->yf<<"\t";
        output_P<<p->zf<<"\t";
        output_P<<p->Ef<<"\t";
        output_P<<p->Elost<<"\t";
        output_P<<p->name<<"\n";
    }
    else
    {
        output_P<<p->x<<"\t";
        output_P<<p->y<<"\t";
        output_P<<p->z<<"\t";
        output_P<<p->t<<"\t";
        output_P<<p->l<<"\t";
        output_P<<p->theta<<"\t";
        output_P<<p->phi<<"\t";
        output_P<<p->e<<"\t";
        output_P<<p->type<<"\t";
        output_P<<p->name<<"\n";
    }

    output_P.close();
}



void Output::particleToFileRaw(PROPOSALParticle *p, string filename, bool isPrimary)
{
    char c  =   '#'; // primary identifier
    char d  =   '&'; // secondary identifier
    char n[10];

    for(size_t i =0; i<p->name.length(); i++)
    {
        n[i]    =   p->name.at(i);
    }

    for(size_t j=p->name.length(); j<10; j++)
    {
        n[j]    =   '\0';
    }

    if(FileExist(filename))
    {
        output_P.open(filename.c_str(),ios_base::app|ios::binary);
    }
    else
    {
        output_P.open(filename.c_str(),ios::binary);
    }

    if(isPrimary)
    {


        output_P.write(&c, sizeof( c ));
        output_P.write(reinterpret_cast <char *> (&p->r), sizeof( p->r ));
        output_P.write(reinterpret_cast <char *> (&p->x), sizeof( p->x ));
        output_P.write(reinterpret_cast <char *> (&p->y), sizeof( p->y ));
        output_P.write(reinterpret_cast <char *> (&p->z), sizeof( p->z ));
        output_P.write(reinterpret_cast <char *> (&p->t), sizeof( p->t ));
        output_P.write(reinterpret_cast <char *> (&p->l), sizeof( p->l ));
        output_P.write(reinterpret_cast <char *> (&p->theta), sizeof( p->theta ));
        output_P.write(reinterpret_cast <char *> (&p->phi), sizeof( p->phi ));
        output_P.write(reinterpret_cast <char *> (&p->e), sizeof( p->e ));
        output_P.write(reinterpret_cast <char *> (&p->type), sizeof( p->type ));
        output_P.write(reinterpret_cast <char *> (&p->xi), sizeof( p->xi ));
        output_P.write(reinterpret_cast <char *> (&p->yi), sizeof( p->yi ));
        output_P.write(reinterpret_cast <char *> (&p->zi), sizeof( p->zi ));
        output_P.write(reinterpret_cast <char *> (&p->Ei), sizeof( p->Ei ));
        output_P.write(reinterpret_cast <char *> (&p->xc), sizeof( p->xc ));
        output_P.write(reinterpret_cast <char *> (&p->yc), sizeof( p->yc ));
        output_P.write(reinterpret_cast <char *> (&p->zc), sizeof( p->zc ));
        output_P.write(reinterpret_cast <char *> (&p->Ec), sizeof( p->Ec ));
        output_P.write(reinterpret_cast <char *> (&p->xf), sizeof( p->xf ));
        output_P.write(reinterpret_cast <char *> (&p->yf), sizeof( p->yf ));
        output_P.write(reinterpret_cast <char *> (&p->zf), sizeof( p->zf ));
        output_P.write(reinterpret_cast <char *> (&p->Ef), sizeof( p->Ef ));
        output_P.write(reinterpret_cast <char *> (&p->Elost), sizeof( p->Elost ));
        output_P.write(reinterpret_cast <char *> (&n), sizeof( n ));

    }
    else
    {
        output_P.write(&d, sizeof( d ));
        output_P.write(reinterpret_cast <char *> (&p->x), sizeof( p->x ));
        output_P.write(reinterpret_cast <char *> (&p->y), sizeof( p->y ));
        output_P.write(reinterpret_cast <char *> (&p->z), sizeof( p->z ));
        output_P.write(reinterpret_cast <char *> (&p->l), sizeof( p->l ));
        output_P.write(reinterpret_cast <char *> (&p->t), sizeof( p->t ));
        output_P.write(reinterpret_cast <char *> (&p->theta), sizeof( p->theta ));
        output_P.write(reinterpret_cast <char *> (&p->phi), sizeof( p->phi ));
        output_P.write(reinterpret_cast <char *> (&p->e), sizeof( p->e ));
        output_P.write(reinterpret_cast <char *> (&p->type), sizeof( p->type ));
        output_P.write(reinterpret_cast <char *> (&n), sizeof( n ));


    }

    output_P.close();
}
