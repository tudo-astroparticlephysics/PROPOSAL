/*
 * Output.cxx
 *
 *  Created on: 29.09.2010
 *      Author: schmitz
 */

#include "CrossSections.h"
#include "PROPOSALParticle.h"
#include "methods.h"
#include <cmath>
#include <sstream>
#include "MMCexception.h"
#include "Output.h"
#include <sys/stat.h>
#include <iomanip>

bool Output::I3flag=false;
std::string Output::version="Muon Propagation Code in Java v. 1.5.5";
long Output::POWOUT=0;
bool Output::AMASIM=false;
int Output::OUTNUM=8;
std::vector<char> Output::c;
bool Output::raw=false;
bool Output::texi=false;
std::fstream Output::ReadWrite;
bool Output::inf=false;
bool Output::outf = false;
bool Output::RecDec=false;
std::ofstream Output::output_P;

using namespace std;

// Constructors

	Output::Output(PROPOSALParticle *part){

                version="Muon Propagation Code in Java v. 1.5.5";

                AMASIM=false;
                // outf = false; <--- currently it is const static
		this->p = part;
                DEBUG = false;
		HIST = -1;
		OUTNUM=8;
		c.resize(16);
                //z0='0', zm='-', zd='.', ze='e';
		POWOUT=(long)roundValue(pow((double) 10, 8));
                //I3flag=false;
                //raw = false;
                texi = false;
                I3hist.resize(HISTSIZE);
                inf=false;
                outf=false;
	
	};

//	Output::Output( ){
//                AMASIM=false;
//                // outf = false; <--- currently it is const static
//			DEBUG = false;
//		HIST = -1;
//		OUTNUM=8;
//			c.resize(16);
//                //z0='0', zm='-', zd='.', ze='e';
//		POWOUT=(long)roundValue(pow((double) 10, 8));
//                //I3flag=false;
//                //raw = false;
//                texi = false;
//                I3hist.resize(HISTSIZE);
//                inf=false;
//                outf=false;
//	}

// Memberfunctions

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
            catch(MMCexception& error)
            {
                cout<<error.what()<<endl;
                return false;
            }

        }

        void Output::setStdout(){
            out.str("");//<---NOT IN JAVA
            out.clear(); //<---NOT IN JAVA
            stdout_backup = cout.rdbuf();
            cout.rdbuf( out.rdbuf() );

        }



        string Output::getStdout(){
            cout.flush();
            cout.rdbuf( stdout_backup );
            return out.str();

        }
/*
	Output::output(int wint, String comp, double de, double ef){
	
		String str;
		history.open(Output.txt);
		err.open(Error.txt);
	
		if(DEBUG || (HIST==1)){
		    str=" . ";
		    if(wint==0){
			if(comp=="conti") str+="lost continuously "+f(de)+" MeV over "+f(ef)+" cm up to "+f(p.r)+" cm\n";
			else str+="muon with energy "+f(de)+" MeV traveled  "+f(ef)+" cm\n";
	 	   }
	 	   else{
			if(wint==1) str+="decayed into "+comp;
			else if(wint==2) str+="ionized by "+comp;
			else if(wint==3) str+="bremsed by "+comp;
			else if(wint==4) str+="photoed by "+comp;
			else if(wint==5) str+="epaired by "+comp;
			str+=" ... lost "+f(de)+" MeV, ef = "+f(ef)+" MeV, xf = "+f(p.r)+" cm\n";
		    }
	 	   if(HIST==1) history>>(str)>>endl;
		    if(DEBUG) err>>(str)>>endl; 
		}

		if(HIST==2 || (RecDec && wint==1)){
	    		if(wint==0 || wint==1) str=comp;
	    		else if(wint==2) str="delta";
	    		else if(wint==3) str="brems";
	    		else if(wint==4) str="munu";
	    		else if(wint==5) str="epair";
	    		else str="";
	    		if(str!=""){
				gens++;
				if(I3flag) 
					I3hist.addElement(new PROPOSALParticle(igen, gens, str, p.x, p.y, p.z, p.theta, p.phi, de, p.t, wint==0?ef:0));
			
				else{
		   		 	str="TR "+gens+" "+igen+" "+str+" "+f(p.x*1.e-2)+" "+f(p.y*1.e-2)+" "+f(p.z*1.e-2);
		    			str+=" "+f(180-p.theta)+" "+f(p.phi<180?p.phi+180:p.phi-180)+" "+f(wint==0?ef*1.e-2:0)+" "+f(de*1.e-3)+" "+f(p.t*1.e9)+"\n";
			
		     			history<<str<<endl;
				}
	    		}
		}
	history.close();
	err.close();

	}

*/
	// -------------------------------------------------------------------------------- //

	double Output::read(){

            //ReadWrite.exceptions ( ifstream::failbit | ifstream::badbit );
            double data;
            if(raw)
            {
                if(ReadWrite.good())
                {
                    ReadWrite.read( reinterpret_cast<char*>( &data ), sizeof data );
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

                // -------------------------------------------------------------------------------- //

        void Output::write(double g)
        {
            if(ReadWrite.good())
            {
                if(raw)
                {
                    ReadWrite.write(reinterpret_cast <char*> (&g), sizeof( g ));
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



        void Output::open(std::string name)
        {
            inf=false;
            outf=false;
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
                inf=true;
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
                outf=true;
                cerr<<"Parametrization tables will be saved to the file "<<name<<endl;
            }


        }


//	std::string Output::f(double d){

            
//	char i, k;
//	int log, dot;
//	long aux, res;
//	double abs;
//	bool end;

//        if(d>0)
//        {
//            k=0;
//            abs=d;
//        }
//        else if(d==0)
//        {
//            return "0";
//        }
//        else
//        {
//            c[0]=zm;
//            k=1;
//            abs=-d;
//        }

//        log=(int)roundValue(std::log(abs)/LOG10);
//	double base = 10;
//	res=(long)roundValue(abs*pow(base, OUTNUM-1-log));
//        cout<<log<<"\t"<<res<<endl;
//        if(res>=POWOUT)
//        {
//            res/=10;
//            log++;
//        }

//        if(log<=-2 || log>=OUTNUM)
//        {
//            dot=k+1;
//            end=true;
//        }
//        else
//        {
//            dot=k+log+1;
//            end=false;
//        }

//	k+=OUTNUM;
//        cout<<k<<endl;
//        for(i=0; i<=OUTNUM; i++)
//        {
//            if(k-i==dot)
//            {
//                c[dot]=zd;
//            }
//            else
//            {
//                c[k-i]=(char)(z0+res%10);
//                res/=10;

//            }
//        }

//	while(c[k]==z0) k--;
//	if(c[k]==zd) k--;
//	k++;

//	if(end){
//            c[k]=ze;
//            k++;
//            if(log<0)
//            {
//                c[k]=zm;
//                k++;
//                log=-log;
//            }
//            if(log<10)
//            {
//                aux=1;
//            }
//            else if(log<100)
//            {
//                aux=2;
//            }
//            else
//            {
//                aux=3;
//            }
//	    k+=aux-1;
//            for(i=0; i<aux; i++)
//            {
//                c[k-i]=(char)(z0+log%10);
//                log/=10;
//            }
//	    k++;
//	}

//	string ResultString;
//        for (int K = 0; K<=k ; K++)
//	{

//		ResultString.append(1,c[K]);
//	}

//	return ResultString;
//	}

        string Output::f(double d){

            stringstream ss;
            ss.precision(OUTNUM);
            ss<<scientific<<d;
            return ss.str();

        }

	void Output::output(int wint, string comp, double de, double ef){
	string str;
        //cout << wint <<"\t"<< de <<endl;
	if(DEBUG || (HIST==1)){
	    str=" . ";
	    if(wint==0){
		if(comp=="conti") str+="lost continuously "+f(de)+" MeV over "+f(ef)+" cm up to "+f(p->r)+" cm\n";
		else str+="muon with energy "+f(de)+" MeV traveled  "+f(ef)+" cm\n";
	    }
	    else{
		if(wint==1) str+="decayed into "+comp;
		else if(wint==2) str+="ionized by "+comp;
		else if(wint==3) str+="bremsed by "+comp;
		else if(wint==4) str+="photoed by "+comp;
		else if(wint==5) str+="epaired by "+comp;
		str+=" ... lost "+f(de)+" MeV, ef = "+f(ef)+" MeV, xf = "+f(p->r)+" cm\n";
	    }
	    if(HIST==1) history<<str<<endl;
	    if(DEBUG) cerr<<str<<endl;
	}
	if(HIST==2 || (RecDec && wint==1)){
	    if(wint==0 || wint==1) str=comp;
	    else if(wint==2) str="delta";
	    else if(wint==3) str="brems";
	    else if(wint==4) str="munu";
	    else if(wint==5) str="epair";
	    else str="";
	    if(str!=""){
		gens++;
		if(I3flag){
			PROPOSALParticle *PROPOSALParticleToPrint = new PROPOSALParticle(igen, gens, str, p->x, p->y, p->z, p->theta, p->phi, de, p->t, wint==0?ef:0);
                        I3hist.push_back(PROPOSALParticleToPrint);
                        //delete PROPOSALParticleToPrint;

		}
		else{
		    	std::ostringstream Print_igen;
			std::ostringstream Print_gens;
			Print_gens <<gens;
			Print_igen <<igen;
		    str="TR "+Print_gens.str()+" "+Print_igen.str()+" "+str+" "+f(p->x*1.e-2)+" "+f(p->y*1.e-2)+" "+f(p->z*1.e-2);
		    str+=" "+f(180-p->theta)+" "+f(p->phi<180?p->phi+180:p->phi-180)+" "+f(wint==0?ef*1.e-2:0)+" "+f(de*1.e-3)+" "+f(p->t*1.e9)+"\n";
		    history<<str<<endl;
		}
	    }
	}
    }

        // split the string into substrings using given delimeter list
        deque<string> Output::splitString(string args, string Delimiters){

            deque<string> Tokens;
            string::size_type Next, EndNext=0;

            while (EndNext != string::npos) {
                // Find next token
                Next = args.find_first_not_of(Delimiters, EndNext);
                EndNext = args.find_first_of(Delimiters, Next);

                // Stop if end of string reached
                if (Next == string::npos) break;

                // Add token to vector.
                Tokens.push_back(args.substr(Next, EndNext - Next));
            }
            return Tokens;
        }


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
                // -------------------------------------------------------------------------------- //

        void Output::init(string name){

            HIST=0;
            p->location(name, 0, 0, 0, 0, 0, 0);
            // history.setLength(0);

        }

        void Output::initDefault(string type){
            HIST=1;
            p->location(type, 0, 0, 0, 0, 0, 0);

            // Figure out what to do with the line below

            //history.setLength(0);
        }

        void Output::initDefault(int igen, int gens, string type, double time, double x, double y, double z, double theta, double phi){
            HIST=3;
            this->igen=igen;
            this->gens=gens;
            p->location(type, time, x, y, z, theta, phi);

            // Figure out what to do with the line below

            //history.setLength(0);

            if(RecDec)
            {
                if(I3flag)
                {
                    I3hist.clear();
                }
            }

        }

        void Output::initF2000(int igen, int gens, string type, double time, double x, double y, double z, double theta, double phi)
        {
            HIST=2;
            this->igen=igen;
            this->gens=gens;
            p->location(type, time, x, y, z, theta, phi);


            // Figure out what to do with the line below

            //history.setLength(0);
            if(I3flag)
            {
                I3hist.clear();
            }
            //cout<<"still alive in Output::initF2000\t"<<endl;

        }


        void Output::close()
        {
            inf=false;
            outf=false;
            ReadWrite.close();

        }

        void Output::Delete(string name){

            Output::close();
            cout<<"... Deleting corrupt file "<<endl;
            if( remove( name.c_str() ) != 0 )
              cout<< "Error deleting file" <<endl;
            else
              cout<< "File successfully deleted" <<endl;
        }


        void Output::PROPOSALParticleToFileAscii(PROPOSALParticle *p, string filename, bool isPrimary)
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
            if(isPrimary){

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



        void Output::PROPOSALParticleToFileRaw(PROPOSALParticle *p, string filename, bool isPrimary)
        {
            char c='#'; // primary identifier
            char d='&'; // secondary identifier
            char n[10];

            for(size_t i =0;i<p->name.length();i++)
            {
                n[i]=p->name.at(i);
            }
            for(size_t j=p->name.length();j<10;j++)
            {
                n[j]='\0';
            }

            if(FileExist(filename))
            {
                output_P.open(filename.c_str(),ios_base::app|ios::binary);
            }
            else
            {
                output_P.open(filename.c_str(),ios::binary);
            }
            if(isPrimary){


                output_P.write(&c, sizeof( c ));
                output_P.write(reinterpret_cast <char*> (&p->r), sizeof( p->r ));
                output_P.write(reinterpret_cast <char*> (&p->x), sizeof( p->x ));
                output_P.write(reinterpret_cast <char*> (&p->y), sizeof( p->y ));
                output_P.write(reinterpret_cast <char*> (&p->z), sizeof( p->z ));
                output_P.write(reinterpret_cast <char*> (&p->t), sizeof( p->t ));
                output_P.write(reinterpret_cast <char*> (&p->l), sizeof( p->l ));
                output_P.write(reinterpret_cast <char*> (&p->theta), sizeof( p->theta ));
                output_P.write(reinterpret_cast <char*> (&p->phi), sizeof( p->phi ));
                output_P.write(reinterpret_cast <char*> (&p->e), sizeof( p->e ));
                output_P.write(reinterpret_cast <char*> (&p->type), sizeof( p->type ));
                output_P.write(reinterpret_cast <char*> (&p->xi), sizeof( p->xi ));
                output_P.write(reinterpret_cast <char*> (&p->yi), sizeof( p->yi ));
                output_P.write(reinterpret_cast <char*> (&p->zi), sizeof( p->zi ));
                output_P.write(reinterpret_cast <char*> (&p->Ei), sizeof( p->Ei ));
                output_P.write(reinterpret_cast <char*> (&p->xc), sizeof( p->xc ));
                output_P.write(reinterpret_cast <char*> (&p->yc), sizeof( p->yc ));
                output_P.write(reinterpret_cast <char*> (&p->zc), sizeof( p->zc ));
                output_P.write(reinterpret_cast <char*> (&p->Ec), sizeof( p->Ec ));
                output_P.write(reinterpret_cast <char*> (&p->xf), sizeof( p->xf ));
                output_P.write(reinterpret_cast <char*> (&p->yf), sizeof( p->yf ));
                output_P.write(reinterpret_cast <char*> (&p->zf), sizeof( p->zf ));
                output_P.write(reinterpret_cast <char*> (&p->Ef), sizeof( p->Ef ));
                output_P.write(reinterpret_cast <char*> (&p->Elost), sizeof( p->Elost ));
                output_P.write(reinterpret_cast <char*> (&n), sizeof( n ));

            }
            else
            {
                output_P.write(&d, sizeof( d ));
                output_P.write(reinterpret_cast <char*> (&p->x), sizeof( p->x ));
                output_P.write(reinterpret_cast <char*> (&p->y), sizeof( p->y ));
                output_P.write(reinterpret_cast <char*> (&p->z), sizeof( p->z ));
                output_P.write(reinterpret_cast <char*> (&p->l), sizeof( p->l ));
                output_P.write(reinterpret_cast <char*> (&p->t), sizeof( p->t ));
                output_P.write(reinterpret_cast <char*> (&p->theta), sizeof( p->theta ));
                output_P.write(reinterpret_cast <char*> (&p->phi), sizeof( p->phi ));
                output_P.write(reinterpret_cast <char*> (&p->e), sizeof( p->e ));
                output_P.write(reinterpret_cast <char*> (&p->type), sizeof( p->type ));
                output_P.write(reinterpret_cast <char*> (&n), sizeof( n ));


            }
            output_P.close();
        }
