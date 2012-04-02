/*
 * PROPOSALParticle.cxx
 *
 *  Created on: 21.06.2010
 *      Author: koehne
 */


#include "methods.h"
#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include <exception>

#include "PROPOSALParticle.h"





using namespace std;

        /// I actually needed to put in the "default" settings into the Constructor.
        /// Those default settings maybe overwritten if you use a special Constructor. e.g. store-Constructor.

	PROPOSALParticle::PROPOSALParticle(){

		
		r=0, x=0, y=0, z=0, t=0;             // coordinates [cm] and age [sec]
		theta=0, phi=0;                      // zenith and azimuth of the momentum in [deg]
                costh=1., sinth=0., cosph=1., sinph=0.;  // cos and sin of theta and phi
		p=0, p2=0, e=0;                      // momentum, momentum square and energy in [MeV]
                m=MMU, l=LMU, c=1;                   // mass [MeV], lifetime [sec] and charge
		name="mu";                           // name of the particle
		low=m;                               // energy below which the particle is lost [MeV]
		type=1;                                 // particle type: 1 for muon, 2 for tau, 3 for electron
		igen=0, gens=1;                         // parent particle id, particle id

		xi=0, yi=0, zi=0, ti=0, Ei=0;        // coordinates at entry point [m,m,m,sec,GeV]
		xf=0, yf=0, zf=0, tf=0, Ef=0;        // coordinates at exit point
		xc=0, yc=0, zc=0, tc=0, Ec=0;        // coordinates at point of closest approach
		Elost=0;                             // energy lost in the detector volume [GeV]


		jt = false;
		df = false;


	}

    /*!
     * initialize particle
     */

    PROPOSALParticle::PROPOSALParticle(Propagate *pr, string name){

			
		r=0, x=0, y=0, z=0, t=0;             // coordinates [cm] and age [sec]
		theta=0, phi=0;                      // zenith and azimuth of the momentum in [deg]
                costh=1., sinth=0., cosph=1., sinph=0.;  // cos and sin of theta and phi
		p=0, p2=0, e=0;                      // momentum, momentum square and energy in [MeV]
                m=MMU, l=LMU, c=1;                   // mass [MeV], lifetime [sec] and charge
		
		low=m;                               // energy below which the particle is lost [MeV]
		type=1;                                 // particle type: 1 for muon, 2 for tau, 3 for electron
		igen=0, gens=1;                         // parent particle id, particle id

		xi=0, yi=0, zi=0, ti=0, Ei=0;        // coordinates at entry point [m,m,m,sec,GeV]
		xf=0, yf=0, zf=0, tf=0, Ef=0;        // coordinates at exit point
		xc=0, yc=0, zc=0, tc=0, Ec=0;        // coordinates at point of closest approach
		Elost=0;                             // energy lost in the detector volume [GeV]


		jt = false;
		df = false;

                //cout<<"In PROPOSALParticle constuctor: values are setted, find out PROPOSALParticle name"<<endl;

	if(StartsWith(name,"tau")){
	    this->name="tau";
	    type=2;
            m=MTAU;
            l=LTAU;
	}
	else if(StartsWith(name,"e")){
	    this->name="e";
	    type=3;
            m=ME;
	    l=-1;
	}
	else if(StartsWith(name,"monopole")){
	    type=4;
	    try{
		m=strtod(name.substr(9).c_str(),NULL)*1.e3;
	    }catch(exception& e){
		m=0;
	    }
            if(m<=0) m=MMON;
	    this->name="monopole-"+output->f(m*1.e-3);
            c=CMON;
	    l=-1;
	}
	else if(StartsWith(name,"stau")){
	    type=5;
	    try{
		m=strtod(name.substr(5).c_str(),NULL)*1.e3;
	    }catch(exception& e){
		m=0;
	    }
            if(m<=0) m=MSTAU;
	    this->name="stau-"+output->f(m*1.e-3);
            l=LSTAU;
	}
	else{
	    this->name="mu";
	    type=1;
            m=MMU;
            l=LMU;
	}

        //cout<<"In PROPOSALParticle constructor: PROPOSALParticle name is\t"<<this->name<<endl;
	low=m; e=m;
        this->propagate_=pr;

        if(PhysicsModel::get_elow()>low){
            low=PhysicsModel::get_elow();
        }
        if(PhysicsModel::get_ebig()<100*low){
            PhysicsModel::set_ebig(max(100*low, BIGENERGY));
        }
        //cout<<"In PROPOSALParticle constructor: initalize scattering"<<endl;
        scattering_ = new Scattering(this);
        //cout<<"In PROPOSALParticle constructor: initalize integral"<<endl;
        integral_ = new Integral(IROMB, IMAXS, IPREC2);
    }

    //----------------------------------------------------------------------------------------------------//

    /*!
     * store particle information
     */

    PROPOSALParticle::PROPOSALParticle(int Igen, int Gens, string Name, double X, double Y, double Z, double Theta, double Phi, double E, double T, double R, PROPOSALParticle *P){
		
		r=0, x=0, y=0, z=0, t=0;             // coordinates [cm] and age [sec]
		theta=0, phi=0;                      // zenith and azimuth of the momentum in [deg]
                costh=1., sinth=0., cosph=1., sinph=0.;  // cos and sin of theta and phi
		p=0, p2=0, e=0;                      // momentum, momentum square and energy in [MeV]
                m=MMU, l=LMU, c=1;                   // mass [MeV], lifetime [sec] and charge
		name="mu";                           // name of the particle
		low=m;                               // energy below which the particle is lost [MeV]
		type=1;                                 // particle type: 1 for muon, 2 for tau, 3 for electron
		igen=0, gens=1;                         // parent particle id, particle id

		xi=0, yi=0, zi=0, ti=0, Ei=0;        // coordinates at entry point [m,m,m,sec,GeV]
		xf=0, yf=0, zf=0, tf=0, Ef=0;        // coordinates at exit point
		xc=0, yc=0, zc=0, tc=0, Ec=0;        // coordinates at point of closest approach
		Elost=0;                             // energy lost in the detector volume [GeV]


		jt = false;
		df = false;


	PROPOSALParticle(Name, X, Y, Z, Theta, Phi, E, T);
	this->r=R;
	this->igen=Igen;
	this->gens=Gens;
	if(P!=NULL){ // orignial p != null..
	    xi=P->xi; yi=P->yi; zi=P->zi; ti=P->ti; Ei=P->Ei;
	    xf=P->xf; yf=P->yf; zf=P->zf; tf=P->tf; Ef=P->Ef;
	    xc=P->xc; yc=P->yc; zc=P->zc; tc=P->tc; Ec=P->Ec;
	    Elost=P->Elost;
	}
    }


    //----------------------------------------------------------------------------------------------------//

    /*!
     * store particle information
     */

    PROPOSALParticle::PROPOSALParticle(int Igen, int Gens, string Name, double X, double Y, double Z, double Theta, double Phi, double E, double T, double R){
	
			
		r=0, x=0, y=0, z=0, t=0;             // coordinates [cm] and age [sec]
		theta=0, phi=0;                      // zenith and azimuth of the momentum in [deg]
                costh=1., sinth=0., cosph=1., sinph=0.;  // cos and sin of theta and phi
		p=0, p2=0, e=0;                      // momentum, momentum square and energy in [MeV]
                m=MMU, l=LMU, c=1;                   // mass [MeV], lifetime [sec] and charge
		name="mu";                           // name of the particle
		low=m;                               // energy below which the particle is lost [MeV]
		type=1;                                 // particle type: 1 for muon, 2 for tau, 3 for electron
		igen=0, gens=1;                         // parent particle id, particle id

		xi=0, yi=0, zi=0, ti=0, Ei=0;        // coordinates at entry point [m,m,m,sec,GeV]
		xf=0, yf=0, zf=0, tf=0, Ef=0;        // coordinates at exit point
		xc=0, yc=0, zc=0, tc=0, Ec=0;        // coordinates at point of closest approach
		Elost=0;                             // energy lost in the detector volume [GeV]


		jt = false;
		df = false;

        string name=Name.length()==0?"?":Name[0]=='a'?Name.substr(1):Name;
        if(name.compare("tau")==0 || name.compare("tau-")==0 || name.compare("tau+")==0){
            if(name.compare("tau+")==0) type=-33;
            else type=-34;
            m=MTAU;
            l=LTAU;
        }
        else if(name.compare("mu")==0 || name.compare("mu-")==0 || name.compare("mu+")==0){
            if(name.compare("mu+")==0) type=-5;
            else type=-6;
            m=MMU;
            l=LMU;
        }
        else if(StartsWith(name,"stau") ){
            if((int)(name.find("stau+"))!=-1) type=-9131; else type=-9132;
            try{
                m=strtod(name.substr(5).c_str(),NULL)*1.e3;
            }catch(exception &e){
                m=0;
            }
            if(m<=0) m=MSTAU;
            l=LSTAU;
        }
        else if(name.compare("e")==0 || name.compare("e-")==0 || name.compare("e+")==0){
            if(name.compare("e+")==0) type=-2;
            else type=-3;
            m=ME;
            l=-1;
        }
        else if((int)(name.find("nu_"))!=-1){
            if(name.compare("nu_e")==0) type=-201;
            else if(name.compare("~nu_e")==0) type=-204;
            else if(name.compare("nu_mu")==0) type=-202;
            else if(name.compare("~nu_mu")==0) type=-205;
            else if(name.compare("nu_tau")==0) type=-203;
            else if(name.compare("~nu_tau")==0) type=-206;
            m=0;
            l=-1;
        }
        else{
            if(name.compare("delta")==0) type=-1002;
            else if(name.compare("brems")==0) type=-1001;
            else if(name.compare("munu")==0) type=-1004;
            else if(name.compare("epair")==0) type=-1003;
            else if(name.compare("hadr")==0) type=-1006;
            else if(name.compare("conti")==0) type=-1111;
            else type=0;
            m=0;
            l=0;
        }
        low=m;
        setEnergy(E);
        location(Name, T, X, Y, Z, Theta, Phi);
	this->r=R;
	this->igen=Igen;
	this->gens=Gens;
    }


    //----------------------------------------------------------------------------------------------------//

    /*!
     * store particle information
     */

    PROPOSALParticle::PROPOSALParticle(string aname, double X, double Y, double Z, double Theta, double Phi, double E, double T){


			
		r=0, x=0, y=0, z=0, t=0;             // coordinates [cm] and age [sec]
		theta=0, phi=0;                      // zenith and azimuth of the momentum in [deg]
                costh=1., sinth=0., cosph=1., sinph=0.;  // cos and sin of theta and phi
		p=0, p2=0, e=0;                      // momentum, momentum square and energy in [MeV]
                m=MMU, l=LMU, c=1;                   // mass [MeV], lifetime [sec] and charge
		name="mu";                           // name of the particle
		low=m;                               // energy below which the particle is lost [MeV]
		type=1;                                 // particle type: 1 for muon, 2 for tau, 3 for electron
		igen=0, gens=1;                         // parent particle id, particle id

		xi=0, yi=0, zi=0, ti=0, Ei=0;        // coordinates at entry point [m,m,m,sec,GeV]
		xf=0, yf=0, zf=0, tf=0, Ef=0;        // coordinates at exit point
		xc=0, yc=0, zc=0, tc=0, Ec=0;        // coordinates at point of closest approach
		Elost=0;                             // energy lost in the detector volume [GeV]


		jt = false;
		df = false;


	string name=aname.length()==0?"?":aname[0]=='a'?aname.substr(1):aname;
        if(name.compare("tau")==0 || name.compare("tau-")==0 || name.compare("tau+")==0){
            if(name.compare("tau+")==0) type=-33;
	    else type=-34;
            m=MTAU;
            l=LTAU;
	}
        else if(name.compare("mu")==0 || name.compare("mu-")==0 || name.compare("mu+")==0){
            if(name.compare("mu+")==0) type=-5;
	    else type=-6;
            m=MMU;
            l=LMU;
	}
	else if(StartsWith(name,"stau") ){
	    if((int)(name.find("stau+"))!=-1) type=-9131; else type=-9132;
	    try{
		m=strtod(name.substr(5).c_str(),NULL)*1.e3;
	    }catch(exception &e){
		m=0;
	    }
            if(m<=0) m=MSTAU;
            l=LSTAU;
	}
        else if(name.compare("e")==0 || name.compare("e-")==0 || name.compare("e+")==0){
            if(name.compare("e+")==0) type=-2;
	    else type=-3;
            m=ME;
            l=-1;
	}
	else if((int)(name.find("nu_"))!=-1){
            if(name.compare("nu_e")==0) type=-201;
            else if(name.compare("~nu_e")==0) type=-204;
            else if(name.compare("nu_mu")==0) type=-202;
            else if(name.compare("~nu_mu")==0) type=-205;
            else if(name.compare("nu_tau")==0) type=-203;
            else if(name.compare("~nu_tau")==0) type=-206;
	    m=0;
	    l=-1;
	}
	else{
            if(name.compare("delta")==0) type=-1002;
            else if(name.compare("brems")==0) type=-1001;
            else if(name.compare("munu")==0) type=-1004;
            else if(name.compare("epair")==0) type=-1003;
            else if(name.compare("hadr")==0) type=-1006;
            else if(name.compare("conti")==0) type=-1111;
	    else type=0;
	    m=0;
	    l=0;
	}
	low=m;
	setEnergy(E);
	location(aname, T, X, Y, Z, Theta, Phi);
    }

    //----------------------------------------------------------------------------------------------------//

    /*!
     * initialize the location and direction of the particle, time in sec, x, y, z in cm, theta and phi in deg
     */

    void PROPOSALParticle::location(string name, double time, double x, double y, double z, double theta, double phi){
	this->name=name;
	r=0;
	t=time;
	this->x=x;
	this->y=y;
	this->z=z;
	this->theta=theta;
	this->phi=phi;
        theta*=(PI/180);
        phi*=(PI/180);
	costh=cos(theta);
	sinth=sin(theta);
	cosph=cos(phi);
	sinph=sin(phi);
    }

    //----------------------------------------------------------------------------------------------------//

    /*!
     * advances the pa>rticle by the given distance
     */

  void PROPOSALParticle::advance(double dr, double ei, double ef){
       //double aux=0;
      long double ax=0, ay=0, az=0;
      long double tho=0, max_=0, rnd1=0, rnd2=0, sx=0, sy=0, sz=0, tx=0, ty=0, tz=0;

	r+=dr;
//        cout<<"PROPOSALParticle::advance starts: "<<"propagate_->get_exactTime() = "<<propagate_->get_exactTime()<<endl;
        if(propagate_->get_exactTime())
        {
            t+=getdt(ei, ef)/propagate_->get_rho();
        }
        else
        {
            t+=dr/SPEED;
        }

//cout<<"PROPOSALParticle::advance: "<<"propagate_->get_molieScat() = "<<propagate_->get_molieScat()<<endl;

        if(propagate_->get_molieScat()){

            tho=(long double)scattering_->gettho(dr, ei, ef);
            //tho=1e-6;
//cout<<"PROPOSALParticle::advance:401 "<<dr<<"\t"<<ei<<"\t"<<ef<<"\t"<<endl;
            max_=1/SQRT2;
            rnd1=(long double)propagate_->get_Standard()->sndrn(RandomDouble(), 0, tho, -max_, max_, false);
            rnd2=(long double)propagate_->get_Standard()->sndrn(RandomDouble(), 0, tho, -max_, max_, false);
//cout<<"PROPOSALParticle::advance:405 "<<rnd1<<"\t"<<rnd2<<"\t"<<endl;

            sx=(rnd1/SQRT3+rnd2)/2;
	    tx=rnd2;
//cout<<"PROPOSALParticle::advance:409 "<<rnd1<<"\t"<<rnd2<<"\t"<<endl;
            rnd1=(long double)propagate_->get_Standard()->sndrn(RandomDouble(), 0, tho, -max_, max_, false);
            rnd2=(long double)propagate_->get_Standard()->sndrn(RandomDouble(), 0, tho, -max_, max_, false);
//cout<<"In PROPOSALParticle::advance:412"<<rnd1<<"\t"<<rnd2<<"\t"<<endl;


            sy=(rnd1/SQRT3+rnd2)/2;
	    ty=rnd2;


                sz=sqrt(max(1.-(sx*sx+sy*sy), (long double)0.));

                tz=sqrt(max(1.-(tx*tx+ty*ty), (long double)0.));


//cout<<sx<<"\t"<<sy<<"\t"<<sz<<"\t";
            ax=sinth*cosph*sz+costh*cosph*sx-sinph*sy;
            ay=sinth*sinph*sz+costh*sinph*sx+cosph*sy;
            az=costh*sz-sinth*sx;

	    x+=ax*dr;
	    y+=ay*dr;
	    z+=az*dr;

//
//           if(x>5){
//            cout<<cosph<<"\t"<<sinph<<"\t"<<costh<<"\t"<<sinth<<"\t"<<ax<<"\t"<<ay<<"\t"<<az<<"\t"<<dr<<"\t"<<sx<<"\t"<<sy<<"\t"<<sz<<"\t"<<rnd1<<"\t"<<rnd2<<"\t"<<x<<endl;
//           }

            ax=sinth*cosph*tz+costh*cosph*tx-sinph*ty;
            ay=sinth*sinph*tz+costh*sinph*tx+cosph*ty;
            az=costh*tz-sinth*tx;



//cout<<ax<<"\t"<<ay<<"\t"<<tz<<"\t";
            costh= az;
            //sinth=sqrt(1-costh)*sqrt(costh +1);
            sinth=sqrt(max(1-costh*costh, (long double)0));
	    if(sinth!=0){
                sinph=ay/sinth;
                cosph=ax/sinth;
	    }

            //theta=acos(costh>1?1:costh<-1?-1:costh)*180/Pi_;
            if(costh>1)
            {
                theta=acos(1)*180./PI;
            }
            else if(costh<-1)
            {
                theta=acos(-1)*180./PI;
            }
            else
            {
                theta = acos(costh)*180./PI;
            }

            //phi=acos(cosph>1?1:cosph<-1?-1:cosph)*180/Pi_;
            if(cosph>1)
            {
                phi=acos(1)*180./PI;
            }
            else if(cosph<-1)
            {
                phi=acos(-1)*180./PI;
            }
            else
            {
                phi=acos(cosph)*180./PI;
            }

            if(sinph<0)
            {
                phi=360.-phi;
            }
            if(phi>=360)
            {
                phi-=360.;
            }
            //cout<<tho<<"\t"<<cosph<<"\t"<<sinph<<"\t"<<phi<<"\t"<<theta<<"\t"<<costh<<"\t"<<sinth<<"\t"<<ax<<"\t"<<ay<<"\t"<<az<<"\t"<<dr<<"\t"<<sx<<"\t"<<sy<<"\t"<<sz<<"\t"<<rnd1<<"\t"<<rnd2<<"\t"<<x<<endl;
            //cout<<sinph<<endl;
//cout<<sinth<<"\t"<<costh<<endl;
//            cout<<sinth<<endl;//"\t"<<costh<<"\t"<<sinph<<"\t"<<cosph<<endl;

	}
	else{
	    x+=sinth*cosph*dr;
	    y+=sinth*sinph*dr;
	    z+=costh*dr;
	}

    }

    //----------------------------------------------------------------------------------------------------//

    /*!
     * sets the energy of the particle
    */

    void PROPOSALParticle::setEnergy(double e){
	this->e=e;
	p2=e*e-m*m;
        p=sqrt(max(p2,0.0));
    }

    //----------------------------------------------------------------------------------------------------//

    /*!
     * function for time delta calculation - interface to Integral
     */

    //IMPORTANT: the ouput lines need to be commented in when class output is running!!!!!!!!!!!!!!!!!!!!!!!!!!!

    double PROPOSALParticle::function(double E){
	double aux;
        aux=propagate_->get_cros()->function(E);
        aux*=e/(p*SPEED);

	return aux;
    }

    //----------------------------------------------------------------------------------------------------//

    /*!
     * time delta, corresponding to the given propagation distance
     */

    double PROPOSALParticle::getdt(double ei, double ef){
	if(jt){
            if(abs(ei-ef)>abs(ei)*HALF_PRECISION)
            {
                double aux=interpolateJ_->interpolate(ei);
                double aux2=aux-interpolateJ_->interpolate(ef);
                if(abs(aux2)>abs(aux)*HALF_PRECISION)
                {
                    return aux2;
                }
	    }
            return interpolateJdf_->interpolate((ei+ef)/2)*(ef-ei);
	}
        else{
            cout<<"PROPOSALParticle::getdt: "<<"ei= "<<ei<<"\t ef= "<<ef<<endl;
            return integral_->integrateWithLog(ei, ef, this);
        }
    }

    //----------------------------------------------------------------------------------------------------//



    /*!
     * 1d parametrization - interface to Interpolate
     */

    double PROPOSALParticle::functionInt(double e){
//cout<<"PROPOSALParticle: FunctionInt"<<endl;
        if(df){
//cout<<"PROPOSALParticle: FunctionInt if"<<endl;
            return function(e);

        }
        else{
//cout<<"PROPOSALParticle: FunctionInt else"<<endl;
//cout<<I->integrateWithLog(e, low, this)<<endl;
            return integral_->integrateWithLog(e, low, this);
        }
    }
    //----------------------------------------------------------------------------------------------------//

    double PROPOSALParticle::alpha(double v){
        double alpha0 = 0.007297352533285885;
        double Q = v*e;
        return alpha0 /( 1- 2*alpha0/(3*PI)*log(Q*Q/ME));



    }

