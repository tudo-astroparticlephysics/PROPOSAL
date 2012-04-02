/*
 * Propagate.cxx
 *
 *  Created on: 05.07.2010
 *      Author: koehne
 */

#include "Propagate.h"
#include <math.h>
#include "methods.h"
#include "algorithm"
#include "stdlib.h"

#include "PROPOSALParticle.h"
#include "Medium.h"

#include "Energy2Loss.h"
#include "Output.h"
#include "StandardNormal.h"
#include "Integral.h"
#include "StandardNormal.h"


#include "IonizContinuous.h"
#include "IonizStochastic.h"


#include "EpairStochastic.h"
#include "EpairContinuous.h"


#include "BremsContinuous.h"
#include "BremsStochastic.h"

#include "Decay.h"

#include "Energy2LossX.h"
#include "Energy2LossE.h"

#include "Photonuclear.h"
#include "PhotoStochastic.h"
#include "PhotoContinuous.h"


#include "Interpolate.h"

int Propagate::g=5;

using namespace std;

    //----------------------------------------------------------------------------------------------------//

    /**
    * Default Constructor
    */

    Propagate::Propagate(){}

    //----------------------------------------------------------------------------------------------------//

    /**
     * initialize all classes necessary for propagation of a muon.
     */

    Propagate::Propagate(string w, double ecut, double vcut){

        init(w, ecut, vcut, "mu", 1.);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initialize all classes necessary for propagation of a muon or tau.
     */

    Propagate::Propagate(string w, double ecut, double vcut, string type){

        init(w, ecut, vcut, type, 1.);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initialize all classes necessary for propagation. string w contains the name of the medium. For the
     * definition of ecut and vcut see help for the class Medium constructor. Set type to "mu" or "tau".
     * The value of rho sets the multiplicative medium density correction factor.
     * To enable parametrization routines, call this->interpolate("all") after this class is created.
     */

    Propagate::Propagate(string w, double ecut, double vcut, string type, double rho){

                init(w, ecut, vcut, type, rho);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * class initializer
     */

    void Propagate::init(std::string w, double ecut, double vcut, std::string type, double rho){
        //cout<<"Init Propagate starts"<<endl;
        this->rho = 1.; // <--- to think about!!! NOT IN THE JAVA-CODE!!!!!
        sdec=false;
        recc=false;
        exactTime=false;
        contiCorr=false;
        molieScat=false;

        dw=false;
        rw=0, hw=0;

        df=false;
        jt=false;

        up = true;

        //cout<<"in init: values are setted. initalize physics objects"<<endl;
        particle_ = new PROPOSALParticle(this, type);
        //cout<<"particle created"<<endl;
        medium_ = new Medium(w, ecut, vcut, rho);
        //cout<<"Medium created"<<endl;
        cros = new CrossSections(particle_, medium_);
        E2Loss_ = new Energy2Loss(cros);
        //cout<<"in init: physics objects are initalized. initalize integral"<<endl;
        integral_.resize(2);
        for(int i=0; i<2; i++){
                integral_.at(i) = new Integral(IROMB, IMAXS, IPREC2);
        }

        StandardN = new StandardNormal(IROMB, IMAXS, IPREC);
        o = new Output(particle_);

    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Propagates the particle of initial energy e to the distance r. Returns the const energy if the
     * particle has survived or the track length to the point of disappearance with a minus sign otherwise.
     */

        double Propagate::propagateTo(double r, double e){
		int wint;
		bool DEBUG, flag;
                double ei, ef=0, efd, efi, aux=0, dr; //ini;
		double rndd, rndi, rnd1, rnd2, rnd3, rnddMin, rndiMin, rndTot;

		DEBUG=o->DEBUG;
  //cout<<"jt in Propagate: jt = "<< jt <<endl;

		if(o->HIST==-1){
                        o->init(particle_->name);
		}

		ei=e;
		ef=ei;
		if(r<0){
			r=0;
		}


                if(e<=particle_->low || r==0){
			flag=false;
		}
		else{
			flag=true;
		}
		if(DEBUG){

                        cerr<<"\nPropagating "<<particle_->name<<" of energy "<<o->f(ei)<<" MeV to a distance of "<<o->f(r)<<" cm"<<endl;
		}
                /////////////
                /////////////

//                fstream EnergyComparisonOut;
//                stringstream Dateiname;
//                if(jt){
//                    Dateiname << "/home/koehne/Interpolation_" << this->g << ".txt";
//                    EnergyComparisonOut.open(Dateiname.str().c_str(), ios::out | ios::app );
//                }
//                else{
//                    EnergyComparisonOut.open("/home/koehne/Integration_Romb.txt", ios::out | ios::app );
//                }

//                EnergyComparisonOut.precision(16);
//                int PrezZahler=1;
//                int PrezMax = 2;

                /////////////Tomas
                /////////////
                while(flag){
//                    cerr<<particle_->name<<endl;

//cout<<"Propagate::propagateTo::while starts"<<endl;
//                    if(PrezZahler==PrezMax)flag=false;
//                    PrezZahler++;
//                       EnergyComparisonOut << ei << endl;


			rndd=-log(RandomDouble());
			rndi=-log(RandomDouble());

			if(DEBUG){
				cerr<<"1. solving the tracking integral ...  rndd = "<<o->f(rndd)<<"  rndi = "<<o->f(rndi)<<" ...  "<<endl;
			}
                        if(particle_->l<0){
				rnddMin = 0;
			}
			else{

                                rnddMin = getpr(ei, rndd, false)/rho;

			}
			if(DEBUG){
				cerr<<" \t \t \t rnddMin = "<<o->f(rnddMin)<<" (d)  "<<endl;
			}
                        rndiMin=getpr(ei, rndi, true);

			if(DEBUG){
				cerr<<"rndiMin = "<<o->f(rndiMin)<<" (i)"<<endl;
			}

			if(DEBUG){
				cerr<<"2. evaluating the energy loss ...  "<<endl;
			}
			if(rndd>=rnddMin || rnddMin<=0){
                                efd=particle_->low;
			}
                        else{
				efd=getef(ei, rndd*rho, false);
			}
			if(DEBUG){
				cerr<<"efd = "<<o->f(efd)<<" MeV  "<<endl;
			}
			if(rndi>=rndiMin || rndiMin<=0){
                                efi=particle_->low;
			}
			else{
				efi=getef(ei, rndi, true);
			}
			if(DEBUG){
				cerr<<"efi = "<<o->f(efi)<<" MeV ...  "<<endl;
			}
			pint=(efi>efd);

			if(pint){
				ef=efi;
			}
			else{
				ef=efd;
			}

			if(DEBUG){
				cerr<<" \t \t \t lost "<<o->f(ei-ef)<<" MeV  ef = "+o->f(ef)+" MeV"<<endl;
			}

			if(DEBUG){
				cerr<<"3. calculating the displacement ...  "<<endl;
			}

                        dr=cros->getdx(ei, ef, rho*(r-particle_->r))/rho;
//cout<<"in propagateTo: "<<dr<<"\t"<<r<<"\t"<<particle_->r<<endl;
			if(DEBUG){
				cerr<<"dr = "<<o->f(dr)<<" cm"<<endl;
                        }
//cout<<"242"<<endl;
                        if(dr<r-particle_->r){
				if(DEBUG){
					cerr<<"4. calculating the local time ...  "<<endl;;
				}
			}
			else{
                                dr=r-particle_->r;
				if(DEBUG){
					cerr<<"4. getting the const energy ...  "<<endl;
				}
                                ef=cros->getef(ei, rho*dr);

				if(DEBUG){
					cerr<<"lost "<<o->f(ei-ef)<<" MeV  ef = "<<o->f(ef)+" MeV"<<endl;
				}
				if(DEBUG){
					cerr<<"5. calculating the local time ...  "<<endl;
				}
			}
//cout<<"262"<<endl;
			if(recc){
                                o->output(0, "a"+particle_->name, ei, dr);
			}

                        particle_->advance(dr, ei, ef);
//cout<<"268"<<endl;
			if(DEBUG){
                                cerr<<"t = "<<o->f(particle_->t)<<" s"<<endl;
			}
                        if(abs(r-particle_->r)<abs(r)*COMPUTER_PRECISION){
                                particle_->r=r;  // computer precision control
				}
//cout<<"275"<<endl;
			if(contiCorr){
                                if(ef!= particle_->low){
                                        ef=StandardN->sndrn(RandomDouble(), ef, sqrt(E2Loss_->e2le->dE2de(ei, ef)), particle_->low, ei, false);
				}

			}
//cout<<"282"<<endl;
                        if(recc){
				o->output(0, "conti", ei-ef, -dr);
			}
                        if(ef==particle_->low || particle_->r==r){
				break;
			}
//cout<<"289"<<endl;
			if(DEBUG){
				cerr<<"5. choosing the cross section ..."<<endl;
			}
			rnd2=RandomDouble();;
			rnd3=RandomDouble();;
                        particle_->setEnergy(ef);
//cout<<"pint = "<<pint<<endl;
			if(pint){
				rnd1=RandomDouble();;

				if(dw){
                                        if(particle_->r>hw){
						double exp=abs(rw);
						double pow_=pow(rnd2, exp);
						if(rw>0){
							rnd2=1-pow_*rnd2;
						}
						else{
							rnd2=pow_*rnd2;
						}

						rw=(1+exp)*pow_;
                                                hw=particle_->r;
						dw=false;
					}
				}
				if(DEBUG){
                                        decayS=cros->get_decay()->decay();
				}

                                ionizS=cros->get_ionization()->get_Stochastic()->dNdx(rnd2);
                                bremsS=cros->get_bremsstrahlung()->get_Stochastic()->dNdx(rnd2);
                                photoS=cros->get_photonuclear()->get_Stochastic()->dNdx(rnd2);
                                epairS=cros->get_epairproduction()->stochastic_->dNdx(rnd2);
				totalS=ionizS+bremsS+photoS+epairS;
				rndTot=rnd1*totalS;
//cout<<ionizS<<"\t"<<bremsS<<"\t"<<photoS<<"\t"<<epairS<<"\t"<<totalS<<endl;
				if(DEBUG){
					cerr<<" . rnd1 = "<<o->f(rnd1)<<" rnd2 = "<<o->f(rnd2)<<
								 " rnd3 = "<<o->f(rnd3)<<" decay = "<<o->f(decayS)<<endl;
				}
				if(DEBUG){
					cerr<<" . ioniz = "<<o->f(ionizS)<<" brems = "<<o->f(bremsS)<<
								 " photo = "<<o->f(photoS)<<" epair = "<<o->f(epairS);
				}
				if(ionizS>rndTot){

                                        aux=cros->get_ionization()->get_Stochastic()->e(rnd3);
					ef-=aux;
					wint=2;
				}
				else if(ionizS+bremsS>rndTot){
                                        aux=cros->get_bremsstrahlung()->get_Stochastic()->e(rnd3);
					ef-=aux;
					wint=3;
				}
				else if(ionizS+bremsS+photoS>rndTot){
                                        aux=cros->get_photonuclear()->get_Stochastic()->e(rnd3);
					ef-=aux;
					wint=4;
				}
				else if(ionizS+bremsS+photoS+epairS>rndTot){
                                        aux=cros->get_epairproduction()->stochastic_->e(rnd3);
					ef-=aux;
					wint=5;
				}
				else{ // due to the parameterization of the cross section cutoffs
					ei=ef;
					continue;
				}
                        }

			else{
                                if(particle_->type==2){
                                        aux=cros->get_decay()->e(rnd2, rnd3, RandomDouble(), o); ef=0; wint=1;
				}
				else{
                                        aux=cros->get_decay()->e(rnd2, rnd3, 0.5, o); ef=0; wint=1;
				}
			}

			if(wint==1){
                                o->output(wint, cros->get_decay()->get_out(), aux, ef);
			}
			else{
                                o->output(wint, medium_->get_E(cros->get_component()), aux, ef);
			}



                        if(ef<=particle_->low){

				break;
			}
			ei=ef;

		}
                ////////////Tomas
//                while(PrezZahler<PrezMax+1){
//                    EnergyComparisonOut << 0 << endl;
//                    PrezZahler++;
//                }
                //////////////

		if(sdec){
                        if(particle_->r!=r && ef!=0 && particle_->l>=0){
                                particle_->setEnergy(particle_->m);
                                particle_->t+=-particle_->l*log(RandomDouble());
                                if(particle_->type==2){
                                        aux=cros->get_decay()->e(RandomDouble(), 0.5, RandomDouble(), o);
				}
				else{
                                        aux=cros->get_decay()->e(RandomDouble(), 0.5, 0.5, o);
				}

				ef=0;
                                o->output(1, cros->get_decay()->get_out(), aux, ef);
			}
		}
                particle_->setEnergy(ef);  // to remember const state of the particle
		o->HIST=-1;  // to make sure user resets particle properties
//cout<<"while_counter="<<while_counter<<endl;

                if(particle_->r==r){
			if(DEBUG){
				cerr<<"PROPOSALParticle reached the border with energy ef = "<<o->f(ef)<<" MeV";
			}

			return ef;
		}
		else{
			if(DEBUG){
                                if(particle_->l<0){
                                        cerr<<"PROPOSALParticle stopped at rf = "<<o->f(particle_->r)<<" cm";
				}
				else{
                                        cerr<<"PROPOSALParticle disappeared at rf = "<<o->f(particle_->r)<<" cm";
				}
			}

                        return -particle_->r;
		}
    }

    //----------------------------------------------------------------------------------------------------//


    /**
     * Propagates the particle of initial energy e to the distance r. Returns the const energy if the
     * particle has survived or the track length to the point of disappearance with a minus sign otherwise.
     * Also calculates particle energy at point rc. Call getPropEc() to get this energy.
     */

	double Propagate::propagateTo(double r, double e, double rc){
		int HIST;
		double result;
		HIST=o->HIST;
		ec=propagateTo(rc, e);
                tc=particle_->t;
		if(ec>0){
			o->HIST=HIST;
			result=propagateTo(r, ec);
		}
		else result=ec;
		return result;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * returns the particle energy at rc if the particle has survived or the
     * distance from the point of decay to rc with a minus sign otherwise.
     */

	double Propagate::getPropEc(){
    	return ec;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * returns the particle time at rc if the particle has survived or at the point of decay otherwise.
     */

	double Propagate::getPropTc(){
		return tc;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * function for energy range calculation - interface to Integral
     */

	double Propagate::function(double E){
		const bool DEBUG=false;
		double aux;
                aux=cros->function(E);

		if(!pint){
                        decayS=cros->get_decay()->decay();
			if(DEBUG){
                                o->err<<" + "+o->f(particle_->e);
			}

			return aux*decayS;
		}
		else{
                        ionizS=cros->get_ionization()->get_Stochastic()->dNdx();
			if(DEBUG){
				o->err<<" \t "+o->f(ionizS);
			}
                        bremsS=cros->get_bremsstrahlung()->get_Stochastic()->dNdx();
			if(DEBUG){
				o->err<<" \t "+o->f(bremsS);
			}
                        photoS=cros->get_photonuclear()->get_Stochastic()->dNdx();
			if(DEBUG){
				o->err<<" \t "+o->f(photoS);
			}
                        epairS=cros->get_epairproduction()->stochastic_->dNdx();
			if(DEBUG){
				o->err<<" \t "+o->f(epairS);
			}
			totalS=ionizS+bremsS+photoS+epairS;

			return aux*totalS;
		}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * call this routine to enable interpolationcros-> To enable everything, set w="all"
     */

	void Propagate::interpolate(string w){


		int i;
                double e_hi=PhysicsModel::get_ebig();
                double e_low=particle_->low;
//cout<<"function interpolate starts"<<endl;
//!!!!!!!!!!!!!!!!!!!!the following line must be commented in when output is running!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                cout<<"Parameterizations apply in the energy range from "<<e_low<<" MeV to "<<e_hi<<" MeV"<<endl;
                if(w.find("all")!=string::npos){
			w+=" crs trackE trackX";
			if(exactTime){
				w+=" trackT";
			}
			if(contiCorr){
				w+=" contiR";
			}
			if(molieScat){
				w+=" molieS gaussR";
			}
		}
		if(w.find("crs")!=string::npos){
			w+=" ionizE ionizS bremsE bremsS photoE photoS epairE epairP epairS";
                        if(cros->get_photonuclear()->get_form()>2){
				w+=" photoP";
			}
		}
		if(w.find("contiR")!=string::npos){
			w+=" gaussR en2ldX en2ldE";
		}

		w=toLowerCase(w);
//cout<<"interpolate line 530\t"<<w<<endl;
//cout<<"ionization cont: 1 pointer"<<endl;
      try{
                cros->get_ionization()->get_Continuous()->set_jt(false);
                if(w.find("ionize")!=string::npos && cros->get_ci()>0){
                        cout<<"Parameterizing ionizE ... ";
                        cros->get_ionization()->get_Continuous()->interpolateJ_ = new Interpolate(NUM1, e_low, e_hi, cros->get_ionization()->get_Continuous(), g, true, false, true, g, false, false, true);
                        cros->get_ionization()->get_Continuous()->set_jt(true);
                        cout<<"done \n";
//cout<<cros->get_ionization()->get_Continuous()->interpolateJ_ <<endl;
		}
//cout<<"ionization stoch: 2 pointers"<<endl;
                cros->get_ionization()->get_Stochastic()->set_jt(false);
                if(w.find("ionizs")!=string::npos && cros->get_ci()>0){
                        cout<<"Parameterizing ionizS ... ";
                        cros->get_ionization()->get_Stochastic()->interpolateJ_ = new Interpolate(NUM1, e_low, e_hi, NUM1, 0, 1, cros->get_ionization()->get_Stochastic(), g, false, false, true, g, false, false, false, g, true, false, false);
                        cros->get_ionization()->get_Stochastic()->interpolateJo_  = new Interpolate(NUM1, e_low, e_hi, cros->get_ionization()->get_Stochastic(), g, false, false, true, g, true, false, false);
                        cros->get_ionization()->get_Stochastic()->set_jt(true);
                        cout<<"done "<<endl;

		}
//cout<<"bremsstrahlung cont: 1 pointer"<<endl;
                cros->get_bremsstrahlung()->get_Continuous()->set_jt(false);
                if(w.find("bremse")!=string::npos && cros->get_cb()>0){
                        cout<<"Parameterizing bremsE ... ";
                        cros->get_bremsstrahlung()->get_Continuous()->interpolateJ_ = new Interpolate(NUM1, e_low, e_hi, cros->get_bremsstrahlung()->get_Continuous(), g, true, false, true, g, false, false, false);
                        cros->get_bremsstrahlung()->get_Continuous()->set_jt(true);
                        cout<<"done \n";
//cout<<cros->get_bremsstrahlung()->get_Continuous()->interpolateJ_<<endl;
		}
//cout<<"bremsstrahlung stoch: 2*"<<medium_->get_numCompontents()<<" pointers"<<endl;
                cros->get_bremsstrahlung()->get_Stochastic()->set_jt(false);
                if(w.find("bremss")!=string::npos && cros->get_cb()>0){
                        cout<<"Parameterizing bremsS ... ";
                        cros->get_bremsstrahlung()->get_Stochastic()->interpolateJ_ = (Interpolate *)calloc(medium_->get_numCompontents(),sizeof(Interpolate));
                        cros->get_bremsstrahlung()->get_Stochastic()->interpolateJo_ = (Interpolate *)calloc(medium_->get_numCompontents(),sizeof(Interpolate));

                        for(i=0; i<medium_->get_numCompontents(); i++){
                                cros->set_component(i);
                                cros->get_bremsstrahlung()->get_Stochastic()->interpolateJ_[i] = Interpolate(NUM1, e_low, e_hi, NUM1, 0, 1, cros->get_bremsstrahlung()->get_Stochastic(), g, false, false, true, g, false, false, false, g, true, false, false);
                                cros->get_bremsstrahlung()->get_Stochastic()->interpolateJo_[i] = Interpolate(NUM1, e_low, e_hi, cros->get_bremsstrahlung()->get_Stochastic(), g, false, false, true, g, true, false, false);

                        }

                        cros->get_bremsstrahlung()->get_Stochastic()->set_jt(true);
                        cout<<"done \n";
		}

                cros->get_photonuclear()->set_jt(false);
                if(w.find("photop")!=string::npos && cros->get_cp()>0){
                        cout<<"Parameterizing photoP ... ";
                        cros->get_photonuclear()->interpolateJ_ = (Interpolate *)calloc(medium_->get_numCompontents(),sizeof(Interpolate));

                        for(i=0; i<medium_->get_numCompontents(); i++){
                            cros->set_component(i);
                            cros->get_photonuclear()->interpolateJ_[i] = Interpolate(NUM1, e_low, e_hi, NUM1, 0., 1., cros->get_photonuclear(), g, false, false, true, g, false, false, false, g, false, false, false);

                        }

                        cros->get_photonuclear()->set_jt(true);
                        cout<<"done \n";
                }

                cros->get_photonuclear()->get_Continuous()->set_jt(false);
                if(w.find("photoe")!=string::npos && cros->get_cp()>0){
                        cout<<"Parameterizing photoE ... ";
                        cros->get_photonuclear()->get_Continuous()->interpolateJ_ = new Interpolate(NUM1, e_low, e_hi, cros->get_photonuclear()->get_Continuous(), g, true, false, true, g, false, false, false);
                        cros->get_photonuclear()->get_Continuous()->set_jt(true);
                        cout<<"done \n";
//cout<<cros->get_photonuclear()->get_Continuous()->interpolateJ_<<endl;
		}
//cout<<"photonuclear stoch: 2*"<<medium_->get_numCompontents()<<" pointers"<<endl;
                cros->get_photonuclear()->get_Stochastic()->set_jt(false);
                if(w.find("photos")!=string::npos && cros->get_cp()>0){
                        cout<<"Parameterizing photoS ... ";
                        cros->get_photonuclear()->get_Stochastic()->interpolateJ_ = (Interpolate *)calloc(medium_->get_numCompontents(),sizeof(Interpolate));
                        cros->get_photonuclear()->get_Stochastic()->interpolateJo_ = (Interpolate *)calloc(medium_->get_numCompontents(),sizeof(Interpolate));

                        for(i=0; i<medium_->get_numCompontents(); i++){
                                cros->set_component(i);
                                cros->get_photonuclear()->get_Stochastic()->interpolateJ_[i] = Interpolate(NUM1, e_low, e_hi, NUM1, 0, 1, cros->get_photonuclear()->get_Stochastic(), g, false, false, true, g, false, false, false, g, true, false, false);
                                cros->get_photonuclear()->get_Stochastic()->interpolateJo_[i] = Interpolate(NUM1, e_low, e_hi, cros->get_photonuclear()->get_Stochastic(), g, false, false, true, g, true, false, false);
//cout<<cros->get_photonuclear()->get_Stochastic()->interpolateJ_.at(i)<<endl;
//cout<<cros->get_photonuclear()->get_Stochastic()->interpolateJo_.at(i)<<endl;
                        }

                        cros->get_photonuclear()->get_Stochastic()->set_jt(true);
                        cout<<"done \n";
		}
//cout<<"epairproduction: "<<medium_->get_numCompontents()<<" pointers"<<endl;
                cros->get_epairproduction()->set_jt(false);
                if(w.find("epairp")!=string::npos && cros->get_ce()>0){
                        cout<<"Parameterizing epairP ... ";
                        cros->get_epairproduction()->interpolateJ_ = (Interpolate *)calloc(medium_->get_numCompontents(),sizeof(Interpolate));

                        for(i=0; i<medium_->get_numCompontents(); i++){
                                cros->set_component(i);
                                cros->get_epairproduction()->interpolateJ_[i] = Interpolate(NUM1, e_low, e_hi, NUM1, 0., 1., cros->get_epairproduction(), g, false, false, true, g, false, false, false, g, false, false, false);
//cout<<cros->get_epairproduction()->interpolateJ_.at(i)<<endl;
                        }

                        cros->get_epairproduction()->set_jt(true);
                        cout<<"done \n";
		}

                cros->get_epairproduction()->get_Continuous()->set_jt(false);
                if(w.find("epaire")!=string::npos && cros->get_ce()>0){
                        cout<<"Parameterizing epairE ... ";
                        cros->get_epairproduction()->get_Continuous()->interpolateJ_ =new Interpolate(NUM1, e_low, e_hi, cros->get_epairproduction()->continuous_, g, true, false, true, g, false, false, false);
                        cros->get_epairproduction()->get_Continuous()->set_jt(true);
                        cout<<"done \n";

                }

                cros->get_epairproduction()->get_Stochastic()->set_jt(false);
                if(w.find("epairs")!=string::npos && cros->get_ce()>0){
                        cout<<"Parameterizing epairS ... ";
                        cros->get_epairproduction()->get_Stochastic()->interpolateJ_= (Interpolate *)calloc(medium_->get_numCompontents(),sizeof(Interpolate));
                        cros->get_epairproduction()->get_Stochastic()->interpolateJo_= (Interpolate *)calloc(medium_->get_numCompontents(),sizeof(Interpolate));
                        for(i=0; i<medium_->get_numCompontents(); i++){
                                cros->set_component(i);
                                cros->get_epairproduction()->get_Stochastic()->interpolateJ_[i] = Interpolate(NUM1, e_low, e_hi, NUM1, 0, 1, cros->get_epairproduction()->stochastic_, g, false, false, true, g, false, false, false, g, true, false, false);
                                cros->get_epairproduction()->get_Stochastic()->interpolateJo_[i] = Interpolate(NUM1, e_low, e_hi, cros->get_epairproduction()->stochastic_, g, false, false, true, g, true, false, false);

                        }

                        cros->get_epairproduction()->get_Stochastic()->set_jt(true);
                        cout<<"done \n";
		}

                this->jt=false;
		if(w.find("tracke")!=string::npos){
                        cout<<"Parameterizing trackE ... ";

                        this->interpolateJ_ = (Interpolate *)calloc(2,sizeof(Interpolate));
                        this->interpolateJdf_ = (Interpolate *)calloc(2,sizeof(Interpolate));
			pint=true;
                        if(abs(-integral_.at(1)->integrateWithLog(particle_->low, particle_->low*10, this))<abs(-integral_.at(1)->integrateWithLog(PhysicsModel::get_ebig(), PhysicsModel::get_ebig()/10, this))){
				up=true;
			}
			else{
				up=false;
			}
			for(i=0; i<2; i++){
                                pint=(i==1);
                                this->df=false;
                                this->interpolateJ_[i] = Interpolate(NUM3, e_low, e_hi, this, g, false, false, true, g, false, false, false);
                                this->df=true;
                                this->interpolateJdf_[i] = Interpolate(NUM3, e_low, e_hi, this, g, false, false, true, g, false, false, false);

			}

                        this->jt=true;
			for(i=0; i<2; i++){
				if(!(up&&(i==1))){
                                        bigLow[i]=interpolateJ_[i].interpolate(particle_->low);
				}
			}
			if(up){
                                cout<<"done \n";
			}
			else{
                                cout<<"down \n";
			}
		}

                cros->set_jt(false);
		if(w.find("trackx")!=string::npos){
                        cout<<"Parameterizing trackX ... ";
                        cros->set_df(false);
                        cros->interpolateJ_ = new Interpolate(NUM3, e_low, e_hi, cros, g, false, false, true, g, false, false, false);
                        cros->set_df(true);
                        cros->interpolateJdf_ = new Interpolate(NUM3, e_low, e_hi, cros, g, false, false, true, g, false, false, false);
                        cros->set_jt(true);
                        cout<<"done \n";

                }

                particle_->jt=false;
		if(w.find("trackt")!=string::npos){
                        cout<<"Parameterizing trackT ... ";
                        particle_->df=false;
                        particle_->interpolateJ_ = new Interpolate(NUM3, e_low, e_hi, particle_, g, false, false, true, g, false, false, false);
                        particle_->df=true;
                        particle_->interpolateJdf_ = new Interpolate(NUM3, e_low, e_hi, particle_, g, false, false, true, g, false, false, false);
                        particle_->jt=true;
                        cout<<"done \n";
		}

                particle_->get_scattering()->set_jt(false);
		if(w.find("molies")!=string::npos){
                        cout<<"Parameterizing molieS ... ";
                        particle_->get_scattering()->set_df(false);
                        particle_->get_scattering()->interpolateJ_ =  new Interpolate(NUM2, e_low, e_hi, particle_->get_scattering(), g, false, false, true, g, false, false, false);
                        particle_->get_scattering()->set_df(true);
                        particle_->get_scattering()->interpolateJdf_ = new Interpolate(NUM2, e_low, e_hi, particle_->get_scattering(), g, false, false, true, g, false, false, false);
                        particle_->get_scattering()->set_jt(true);
                        cout<<"done"<<endl;
;
                }
;
                StandardN->set_jt(false);
		if(w.find("gaussr")!=string::npos){
                        cout<<"Parameterizing gaussR ... ";
                        StandardN->interpolateJ_ = new Interpolate(NUM2, -5, 5, StandardN, g, true, false, false, g, true, false, false);
                        StandardN->set_jt(true);
                        cout<<"done \n";
;
                }

                E2Loss_->e2lx->set_jt(false);
		if(w.find("en2ldx")!=string::npos){
                        cout<<"Parameterizing en2ldX ... ";
                        E2Loss_->e2lx->interpolateJ_ = new Interpolate(NUM2, e_low, e_hi, E2Loss_->e2lx, g, false, false, true, g, false, false, false);
                        E2Loss_->e2lx->set_jt(true);
                        cout<<"done \n";
                }

                E2Loss_->e2le->set_jt(false);
		if(w.find("en2lde")!=string::npos){
                        cout<<"Parameterizing en2ldE ... ";
                        E2Loss_->e2le->set_df(false);
                        E2Loss_->e2le->interpolateJ_ = new Interpolate(NUM2, e_low, e_hi, E2Loss_->e2le, g, false, false, true, g, false, false, false);
                        E2Loss_->e2le->set_df(true);
                        E2Loss_->e2le->interpolateJdf_ = new Interpolate(NUM2, e_low, e_hi, E2Loss_->e2le, g, false, false, true, g, false, false, false);
                        E2Loss_->e2le->set_jt(true);
                        cout<<"done \n";
                }

                cout<<"Finished parameterizations"<<endl;
                }catch(int a){
                    throw 0;
                }

    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * call this routine to enable interpolations and their save and reread. To enable everything, set w="all"
     */

	void Propagate::interpolate(string w, string filename){
                stringstream name;
		bool flag;
                name<<filename<<".";
                name<<particle_->name<<"_"<<replaceAll((toLowerCase((medium_->get_name()+"_"+w))),' ', '-');
                name<<"_"<<Output::f(medium_->get_ecut())<<"_"<<Output::f(medium_->get_vcut())<<"_";
                name<<cros->get_bremsstrahlung()->get_form();
                name<<cros->get_photonuclear()->get_form();
                name<<cros->get_photonuclear()->get_bb();
                name<<cros->get_photonuclear()->get_shadow();
                if(cros->get_lpm()){
                        name<<"t";
		}
		else{
                        name<<"f";
		}
		if(exactTime){
                        name<<"t";
		}
		else{
                        name<<"f";
		}
		if(contiCorr){
                        name<<"t";
		}
		else{
                        name<<"f";
		}
		if(molieScat){
                        name<<"t";
		}
		else{
                        name<<"f";
		}
                if(PhysicsModel::get_elow()==particle_->low){
                        name<<"_l"<<Output::f(PhysicsModel::get_elow());
		}
                if(PhysicsModel::get_ebig()!=BIGENERGY){
                        name<<"_b"<<Output::f(PhysicsModel::get_ebig());
		}
                if(cros->get_ci()!=1 || cros->get_cb()!=1 || cros->get_cp()!=1 || cros->get_ce()!=1 || cros->get_cd()!=1 || medium_->get_rho()!=1){
                        name<<"_"<<Output::f(cros->get_ci());
                        name<<","<<Output::f(cros->get_cb());
                        name<<","<<Output::f(cros->get_cp());
                        name<<","<<Output::f(cros->get_ce());
                        name<<","<<Output::f(cros->get_cd());
                        name<<","<<Output::f(medium_->get_rho());
		}
                if(Output::raw){
                        name<<"_raw";
		}
		else{
                        name<<"_ascii";
		}
                name<<".data";

		do {
//                    if(Output::texi){
//                            return;
//                    }
//                    flag=false;
//                    try{
//                        Output::open(name.str());
//                        interpolate(w);
//                        Output::close();
//                    }catch (int mmcException){
//                        flag=true;
//                        Output::Delete(name.str());
//                    }
                    if(Output::texi){
                            return;
                    }

                    flag=false;
                        try{
                        Output::open(name.str());
                        interpolate(w);
                        Output::close();
                    } catch (int a){
                        cout<<"EXCEPTION"<<endl;
                        flag=true;
                        Output::Delete(name.str());
                    }
		} while(flag);
    }

    //----------------------------------------------------------------------------------------------------//

	/**
     * returns the value of the tracking integral
     */

	double Propagate::getpr(double ei, double rnd, bool pint){
//cout<<"Propagate::getpr starts:\t jt="<<jt<<"\t pint="<<pint<<endl;
		if(jt){
			if(pint){
                                storeDif[1]=interpolateJ_[1].interpolate(ei);
			}
			else{
                                storeDif[0]=interpolateJ_[0].interpolate(ei);
			}


			if(up&&pint){
				if(pint){
                                        return max(storeDif[1], 0.0);
				}
				else{
                                        return max(storeDif[0], 0.0);
				}
                        }
			else{
				if(pint){
                                        return max(bigLow[1]-storeDif[1], 0.0);
				}
				else{
                                        return max(bigLow[0]-storeDif[0], 0.0);
				}
			}
		}
                else{

//cout<<"INTEGRAL"<<endl;
                        this->pint=pint;
			if(pint){
                                return integral_.at(1)->integrateWithLog(ei, particle_->low, this, -rnd);
			}
			else{
                                return integral_.at(0)->integrateWithLog(ei, particle_->low, this, -rnd);
			}
		}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * final energy, corresponding to the given value rnd of the tracking integral                cout <<"vmax: "<<vMax_<<endl;
     */

	double Propagate::getef(double ei, double rnd, bool pint){
		if(jt){
			if(pint){
                                if(abs(rnd)>abs(storeDif[1])*HALF_PRECISION){
					double aux;
					if(up&&pint){
						if(pint){
                                                        aux=interpolateJ_[1].findLimit(storeDif[1]-rnd);
						}
						else{
                                                        aux=interpolateJ_[0].findLimit(storeDif[0]-rnd);
						}
					}
					else{
						if(pint){
                                                        aux=interpolateJ_[1].findLimit(storeDif[1]+rnd);
						}
						else{
                                                        aux=interpolateJ_[0].findLimit(storeDif[0]+rnd);
						}
					}
                                        if(abs(ei-aux)>abs(ei)*HALF_PRECISION){
                                                return min(max(aux, particle_->low), ei);
					}
				}
			}
			else{
                                if(abs(rnd)>abs(storeDif[0])*HALF_PRECISION){
					double aux;
					if(up&&pint){
						if(pint){
                                                        aux=interpolateJ_[1].findLimit(storeDif[1]-rnd);
						}
						else{
                                                        aux=interpolateJ_[0].findLimit(storeDif[0]-rnd);
						}
					}
					else{
						if(pint){
                                                        aux=interpolateJ_[1].findLimit(storeDif[1]+rnd);
						}
						else{
                                                        aux=interpolateJ_[0].findLimit(storeDif[0]+rnd);
						}
					}
                                        if(abs(ei-aux)>abs(ei)*HALF_PRECISION){
                                                return min(max(aux, particle_->low), ei);
					}
				}
			}
			if(pint){
                                return min(max(ei+rnd/interpolateJdf_[1].interpolate(ei+rnd/(2*interpolateJdf_[1].interpolate(ei))), particle_->low), ei);
			}
			else{
                                return min(max(ei+rnd/interpolateJdf_[0].interpolate(ei+rnd/(2*interpolateJdf_[0].interpolate(ei))), particle_->low), ei);
			}
		}
		else{
                        this->pint=pint;
			if(pint){
                                return integral_.at(1)->getUpperLimit();
			}
			else{
                                return integral_.at(0)->getUpperLimit();
			}

                }
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * 1d parametrization - interface to Interpolate
     */

	double Propagate::functionInt(double e){
		if(df){
			return function(e);
		}
		else{
			if(up&&pint){
				if(pint){
                                        return integral_.at(1)->integrateWithLog(e, particle_->low, this);
				}
				else{
                                        return integral_.at(0)->integrateWithLog(e, particle_->low, this);
				}
			}
			else{
				if(pint){
                                        return -integral_.at(1)->integrateWithLog(e, PhysicsModel::get_ebig(), this);
				}
				else{
                                        return -integral_.at(0)->integrateWithLog(e, PhysicsModel::get_ebig(), this);
				}
			}
		}
    }

// Getter
        StandardNormal* Propagate::get_Standard() { return StandardN; };
        CrossSections* Propagate::get_cros() { return cros;}
        Medium* Propagate::get_Medium() { return medium_;}
//Setter
        void Propagate::set_exactTime(bool eTime){exactTime=eTime;}
        void Propagate::set_molieScat(bool mScat){molieScat=mScat;}
        void Propagate::set_contiCorr(bool cCorr){contiCorr=cCorr;};
        void Propagate::set_pint(bool newPint){pint=newPint;}
        void Propagate::set_rho(double newRho){rho=newRho;}
