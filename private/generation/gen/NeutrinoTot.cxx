#include "generation/gen/NeutrinoTot.h"
#include <algorithm>
#include <cmath>
#include "PROPOSAL/Output.h"
#include <sstream>
using namespace std;

/**
 * Class contains functions for calculation of neutrino interaction cross sections.
 */


    //----------------------------------------------------------------------------------------------------//

    /**
     * Class constructor.
     */

    NeutrinoTot::NeutrinoTot(){
        jt=false;
        N = new NeutrinoInt();
        I.resize(4);
        for(int i=0; i<4; i++) I.at(i) = new Integral(IROMB, IMAXS, IPREC);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Charged and Neutral current interaction for neutrinos and antineutrinos
     */

    double NeutrinoTot::dSdy(double E, bool cc, bool nu){
        this->E=E;
        if(jt) return max(Jo.at((cc?1:0)+(nu?2:0))->interpolate(E), 0.);
        this->cc=cc;
        this->nu=nu;
        return I.at((cc?1:0)+(nu?2:0))->integrateOpened(0, 1, this);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Charged and Neutral current interaction for neutrinos and antineutrinos
     */

    double NeutrinoTot::dSdy(double E, bool cc, bool nu, double rnd){
        this->E=E;
        if(jt){
            rns=rnd;
            int i=(cc?1:0)+(nu?2:0);
            return H[i]=max(Jo.at(i)->interpolate(E), 0.);
        }
        this->cc=cc;
        this->nu=nu;
        return I.at((cc?1:0)+(nu?2:0))->integrateOpened(0, 1, this, rnd);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Returns the energy transferred to the hardonic state
     */

    double NeutrinoTot::e(bool cc, bool nu){
        if(jt){
            int i=(cc?1:0)+(nu?2:0);
            return E*J.at(i)->findLimit(E, rns*H[i]);
        }
        return E*I.at((cc?1:0)+(nu?2:0))->getUpperLimit();
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Glashow resonance cross section
     */

    double NeutrinoTot::dSdy(double E){
        return dSdy(E, 0);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Glashow resonance cross section
     */

    double NeutrinoTot::dSdy(double E, double rnd){
        this->E=E;
        this->rnd=rnd;
        double aux=1.e-6*MW*MW;
        double s=2.e-3*ME*E;
        aux*=aux/((s-aux)*(s-aux)+1.e-6*GW*GW*aux);
        aux*=9.e12*GF*GF*s/(3*PI);
        s=RE*1.e-3*ME/ALPHA;
        return N->RZ*s*s*aux;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Glashow resonance energy lost by neutrino
     */

    double NeutrinoTot::e(double m){
        double M=sqrt(2.e-3*ME*E+1.e-6*ME*ME);
        double e=(M*M-1.e-6*m*m)/(2*M);
        return E-(e*((E+1.e-3*ME)/M)+e*(E/M)*(2*rnd-1));
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Neutrino oscillation probability (mu->tau); units: E [GeV] L [m]
     */

    double NeutrinoTot::Pmt(double E, double L){
        double aux;
        aux=DE2*L/E;
        aux*=(1.e-12*1.e2/1.e3)*(ALPHA/(RE*ME))/4;
        aux=sin(aux);
        aux*=aux;
        return ST2*aux;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * returns the neutrino cross section dSdy - interface to Integral.
     */

    double NeutrinoTot::function(double y){
        return N->RA*N->dSdy(y, E, cc, nu);
        //return N->dSdy(y, E, cc, nu);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Prints the cross sections in [cm^2] as a function of energy in [GeV].
     */

//    static void NeutrinoTot::main(String[] args){
//        NeutrinoTot T = new NeutrinoTot();
//        Output::raw=true;
//        T.interpolate(".atmflux");
//        for(double E=1.e1; E<1.e12; E*=1.05){
//            cout<<Output::f(E)<<" "<<Output::f(T.dSdy(E, true, true)/T.N.RA)<<" "<<Output::f(T.dSdy(E, true, false)/T.N.RA)<<" "<<Output::f(T.dSdy(E, false, true)/T.N.RA)<<" "<<Output::f(T.dSdy(E, false, false)/T.N.RA)<<" "<<Output::f(T.dSdy(E)/T.N.RZ))<<endl;
//        }
//    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Parametrizes the integral of this class.
     */

    void NeutrinoTot::interpolate(string name){
        bool flag;
        stringstream ss;
        ss<<name;
        ss<<".gen_neutrino";

        if(nlow_!=ME) ss<<"_l"<<Output::f(nlow_);
        if(ebig_!=BIGENERGY) ss<<"_b"<<Output::f(ebig_);

        if(Output::raw) ss<<"_raw"; else ss<<"_ascii";
        ss<<".data";
        name=ss.str();
        do {
            if(Output::texi) return;
            flag=false;
            try{
                if(!Output::FileExist(name)){
                    CteqPDF::Ctq->interpolate();
                    N->interpolate();
                }
                Output::open(name);

                interpolate();
                cerr<<"Finished parameterizations"<<endl;

                Output::close();
            }catch(int a){
                throw 0;
                flag=true;
                Output::Delete(name);
            }
        } while(flag);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Parametrizes the integral of this class.
     */

    void NeutrinoTot::interpolate(){
        int g=5;
        double e_hi=ebig_*1.e-3;
        double e_low=nlow_*1.e-3;

        jt=false;
        H = new double[4];
        cerr<<"Parameterizing neutrT ... "<<endl;
        J.resize(4);
        Jo.resize(4);
        for(int i=0; i<4; i++){
            switch(i){
            case 0: cc=false; nu=false; break;
            case 1: cc=true;  nu=false; break;
            case 2: cc=false; nu=true;  break;
            case 3: cc=true;  nu=true;  break;
            }
            J.at(i) = new Interpolate(NUM1, e_low, e_hi, NUM1, 0, 1, this, g, false, false, true, g, false, false, false, g, true, false, false);
            Jo.at(i) = new Interpolate(NUM1, e_low, e_hi, this, g, false, false, true, g, true, false, false);
        }
        jt=true;
        cerr<<"done"<<endl;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * 2d parametrization - interface to Interpolate
     */

    double NeutrinoTot::functionInt(double E, double y){
        this->E=E;
        return I.at((cc?1:0)+(nu?2:0))->integrateOpened(0, y, this);
    }

    //----------------------------------------------------------------------------------------------------//


    /**
     * 1d parametrization - interface to Interpolate
     */

    double NeutrinoTot::functionInt(double E){
        return dSdy(E, cc, nu);
    }



