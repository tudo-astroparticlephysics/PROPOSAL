#include "generation/gen/NeutrinoInt.h"
#include <algorithm>
#include <cmath>
#include <iostream>

using namespace std;

/**
 * Class contains functions for calculation of neutrino interaction differential cross sections.
 */


    //----------------------------------------------------------------------------------------------------//

    /**
     * Class constructor. Sets up the values of cross section parameters and calls the CteqPDF constructor.
     */

    NeutrinoInt::NeutrinoInt(){

        name="Standard Rock";
        jt=false;

        double d12, d22, d32, d42;
        d12=1/2.-2*XW/3; d12*=d12;
        d22=-1/2.+XW/3;  d22*=d22;
        d32=-2*XW/3;     d32*=d32;
        d42=XW/3;        d42*=d42;
        d1=d12+d22+d32+d42;
        d2=d12+d32-d22-d42;
        d3=d12+d22-d32-d42;
        F = new CteqPDF();
        m = new Medium(name, -1, 1, 1);
        RA=m->get_molDensity()*m->get_sumNucleons()/m->get_massDensity();
        RZ=m->get_molDensity()*m->get_sumCharge()/m->get_massDensity();
        I = new Integral(IROMB, IMAXS, IPREC);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * returns the neutrino (nu=true) and anti neutrino (nu=false) charged (cc=true) and neutral (cc=false) current
     * interaction cross sections as functions of q(x,y) and y and energy E in [GeV].
     */

    double NeutrinoInt::dS2dxdy(double q, double y, double E, bool cc, bool nu){
        int inu;
        inu=nu?1:-1;

        double aux, sum, x, Q2, qa, F2, F3;
        Q2=q*q;
        qa=1.e-3*m->get_MM()*y*E;
        x=q*q/(2*qa);
        if(x>1) return 0;

        if(cc){
            F2=F->PDF(2, x, q);
            F3=F->PDF(1, x, q)+inu*F->PDF(3, x, q);
        }
        else{
            F2=F->PDF(2, x, q)*d1-F->PDF(3, x, q)*d2;
            F3=F->PDF(1, x, q)*d3;
        }

        aux=1.e-3*(cc?MW:MZ);
        aux*=aux;
        sum=aux/(Q2+aux);
        sum*=sum*(1.e9*GF*GF*m->get_MM()*E/PI);
        aux=y-y*y/2;
        sum*=(1-aux)*F2+inu*aux*F3;
        aux=1.e-3*RE*ME/ALPHA;
        aux*=aux;

        aux*=sum*q/qa;
        return max(aux, 0.);  // multiply by (m.No*m.totA) to get 1/length
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * returns the neutrino cross sections integrated over x as functions of y and energy E in [GeV].
     */

    double NeutrinoInt::dSdy(double y, double E, bool cc, bool nu){
        if(jt) return max(J.at((cc?1:0)+(nu?2:0))->interpolate(E, y), 0.);
        this->y=y;
        this->E=E;
        this->cc=cc;
        this->nu=nu;
        double qmax=sqrt(2.e-3*m->get_MM()*y*E);
        return I->integrateWithLog(min(F->Qn, qmax), min(F->Qm, qmax), this);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * returns the neutrino cross section dS2dxdy - interface to Integral.
     */

    double NeutrinoInt::function(double x){
        return dS2dxdy(x, y, E, cc, nu);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Parametrizes the integral of this class.
     */

    void NeutrinoInt::interpolate(){
        int g=5;
        double e_hi=ebig_*1.e-3;
        double e_low=nlow_*1.e-3;

        jt=false;
        cerr<<"Parameterizing neutrI ... "<<endl;
        J.resize(4);
        for(int i=0; i<4; i++){
            switch(i){
            case 0: cc=false; nu=false; break;
            case 1: cc=true;  nu=false; break;
            case 2: cc=false; nu=true;  break;
            case 3: cc=true;  nu=true;  break;
            }
            J.at(i) = new Interpolate(NUM1, e_low, e_hi, NUM1, 0, 1, this, g, false, false, true, g, false, false, false, g, true, false, false);
        }
        jt=true;
        cerr<<"done"<<endl;
    }

    //----------------------------------------------------------------------------------------------------//


    /**
     * 2d parametrization - interface to Interpolate
     */

    double NeutrinoInt::functionInt(double E, double y){
        return dSdy(y, E, cc, nu);
    }


