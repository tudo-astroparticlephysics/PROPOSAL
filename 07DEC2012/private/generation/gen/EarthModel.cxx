#include "generation/gen/EarthModel.h"
#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;
/**
 * Earth density model implementation.
 */



    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize class with atmospheric model, ground elevation z0 in [km], bedrock elevation b0 in [km], and detector depth D.
     */

    EarthModel::EarthModel(int model, double z0, double b0, double D){

        z0=0;
        b0=0;
        D=0;
        Rs=R0;
        ts=0.01;
        num=0;
        jt=false;
        this->z0=z0;
        this->b0=b0;
        this->D=D;
        Rs=R0+z0-D;
        A = new Atmosphere(model, 0);
        I = new Integral(IROMB, IMAXS, IPREC);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates density in [g/cm^3] as a function of distance to the center R in [km] according to the preliminary Earth model.
     */

    double EarthModel::rho(double R){
        double aux=0, x;
        x=R/R0;
        if (R<0) aux=0;
        else if (R<1221.5) aux=13.0885-8.8381*x*x;
        else if (R<3480.0) aux=12.5815-x*(1.2638+x*(3.6426+x*5.5281));
        else if (R<5701.0) aux=7.9565-x*(6.4761-x*(5.5283-x*3.0807));
        else if (R<5771.0) aux=5.3197-1.4836*x;
        else if (R<5971.0) aux=11.2494-8.0298*x;
        else if (R<6151.0) aux=7.1089-3.8045*x;
        else if (R<6346.6) aux=2.691+0.6924*x;
        else if (R<6356.0) aux=2.9;
        else if (R<=R0+b0) aux=2.6;
        else{ for(int i=num-1; i>=0; i--) if((mRo[i]>0 && R<mRa[i]) || i==0){ aux=mRo[i]; break; } }
        if(R>R0+z0) if(mRa[0]<0) aux=A->dXdh(R-R0)*1.e-5;  // not used by MMC
        return aux;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates distance to the center r in [km] as a function of path length in [km]. sett must be called first.
     */

    double EarthModel::r(double t){
        return sqrt(t*t+2*pz*Rs*t+Rs*Rs);
    }

    //----------------------------------------------------------------------------------------------------//


    /**
     * Calculates and stores path length in [km] to the entry and exit points.
     */

    void EarthModel::sett(double px, double py, double pz){
        this->px=px;
        this->py=py;
        this->pz=pz;
        double r=R0+z0;
        double aux=sqrt(r*r-Rs*Rs*(px*px+py*py));
        ti=-pz*Rs-aux;
        tf=-pz*Rs+aux;
        if(jt) Xmid=J->interpolate(abs(pz), 1.);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates path length in [km] from the entry point. sett must be called first.
     */

    double EarthModel::r(double x, double y, double z){
        return -pz*Rs+(px*x+py*y+pz*(Rs+z));
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates mass overburden in [g/cm^2] as a function of path length in [km] to the entry and exit points.
     */

    double EarthModel::X(double t1, double t2){
        if(jt) return jInt(t2)-jInt(t1);
        double tm1=min(t2, max(t1, -100.)), tm2=max(t1, min(t2, 100.));
        return 1.e5*(I->integrateOpened(t1, tm1, this)+I->integrateOpened(tm1, tm2, this)+I->integrateOpened(tm2, t2, this));
    }

    //----------------------------------------------------------------------------------------------------//


    /**
     * Calculates mass overburden in [g/cm^2] as a function of path length in [km] to the entry and exit points.
     * Also calculates path length corresponding to the mass overburden X0 in [g/cm^2] from the entry point.
     */

    double EarthModel::X(double t1, double t2, double X0){
        if(jt){
            this->t1=t1; this->t2=t2;
            double X1=jInt(t1);
            this->X0=X0+X1;
            return jInt(t2)-X1;
        }
        return 1.e5*I->integrateOpened(t1, t2, this, -X0*1.e-5);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Returns the path length in [km] corresponding to the mass overburden X0 in [g/cm^2] set by a call to X(t1, t2, X0).
     */

    double EarthModel::t(){
        if(jt) return min(max(fInt(X0), t1), t2);
        else return I->getUpperLimit();
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Returns the density in [g/cm^3] as a function of the path length in [km] - interface to Integral.
     */

    double EarthModel::function(double x){
        return rho(r(x));
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Returns the density in [g/cm^3] as a function of the path length in [km], gives density of the internal medium.
     */

    double EarthModel::intRho(double x){
        return rho(min(r(x), R0+b0));
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Prints the density profile.
     */

//    static void EarthModel::main(string[] args){
//        EarthModel E = new EarthModel(0, 2.834, 0.024, 1.730);
//        E.num=1;
//        E.mRa = new double[1];
//        E.mRo = new double[1];
//        E.mRa[0]=E.R0+E.z0;
//        E.mRo[0]=0.917;
//        Output.out.println("Radius [km], density [g/cm^3]:");
//        for(double R=0; R<7.e3; R+=10){
//            Output.out.println(Output.f(R)+" "+Output.f(E.rho(R)));
//        }
//        Output.out.println("\nzenith angle, mass overburden [g/cm^2]:");
//        for(double th=0; th<=180; th+=1){
//            double pz=Math.cos(th*Pi/180);
//            E.sett(Math.sqrt(1-pz*pz), 0, pz);
//            Output.out.println(Output.f(th)+" "+Output.f(E.X(E.ti, 0)));
//        }
//    }

    //----------------------------------------------------------------------------------------------------//

    double EarthModel::jInt(double t){
        double x;
        x=(t-ti)/(tf-ti);
        if(x<0.5) return J->interpolate(abs(pz), (2-2*ts)*x+ts);
        else return 2*Xmid-J->interpolate(abs(pz), (2-2*ts)*(1-x)+ts);
    }

    //----------------------------------------------------------------------------------------------------//

    double EarthModel::fInt(double X){
        double x;
        if(X<Xmid) x=(J->findLimit(abs(pz), X)-ts)/(2-2*ts);
        else x=1-(J->findLimit(abs(pz), 2*Xmid-X)-ts)/(2-2*ts);
        return ti+(tf-ti)*x;
    }

    //----------------------------------------------------------------------------------------------------//


    /**
     * Parametrizes the integral of this class.
     */

    void EarthModel::interpolate(){
        int g=3, g2=2;
        jt=false;
        cerr<<"Parameterizing earthM ... "<<endl;
        J = new Interpolate(NUM1, 0, 1, NUM3, ts, 1, this, g, true, false, false, g2, true, false, true, g2, false, false, false);
        jt=true;
        cerr<<"done"<<endl;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * 2d parametrization - interface to Interpolate
     */

    double EarthModel::functionInt(double x, double t){
        sett(sqrt(1-x*x), 0, x);
        return X(ti, ti+(tf-ti)*(t-ts)/(2-2*ts));
    }


