#include "generation/gen/CosCorr.h"
#include <cmath>


/**
 * This class calculates atmosphere-dependent quantities for a fixed first interaction depth.
 */


    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize class with atmospheric model and ground elevation z0 in [km].
     * h0 is the average production height in km or in g/cm^2 if negative.
     */

    CosCorr::CosCorr(int model, double h0, double z0){

        st=1;
        h0=19;
        X0=85;

        R=EarthModel::R0;
        A = new Atmosphere(model, z0);
        I = new Integral(IROMB, IMAXS, IPREC);
        if(h0<0){
            this->X0=-h0;
            this->h0=A->h(X0);
        }
        else{
            this->h0=h0;
            this->X0=A->X(h0);
        }
        R+=z0;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Prints zenith angle profile of production height, cos*, total X, track length.
     */

//    static void CosCorr::main(String[] args){
//        int m=1;
//        double aux, x, z0=0, h0=15.5;

//        for(int n=0; n<args.length; n++){
//            if(args[n].equals("-help") || args[n].equals("-h") || args[n].equals("--help")){
//                Output.out.println("\n"+
//"This class calculates muon production height, cos(th*), etc.\n"+
//"                       -m=[0-3] choose atmosphere model\n"+
//"                       -z0=[ground elevation in km]\n"+
//"                       -h0=[average production height in km\n"+
//"                                   or in g/cm^2 if negative]\n");
//                return;
//                }
//            else if(args[n].startsWith("-m=")){
//                try{
//                    m=(int)Double.parseDouble(args[n].substring(3));
//                }catch(Exception error){
//                    m=0;
//                }
//            }
//            else if(args[n].startsWith("-z0=")){
//                try{
//                    z0=Double.parseDouble(args[n].substring(4));
//                }catch(Exception error){
//                    z0=0;
//                }
//            }
//            else if(args[n].startsWith("-h0=")){
//                try{
//                    h0=Double.parseDouble(args[n].substring(4));
//                }catch(Exception error){
//                    h0=15.5;
//                }
//            }
//        }

//        if(m<0 || m>3) m=0;
//        Output.err.println("Choosing m="+m+" z0="+Output.f(z0)+" h0="+Output.f(h0));
//        CosCorr C = new CosCorr(m, h0, z0);

//        Output.err.println("X0="+Output.f(C.X0)+" h0="+Output.f(C.A.h(C.X0)));
//        Output.err.println("cos, production height, cos, total X, track length");

//        double cos0=C.A.dXdh(C.geth(1, C.X0));
//        for(int i=0; i<=100; i++){
//            x=i/100.;
//            aux=C.geth(x, C.X0);
//            Output.out.println(Output.f(x)+" "+Output.f(aux)+" "+Output.f(C.A.dXdh(aux)/cos0)
//                               +" "+Output.f(C.sum)+" "+Output.f(C.getx(x, aux)));
//        }
//    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates altitude in [km] corresponding to depth X0 in [g/cm^2] as a function of x=cos(zenith angle).
     */

    double CosCorr::geth(double x, double X0){
        double aux;
        st=sqrt(1-x*x);
        sum=I->integrateWithSubstitution(1, -1, this, -2, -X0);
        if(X0<=sum){
            aux=-I->getUpperLimit();
            sum+=I->integrateOpened(-1, 0, this);
        }
        else{
            sum+=I->integrateOpened(-1, 0, this, -X0);
            aux=-I->getUpperLimit();
        }
        return aux;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates atmospheric depth in [g/cm^2] as a function of x=cos(zenith angle).
     */

    double CosCorr::getX(double x){
        st=sqrt(1-x*x);
        return I->integrateOpened(-1, 0, this)+I->integrateWithSubstitution(1, -1, this, -2);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates atmospheric depth in [g/cm^2] above altitude h in [km] as a function of x=cos(zenith angle).
     */

    double CosCorr::getX(double x, double h){
        st=sqrt(1-x*x);
        if(h<1) return I->integrateOpened(-1, -h, this)+I->integrateWithSubstitution(1, -1, this, -2);
        else  return I->integrateWithSubstitution(1, -h, this, -2);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates atmospheric density in [g/cm^2/km] above altitude h in [km] as a function of x=cos(zenith angle).
     */

    double CosCorr::dXdh(double x, double h){
        st=sqrt(1-x*x);
        return function(-h);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates path length in [km] to altitude h in [km] as a function of x=cos(zenith angle).
     */

    double CosCorr::getx(double x, double h){
        double aux=R*x;
        return sqrt(aux*aux+2*R*h+h*h)-aux;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates derivative of mass overburden wrt path distance in [km] as a function of x=cos(zenith angle).
     */

    double CosCorr::function(double x){
        double xs=-x;
        double aux=st/(1+xs/R);
        return A->dXdh(xs)/sqrt(1-aux*aux);
    }


