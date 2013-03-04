#include "generation/gen/ExpCorr.h"
#include <cmath>

/**
 * This class calculates atmosphere-dependent quantities for a given average first interaction depth.
 */


    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize class with atmospheric model and ground elevation z0 in [km].
     * h0 is the average production height in km or in g/cm^2 if negative.
     */

    ExpCorr::ExpCorr(int model, double h0, double z0){

        func=0;
        xs=1;
        Norm=1;
        h0=19.279;
        X0=114.8;

        C = new CosCorr(model, h0, z0);
        I = new Integral(IROMB, IMAXS, IPREC2);
        if(h0<0){
            this->X0=-h0;
            Norm=getExp(0, 1);
            this->h0=getExp(1, 1);
        }
        else{
            this->h0=h0;
            this->X0=C->A->X(h0-z0);
            Norm=getExp(0, 1);
        }
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Prints zenith angle profile of production height, cos*, total X, track length.
     */

//    static void ExpCorr::main(String[] args){
//        int w=0, m=1;
//        double h0=15.5, z0=0;
//        ExpCorr E;

//        for(int n=0; n<args.length; n++){
//            if(args[n].equals("-help") || args[n].equals("-h") || args[n].equals("--help")){
//                Output.out.println("\n"+
//"This class calculates mean production height, cos(th*), etc.\n"+
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
//                    m=1;
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

//        switch(w){
//        case 1: /* Standard US Atmosphere elevation 2.834 m*/
//            E = new ExpCorr(1, 15.5, 2.834);
//            break;
//        case 2: /* South Pole Atmosphere Oct 01 */
//            E = new ExpCorr(2, -114.8, 2.834);
//            break;
//        case 3: /* South Pole Atmosphere Jul 01 */
//            E = new ExpCorr(3, -114.8, 2.834);
//            break;
//        case 0: /* Standard US Atmosphere parameters: 1, 15.5, 0 */
//        default:
//            E = new ExpCorr(m, h0, z0);
//        }
//        Output.err.println("X0="+Output.f(E.X0)+" h0="+Output.f(E.getExp(1, 1)));
//        Output.err.println("cos, production height, cos, total X, track length");

//        for(int i=0; i<=100; i++){
//            double x=i/100.;
//            Output.out.println(Output.f(x)+" "+Output.f(E.getExp(1, x))
//                               +" "+Output.f(1/E.getExp(2, x))+" "+
//                               Output.f(E.C.getX(x))+" "+Output.f(E.getExp(3, x)));
//        }
//    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Common function for calculation of atmospheric averages:
     * <ol>
     * <li> func=1: muon production height
     * <li> func=2: 1/cos~1/rho
     * <li> func=3: muon track length
     * <li> all others: normalization.
     * </ol>
     */

    double ExpCorr::getExp(int func, double x){
        this->func=func;
        xs=x;
        double aux=0;
        aux=I->integrateOpened(0, X0, this)+I->integrateWithSubstitution(X0, -1, this, 2);
        return aux/Norm;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Common function for calculation of atmospheric averages - interface to Integral.
     */

    double ExpCorr::function(double X){
        double aux=1;
        double rho=exp(-X/X0);
        switch(func){
        case 0: aux=1; break;                          // normalization
        case 1: aux=C->geth(xs, X); break;              // muon production height
        case 2: aux=C->A->dXdh(C->geth(1, X))/C->A->dXdh(C->geth(xs, X));
            break;                                     // 1/cos~1/rho
        case 3: aux=C->getx(xs, C->geth(xs, X)); break;  // muon track length
        default:
            break;
        }
        return rho*aux;
    }


