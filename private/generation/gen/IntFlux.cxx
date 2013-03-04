#include "generation/gen/IntFlux.h"
#include <cmath>
#include <algorithm>
#include <sstream>

using namespace std;

int IntFlux::sM;
double IntFlux::Efr;
double IntFlux::mF;
double IntFlux::gD;
double IntFlux::gEcut;
EarthModel* IntFlux::eM;
NeutrinoTot* IntFlux::nT;
std::vector<Interpolate*> IntFlux::Je;
bool IntFlux::je;
/**
 * Class defines Gaisser-like muon and muon- and electron-neutrino energy spectra with cos*, muon energy loss, and decay corrections.
 */



    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize class with particle type (all mu/mu-/mu+/all nu_mu/nu_mu/~nu_mu/all nu_e/nu_e/~nu_e/),
     * model and ground elevation z0 in [km]. flag switches between two cos* calculation algorithms.
     * h0 is the average production height in km or in g/cm^2 if negative.
     */

    void IntFlux::init(int type, int model, int flag, double h0, double z0){

        z0=0;
        X0=1.148;
        R=EarthModel::R0;
        P=0;
        M=0;
        S=0;
        dG=0;
        xx=0, xs=0, xo=0, xl=0, xn=0;

        double temp1[][5]={
            /* Standard US Atmosphere average parameters */
            { 0.00851164, 0.124534, 0.059761, 2.32876, 19.279 },
            { 0.102573, -0.068287, 0.958633, 0.0407253, 0.817285 },
            { -0.017326, 0.114236, 1.15043, 0.0200854, 1.16714 },
            { 1.3144, 50.2813, 1.33545, 0.252313, 41.0344 }
        };
        double temp2[][5]={
            {0.701, 2.715, 1, 1, 1},
            {0.340, 2.720, 1, 1, 1},
            {0.367, 2.712, 1, 1, 1},
            {0.646, 2.684, 1, 1, 1},
            {0.352, 2.680, 1, 1, 1},
            {0.310, 2.696, 1, 1, 1},
            {0.828, 2.710, 1, 1, 1},
            {0.465, 2.719, 1, 1, 1},
            {0.472, 2.741, 1, 1, 1}
            // { 0.701459, 2.71461, 1, 1, 1 }; all corrections
            // { 0.590968, 2.69506, 1, 0, 0 }; just cos* correction
            // { 0.594662, 2.69671, 0, 0, 0 }; original formula
        };

        for(int i =0;i<4;i++){
            for(int j =0;j<5;j++){
                xpar[i][j]=temp1[i][j];
            }
        }
        for(int i =0;i<9;i++){
            for(int j =0;j<5;j++){
                fpar[i][j]=temp2[i][j];
            }
        }


        D=0;
        Rc=0;

        sM=0;
        Efr=10;
        mF=0;
        gD=0.2;
        gEcut=0;
        jt=false;
        je=false;

        I = new Integral(IROMB, IMAXS, IPREC);
        P=type-1;
        M=model;
        S=flag;
        stringstream ss;

        if(M>=0){
            eC = new ExpCorr(model, h0, z0);
            X0=eC->X0;
            ss<<model<<"_"<<Output::f(h0)<<"_"+Output::f(z0);
            name=ss.str();
        }
        R+=z0;

    }

    IntFlux::IntFlux(int type, int model, int flag, double h0, double z0){

        init(type,model,flag,h0,z0);

    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize the class for spectrum index correction dG and prompt to pion muon component ratio Rc.
     */

    IntFlux::IntFlux(int type, int model, int flag, double h0, double z0, double dG, double Rc, double D){

        init(type,model,flag,h0,z0);
        this->dG=dG;
        this->Rc=Rc;
        this->D=D;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Parametrize h(x) c(x) o(x) l(x).
     */

    void IntFlux::interpolate(string name){
        if(M<0) return;
        int g=5;
        bool flag;
        stringstream ss;
        ss<<".gen_hcolxc_"<<this->name;

        if(Output::raw) ss<<"_raw"; else ss<<"_ascii";
        ss<<".data";
        name=ss.str();
        do {
            if(Output::texi) return;
            flag=false;
            try{
                Output::open(name);

                je=false;
                cerr<<"Parameterizing hcolxC ... "<<endl;
                Je.resize(4);
                for(int i=0; i<4; i++){
                    hcol=i;
                    Je.at(i) = new Interpolate(NUM1, 0, 1, this, g, false, false, false, g, false, false, false);
                }
                je=true;
                cerr<<"done"<<endl;
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
     * Calculates 3 example distrbutions.
     */

//   static void IntFlux::main(String[] args){
//        int w=1, p=0, f=0, m=-1;
//        double h0=-114.8, z0=0;

//        for(int n=0; n<args.length; n++){
//            if(args[n].equals("-help") || args[n].equals("-h") || args[n].equals("--help")){
//                Output.out.println("\n"+
//"This program calculates Energy Integral of atm. lepton fluxes\n"+
//"                       -w=[1-3] which program to run\n"+
//"                       -t=[1-9] particle mu-+/num-+/nue-+\n"+
//"                       -f=[0-1] method of cos* correction\n"+
//"                       -m=[-1/0-3] choose atmosphere model\n"+
//"                       -z0=[ground elevation in km]\n"+
//"                       -h0=[average production height in km\n"+
//"                                   or in g/cm^2 if negative]\n");
//                return;
//                }
//            else if(args[n].startsWith("-t=")){
//                try{
//                    p=(int)Double.parseDouble(args[n].substring(3));
//                }catch(Exception error){
//                    p=1;
//                }
//            }
//            else if(args[n].startsWith("-w=")){
//                try{
//                    w=(int)Double.parseDouble(args[n].substring(3));
//                }catch(Exception error){
//                    w=1;
//                }
//            }
//            else if(args[n].startsWith("-f=")){
//                try{
//                    f=(int)Double.parseDouble(args[n].substring(3));
//                }catch(Exception error){
//                    f=0;
//                }
//            }
//            else if(args[n].startsWith("-m=")){
//                try{
//                    m=(int)Double.parseDouble(args[n].substring(3));
//                }catch(Exception error){
//                    m=-1;
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
//                    h0=-114.8;
//                }
//            }
//        }

//        if(h0==0) m=-1;
//        if(w<1 || w>3) w=1;
//        if(p<1 || p>9) p=1;
//        if(f<0 || f>1) f=0;
//        if(m<-1 || m>3) m=-1;
//        Output.err.println("Choosing w="+w+" t="+p+" f="+f+" m="+m+" z0="+Output.f(z0)+" h0="+Output.f(h0));


//        double x;
//        IntFlux F = new IntFlux(p, m, f, h0, z0);


//        switch(w){
//        case 1:
//            for(int i=0; i<=1000; i++){
//                x=i/1000.;
//                Output.out.println(Output.f(x)+" "+Output.f(F.getIntFlux(600, x)));
//            }
//            break;
//        case 2:
//            x=1;
//            Output.err.println("At cos="+Output.f(x)+", total flux is "+Output.f(F.getIntFlux(600, x))
//                               +", mass overburden is "+Output.f((F.o(x)-F.X0)));
//            for(int i=0; i<=1000000; i++){
//                Output.out.println(Output.f(F.getE(600, x, Math.random())));
//            }
//            break;
//        case 3:
//            for(double e=5;e<5.e5;e*=1.01) Output.out.println(Output.f(e)+" "+Output.f(F.getfl(e, 1)));
//        default:
//        }
//    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates integral flux above energy E in [GeV] at x=cos(zenith angle).
     */

    double IntFlux::getIntFlux(double E, double x){
        return getE(E, x, -1);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates particle energy above threshold E0 in [GeV] at x=cos(zenith angle) as a function of a random number rnd.
     */

    double IntFlux::getE(double E0, double x, double rnd){
        if(rnd<0){
            if(jt) return J->interpolate(x, E0);
            flSet(x);
            return I->integrateWithSubstitution(1, -E0, this, -2);
        }
        else{
            if(jt) return J->findLimit(x, rnd*J->interpolate(x, E0));
            flSet(x);
            I->integrateWithSubstitution(1, -E0, this, -2, rnd);
            return -I->getUpperLimit();
        }
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates differential flux for energy E in [GeV] at x=cos(zenith angle).
     */

    double IntFlux::getfl(double E, double x){
        flSet(x);
        return function(-E);
    }

    //----------------------------------------------------------------------------------------------------//

    void IntFlux::flSet(double x){
        if(mF>0 && P>=3){
            eM->sett(sqrt(1-x*x), 0, -x);
            xn=eM->X(eM->ti, eM->tf)*mF;
        }
        double s=(1-D/R)*sqrt(1-x*x);
        x=sqrt(1-s*s);

        xx=x;
        switch(S){
        case 1:
            xs=sqrt(1-x*x)/(1+h(x)/R);
            xs=sqrt(1-xs*xs);
            break;
        case 0:
        default:
            xs=c(x);
        }
        if(P<3){
            xo=fpar[P][3]*o(x)-X0; if(xo<0) xo=0;
            xl=fpar[P][4]*l(x)*1.e2*MMU/(LMU*SPEED);
        }
    }

    //----------------------------------------------------------------------------------------------------//

    double IntFlux::h(double x){
        if(M>=0) return je?max(Je.at(0)->interpolate(x), 0.):eC->getExp(1, x);
        double s=sqrt(1-x*x);
        double y=s/(1+xpar[0][0]);
        return xpar[0][4]*(1+xpar[0][2]*pow(s, xpar[0][3]))/
            pow(1-(y*y), xpar[0][1]);
    }

    //----------------------------------------------------------------------------------------------------//

    double IntFlux::c(double x){
        if(M>=0) return je?max(Je.at(1)->interpolate(x), 0.):1/eC->getExp(2, x);
        return sqrt((x*x+fpar[P][2]*
                          (xpar[1][0]*xpar[1][0]+
                           xpar[1][1]*pow(x, xpar[1][2])+
                           xpar[1][3]*pow(x, xpar[1][4])))/
                         (1+fpar[P][2]*(xpar[1][0]*xpar[1][0]+
                                        xpar[1][1]+xpar[1][3])));
    }

    //----------------------------------------------------------------------------------------------------//

    double IntFlux::o(double x){
        if(M>=0) return je?max(Je.at(2)->interpolate(x), 0.):eC->C->getX(x)/1.e2;
        return 1/(xpar[2][0]+xpar[2][1]*pow(x, xpar[2][2])+
                  xpar[2][3]*pow(1-x*x, xpar[2][4]));
    }

    //----------------------------------------------------------------------------------------------------//

    double IntFlux::l(double x){
        if(M>=0) return je?max(Je.at(3)->interpolate(x), 0.):eC->getExp(3, x);
        return 1e3/(xpar[3][0]+xpar[3][1]*pow(x, xpar[3][2])+
                                  xpar[3][3]*pow(1-x*x, xpar[3][4]));
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Energy distribution as a function of energy x in [GeV]; zenith angle is set in a private variable - interface to Integral.
     */

    double IntFlux::function(double x){
        double ei, ef, xn, sm;
        ef=-x;

        if(mF>0 && P>=3){
            bool nu=(P==4 || P==7);
            double tI=(P==8)?nT->dSdy(ef):0;
            tI+=nT->dSdy(ef, true, nu);
            tI+=nT->dSdy(ef, false, nu);
            xn=this->xn*tI;
            if(xn>1) xn=1;
        }
        else xn=1;

        if(sM>0){
            const double Nmax=2465;
            double N;
            switch(sM){
            case 1: N=2445; break;
            case 2: N=2300; break;
            case 3: N=2115; break;
            default: N=min((double)sM, Nmax-1);
            }
            sm=exp(-(1.15+14.9*pow(1-N/Nmax, 1.12))/(0.97+ef*Efr));
        }
        else sm=1;

        if(gEcut>0) sm/=1+exp(-log(ef/(gEcut*Efr))/(gD*LOG10));

        switch((int)floor(P/3)){
        case 1:
            {
                const double E1=121;
                const double E2=897;
                const double pK=0.213;
                const double A=2.85e-2;
                ei=ef;
                return sm*xn*A*fpar[P][0]*pow(ei,-(fpar[P][1]+dG))*(1/(1+6*ei*xs/E1)+pK/(1+1.44*ei*xs/E2));
            }
        case 2:
            {
                const double E1=121;
                const double E2=897;
                const double E3=194;
                const double A=2.4e-3;

                double xi, a, b;
                ei=ef;

                if(xx<0.3){
                    a=0.11-2.4*xx;
                    b=-0.22+0.69*xx;
                }
                else{
                    a=-0.46-0.54*xx;
                    b=-0.01+0.01*xx;
                }
                xi=a+b*log(ei)/LOG10;
                return sm*xn*A*fpar[P][0]*pow(ei,-(fpar[P][1]+dG))*(0.05/(1+1.5*ei*xs/E2)+0.185/(1+1.5*ei*xs/E3)+
                                                                   11.4*pow(ei, xi)/(1+1.21*ei*xs/E1));
            }
        case 0:
        default:
            {
                const double E1=115;
                const double E2=850;
                const double pK=0.054;
                const double a=0.260;
                const double b=0.360e-3;
                const double as=-0.522;
                const double bs=8.893e-3;
                const double A=0.14;

                double zfi, zff, p, aux, aux1, aux2, auf1, auf2, f1, f2, f3, f4, de;
                ei=((a+b*ef)*exp(b*xo)-a)/b;

                aux1=(bs*a-as*b);
                aux1*=2*aux1;
                aux2=2*b*b*b;
                zfi=(bs*ei*b*(-2*bs*a+4*as*b+bs*ei*b)+aux1*log(a+ei*b))/aux2;
                zff=(bs*ef*b*(-2*bs*a+4*as*b+bs*ef*b)+aux1*log(a+ef*b))/aux2;

                aux=1.1*xs;
                auf1=aux/E1;
                auf2=aux/E2;
                aux1=1/(1+auf1*ei);
                aux2=1/(1+auf2*ei);
                aux=aux1+pK*aux2;

                f1=pow(ei, -fpar[P][1]-1)*aux;
                f2=pow(ei, -fpar[P][1]-1)*(-(fpar[P][1]+1)*aux/ei-auf1*aux1*aux1-pK*auf2*aux2*aux2);
                f3=(-fpar[P][1]-1)*f2/ei+pow(ei, -fpar[P][1]-1)
                    *((fpar[P][1]+1)*aux/(ei*ei)
                      +((fpar[P][1]+1)/ei)*(auf1*aux1*aux1+pK*auf2*aux2*aux2)
                      +2*auf1*auf1*aux1*aux1*aux1+2*pK*auf2*auf2*aux2*aux2*aux2);
                aux=f2/f1;
                f4=f3/f1-aux*aux;

                aux1=as+bs*ef;
                aux2=as+bs*ei;
                de=exp(b*xo)-(f2/(2*f1))*(aux1*aux1-aux2*aux2)/(a+b*ef)-(f4/2)*(zff-zfi);

                ei=ei-(f2/(2*f1))*(zff-zfi);

                p=exp(-xl/ei);
                return sm*A*fpar[P][0]*p*de*pow(ei,-(fpar[P][1]+dG))*(Rc<0?-Rc:(1/(1+1.1*ei*xs/E1)+pK/(1+1.1*ei*xs/E2)+Rc));
            }
        }
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Parametrizes the integral of this class.
     */

    void IntFlux::interpolate(double E0){
        int g=5;
        double e_hi=BIGENERGY*1.e-3;
        double e_low=E0;
        jt=false;
        J = new Interpolate(NUM1, 0, 1, NUM1, e_low, e_hi, this, g, false, false, false, g, false, false, true, g, false, false, true);
        jt=true;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * 2d parametrization - interface to Interpolate
     */

    double IntFlux::functionInt(double x, double E){
        flSet(x);
        return I->integrateWithSubstitution(1, -E, this, -2);
    }

    //----------------------------------------------------------------------------------------------------//



    /**
     * 1d parametrization - interface to Interpolate
     */

    double IntFlux::functionInt(double x){
        switch(hcol){
        case 0: return h(x);
        case 1: return c(x);
        case 2: return o(x);
        case 3: return l(x);
        default: return 0;
        }
    }



