#include "generation/gen/AtmFlux.h"
#include <cmath>
#include <algorithm>
#include "PROPOSAL/methods.h"
#include <sstream>
#include "PROPOSAL/Decay.h"


using namespace std;

/**
 * This is the main lepton generator class.
 */

AtmFlux::AtmFlux(){

    R=400;
    tF=0;
    nF=1, mF=0, xn=1;
    tT=1;
    t0=0;
    evt=0;
    OSC=true;
    GEN=false;
    PROP=false;

    R0=eM->R0*1.e3;

    lept=0;
    f2k=true;
    tdir="";
    dtct="";
    string lala[]={"mu-", "mu+", "nu_mu", "~nu_mu", "nu_e", "~nu_e", "nu_tau", "~nu_tau", "e-", "e+", "tau-", "tau+", "hadr"};
    for (int i=0; i<13; ++i){
        name[i]=lala[i];
    }


    D=1730;
    jt=false;

}

    //----------------------------------------------------------------------------------------------------//

    /**
     * Class initializer, command-line option parser. Call with "-help" to list all options.
     */
long AtmFlux::initAtmFlux(deque<std::string> args, bool f2k){
        this->f2k=f2k;
        string param="AtmFlux";
        int i;
        long num=10;

        double Emu=300;
        double Enm=100;
        double Ene=100;
        double Ent=100;
        double Amu=1, Anm=1, Ane=1;
        double Gmu=0, Gnm=0, Gne=0;
        double ant=0, anm=0, ane=0;
        double gnt=2, gnm=2, gne=2;
        double prompt=0;

        int w=1, f=0, m=-1;
        double h0=-114.8, z0=2.834, b0=0.024;
        int sM=0;
        double Efr=10;
        double gE=0, gD=0.2;

        long seed=0;
        bool SEED=false;
        stringstream paramss;
        stringstream pbadss;
        stringstream enamess;
        int bnum=0;
        string pbad="";
        bool pflag;

        for(i=0; i<args.size(); i++){
            pflag=true;
            if(args.at(i).compare("-help")==0 || args.at(i).compare("-h")==0 || args.at(i).compare("--help")==0){
                cerr<<"\n"<<
"This program generates atmospheric lepton fluxes at 0.1-4000 TeV\n"<<
"                       -Emu=[100-1.e5 GeV] muon low energy threshold\n"<<
"                       -Enm=[GeV] muon neutrino low energy threshold\n"<<
"                       -Ene=[GeV] electron neutrino energy threshold\n"<<
"                       -Ent=[GeV] tau neutrino  low energy threshold\n"<<
"                       -A[mu/nm/ne]=[normalization  correction]\n"<<
"                       -G[mu/nm/ne]=[spectral index correction]\n"<<
"                       -a[nm/ne/nt]=[    for the ad-hoc component    ]\n"<<
"                       -g[nm/ne/nt]=[1e-7/(GeV sr s cm^2) a(E/GeV)^-g]\n"<<
"                       -prompt=[R_c]  ratio of prompt to pion muons\n"<<
"                       -R=[radius] of the detector in meters\n"<<
"                       -D=[depth]  of the detector in meters\n"<<
"                       -N=[number] of events to generate\n"<<
"                       -f=[0-1]  method of cos* correction\n"<<
"                       -m=[-1/0-3] choose atmosphere model\n"<<
"                       -z0=[ground  elevation in km]\n"<<
"                       -b0=[bedrock elevation in km]\n"<<
"                       -h0=[average production height in km\n"<<
"                                   or in g/cm^2 if negative]\n"<<
"                       -nF=[neutrino cross sections multiplic. factor]\n"<<
"                       -mF=[max. mass overburden/neutrino int. length]\n"<<
"                       -sM=[0-4]  Solar modulation:  none/min/mid/max\n"<<
"                       -Efr=[cr primary to neutrino energy fraction]\n"<<
"                       -gE=[GeV]  artificial geomagnetic cutoff value\n"<<
"                       -gD=[decades]  span of the response function\n"<<
"                       -noosc   disable mu < - > tau oscillations\n"<<
"                       -nlow=[GeV]   lowest neutrino energy\n"<<
"                       -gen   use as generator  only\n"<<
"                       -prop  use as propagator only\n"<<
"                       -seed=[integer] random number generator seed\n"<<
"mmc help lines follow:"<<endl;

                //???????????????????????????????????????????????????????????????????????????????????????????
                Amanda *a =new Amanda();
                a->setp(args);
                return -1;
            }
            else if(args.at(i).compare("-gen")==0){
                GEN=true;
            }
            else if(args.at(i).compare("-prop")==0){
                PROP=true;
            }
            else if(args.at(i).compare("-noosc")==0){
                OSC=false;
            }
            else if(StartsWith(args.at(i),"-Emu=")){
                    Emu=atof(args.at(i).substr(5).c_str());

                   // Emu=600;

            }
            else if(StartsWith(args.at(i),"-Enm=")){

                     Enm=atof(args.at(i).substr(5).c_str());
                    //Enm=600;

            }
            else if(StartsWith(args.at(i),"-Ene=")){

                     Ene=atof(args.at(i).substr(5).c_str());
                   // Ene=600;

            }
            else if(StartsWith(args.at(i),"-Ent=")){

                     Ent=atof(args.at(i).substr(5).c_str());
                    //Ent=600;

            }
            else if(StartsWith(args.at(i),"-Amu=")){

                     Amu=atof(args.at(i).substr(5).c_str());
                     //Amu=1;

            }
            else if(StartsWith(args.at(i),"-Anm=")){
                    Anm=atof(args.at(i).substr(5).c_str());
                    //Anm=1;

            }
            else if(StartsWith(args.at(i),"-Ane=")){
                    Ane=atof(args.at(i).substr(5).c_str());
                    //Ane=1;

            }
            else if(StartsWith(args.at(i),"-Gmu=")){
                    Gmu=atof(args.at(i).substr(5).c_str());
                    //Gmu=0;

            }
            else if(StartsWith(args.at(i),"-Gnm=")){
                    Gnm=atof(args.at(i).substr(5).c_str());
                    //Gnm=0;

            }
            else if(StartsWith(args.at(i),"-Gne=")){
                    Gne=atof(args.at(i).substr(5).c_str());
                    //Gne=0;

            }
            else if(StartsWith(args.at(i),"-ant=")){
                    ant=atof(args.at(i).substr(5).c_str());
                    //ant=0;

            }
            else if(StartsWith(args.at(i),"-anm=")){
                    anm=atof(args.at(i).substr(5).c_str());
                    //anm=0;

            }
            else if(StartsWith(args.at(i),"-ane=")){
                    ane=atof(args.at(i).substr(5).c_str());
                    //ane=0;

            }
            else if(StartsWith(args.at(i),"-gnt=")){
                    gnt=atof(args.at(i).substr(5).c_str());
                    //gnt=2;

            }
            else if(StartsWith(args.at(i),"-gnm=")){
                gnm=atof(args.at(i).substr(5).c_str());
                    //gnm=2;

            }
            else if(StartsWith(args.at(i),"-gne=")){
                    gne=atof(args.at(i).substr(5).c_str());
                    //gne=2;

            }
            else if(StartsWith(args.at(i),"-prompt=")){
                    prompt=atof(args.at(i).substr(8).c_str());
                    //prompt=0;

            }
            else if(StartsWith(args.at(i),"-R=")){
                    R=atof(args.at(i).substr(3).c_str());
                    //R=10;

            }
            else if(StartsWith(args.at(i),"-D=")){
                    D=atof(args.at(i).substr(3).c_str());
                    //D=0;

            }
            else if(StartsWith(args.at(i),"-N=")){
                    num=(long)atof(args.at(i).substr(3).c_str());
                    //num=1;

            }
            else if(StartsWith(args.at(i),"-f=")){
                    f=(int)atof(args.at(i).substr(3).c_str());
                    //f=0;

            }
            else if(StartsWith(args.at(i),"-m=")){
                    m=(int)atof(args.at(i).substr(3).c_str());
                    //m=-1;

            }
            else if(StartsWith(args.at(i),"-z0=")){
                    z0=atof(args.at(i).substr(4).c_str());
                    //z0=0;

            }
            else if(StartsWith(args.at(i),"-b0=")){
                    b0=atof(args.at(i).substr(4).c_str());
                    //b0=0;

            }
            else if(StartsWith(args.at(i),"-h0=")){
                    h0=atof(args.at(i).substr(4).c_str());
                    //h0=-114.8;

            }
            else if(StartsWith(args.at(i),"-nF=")){
                    nF=atof(args.at(i).substr(4).c_str());
                    //nF=1;

            }
            else if(StartsWith(args.at(i),"-mF=")){
                    mF=atof(args.at(i).substr(4).c_str());
                    //mF=0;

            }
            else if(StartsWith(args.at(i),"-sM=")){
                    sM=(int)atof(args.at(i).substr(4).c_str());
                    //sM=0;

            }
            else if(StartsWith(args.at(i),"-Efr=")){
                    Efr=atof(args.at(i).substr(5).c_str());
                    //Efr=10;

            }
            else if(StartsWith(args.at(i),"-gE=")){
                    gE=atof(args.at(i).substr(4).c_str());
                    //gE=0;

            }
            else if(StartsWith(args.at(i),"-gD=")){
                    gD=atof(args.at(i).substr(4).c_str());
                    //gD=0.2;

            }
            else if(StartsWith(args.at(i),"-nlow=")){
                    nlow_=atof(args.at(i).substr(6).c_str())*1.e3;
                    //nlow=Me;

            }
            else if(args.at(i).compare("-tau")==0){
                    args.at(i)="-ignore-tau";
            }
            else if(StartsWith(args.at(i),"-seed=")){
                    seed=(long)atof(args.at(i).substr(6).c_str());
                    seed=0;
                    SEED=true;
            }
            else if(StartsWith(args.at(i),"-tdir=")){
                    tdir=args.at(i).substr(6)+"/";
            }
            else{
                bnum++;
                pbadss<<(bnum>1?",":"")<<" \""<<args.at(i)<<"\"";
                pbad+=pbadss.str();
                pbadss.clear();
                pbadss.str("");
                pflag=false;
            }
            if(pflag){
                paramss<<" "<<args.at(i);

                param+=paramss.str();
                paramss.clear();
                paramss.str("");
            }
        }

        cerr<<Output::version<<endl;
        cerr<<"Running \""<<param<<"\""<<endl;
        if(bnum>0){
            paramss<<" (not used:"<<pbad<<")";
            param+=paramss.str();
            paramss.clear();
            paramss.str("");

            bnum==1?pbadss<<" is":pbadss<<"s"<<" are";
            pbad+=pbadss.str();
            pbadss.clear();
            pbadss.str("");
            cerr<<"Warning: Parameter"<<pbad<<" not recognized"<<endl;
        }

        if(Emu<MMU*1.e-3) Emu=MMU*1.e-3;
        if(Enm<ME*1.e-3) Enm=ME*1.e-3;
        if(Ene<ME*1.e-3) Ene=ME*1.e-3;
        if(Ent<ME*1.e-3) Ent=ME*1.e-3;

        if(Amu<0) Amu=0;
        if(Anm<0) Anm=0;
        if(Ane<0) Ane=0;
        if(ant<0) ant=0;
        if(anm<0) anm=0;
        if(ane<0) ane=0;

        if(Gmu<-1) Gmu=0;
        if(Gnm<-1) Gnm=0;
        if(Gne<-1) Gne=0;
        if(gnt<=1) gnt=2;
        if(gnm<=1) gnm=2;
        if(gne<=1) gne=2;

        if(Amu==0 && Anm==0 && Ane==0 && ant==0 && anm==0 && ane==0) Amu=1;
        if(GEN) PROP=false;

        if(h0==0) m=-1;
        if(f<0 || f>1) f=0;
        if(m<-1 || m>3) m=-1;
        cerr<<"Choosing f="<<f<<" m="<<m<<" z0="<<Output::f(z0)<<" km  b0="<<
                           Output::f(b0)<<" km  h0="<<Output::f(h0)<<" "<<(h0>0?"km":"g/cm^2")<<endl;
        cerr<<"Detector with radius R="<<Output::f(R)<<" m is at depth="<<Output::f(D)<<" m"<<endl;
        cerr<<"Energy thresholds: Emu="<<Output::f(Emu)<<" GeV, Enm="<<Output::f(Enm)<<" GeV, Ene="<<Output::f(Ene)<<" GeV, Ent="<<Output::f(Ent)<<" GeV"<<endl;

        if(nF<=0) nF=1;
        if(nF>1) cerr<<"Use nF>1 to speed up the calculation for energy thresholds << 10 TeV"<<endl;
        else if(nF<1) cerr<<"Setting nF value to "<<Output::f(nF)<<" < 1. Your time would probably be much better spent elsewhere"<<endl;

        if(mF<1) mF=0;
        if(GEN || PROP) mF=0;
        if(mF>0) cerr<<"Using re-weighting technique. This is an approximation! Only use mF>>1, e.g. mF=10"<<endl;

        if(sM<0) sM=0;
        if(PROP) sM=0;
        if(Efr<1) Efr=10;
        if(sM>0) cerr<<"Solar modulation enabled with parameters: sM="<<sM<<", Efr="<<Output::f(Efr)<<" GeV"<<endl;
        if(gE<0) gE=0;
        if(gD<=0) gD=0.2;
        if(gE>0) cerr<<"Geomagnetic cutoff enabled with parameters: gE="<<Output::f(gE)<<
                                    " GeV, gD="<<Output::f(gD)<<(sM>0?"":", Efr=")<<Output::f(Efr)+" GeV"<<endl;

        Decay::flag=true;
        Output::RecDec=true;
        if(!GEN){
            mP = new Amanda();
            tP = new Amanda();
            mP->dtct=dtct;
            tP->dtct=dtct;
            eM = new EarthModel(m<0?0:m, z0, b0, D*1.e-3);
            nT = new NeutrinoTot();
        }

        //rand = new Random();
        R0+=z0*1.e3;

        if(f2k) cerr<<"V 2000.1.2"<<endl;

        string ename="";

        if(!GEN){
            mP->setup(args);
            if(f2k) mP->history();
            deque<string> argt;
            argt.resize(args.size()+1);
            for(i=0; i<args.size(); i++) argt.at(i)=args.at(i);
            argt.at(args.size()-1)="-tau";
            tP->setup(argt);
            if(f2k) tP->history();
            switch(mP->gdet){
            case 1:
                length=mP->height;
                radius=mP->length/2;
                break;
            case 2:
                length=mP->radius*2;
                radius=mP->radius;
                break;
            case 0:
            default:
                length=mP->length;
                radius=mP->radius;
                break;
            }

            eM->num=mP->medianum;
            eM->mRa = new double[eM->num];
            eM->mRo = new double[eM->num];
            eM->mRa[0]=R0*1.e-3;
            eM->mRo[0]=mP->srho[0];
            enamess<<"_"<<eM->mRo[0];
            ename=enamess.str();
            enamess.clear();
            enamess.str("");
            for(i=1; i<eM->num; i++){
                if(mP->medt[i]==1 && abs((mP->sphr[i]-mP->sphz[i])-(R0-D))<R0*HALF_PRECISION){
                    eM->mRa[i]=(mP->sphz[i]+R0-D)*1.e-3;
                    eM->mRo[i]=mP->srho[i];
                    enamess<<"_"<<Output::f(mP->sphz[i])<<"-"<<Output::f(eM->mRo[i]);
                    ename+=enamess.str();
                    enamess.clear();
                    enamess.str("");
                }
                else eM->mRo[i]=0;
               // cout<<i<<"\t eM->mRo[i]="<<eM->mRo[i]<<endl;
            }
        }

        string name=tdir+dtct;

        if(!PROP){
            I = new Integral(IROMB, IMAXS, IPREC);
            int type[] = {2, 3, 5, 6, 8, 9};
            tf = new double[6];
            af = new double[6];
            sF = new double[6];
            Ecut = new double[8];
            A = new double[6];
            G = new double[6];
            a = new double[6];
            g = new double[6];
            Ecut[0]=Emu; A[0]=Amu; G[0]=Gmu; a[0]=ant*1.e-7; g[0]=gnt;
            Ecut[1]=Emu; A[1]=Amu; G[1]=Gmu; a[1]=ant*1.e-7; g[1]=gnt;
            Ecut[2]=Enm; A[2]=Anm; G[2]=Gnm; a[2]=anm*1.e-7; g[2]=gnm;
            Ecut[3]=Enm; A[3]=Anm; G[3]=Gnm; a[3]=anm*1.e-7; g[3]=gnm;
            Ecut[4]=Ene; A[4]=Ane; G[4]=Gne; a[4]=ane*1.e-7; g[4]=gne;
            Ecut[5]=Ene; A[5]=Ane; G[5]=Gne; a[5]=ane*1.e-7; g[5]=gne;
            Ecut[6]=Ent; Ecut[7]=Ent;
            iF.resize(6);
            for(i=0; i<6; i++) iF.at(i) = new IntFlux(type[i], m, f, h0, z0, G[i], prompt, D*1.e-3);
            iF[0]->interpolate(name);
        }

        CteqPDF::tdir=tdir;
        if(!GEN) nT->interpolate(name);
        stringstream namess;

        namess<<".gen_";
        if(m==-1) namess<<"t"; else namess<<m<<"_"<<f<<"_"<<Output::f(h0);
        if(!PROP && !(Amu==0 && Anm==0 && Ane==0)){
            namess<<"_"<<Output::f(Ecut[0])<<"_"<<Output::f(Ecut[2])<<"_"<<Output::f(Ecut[4]);
            namess<<"_"<<Output::f(G[0])<<"_"<<Output::f(G[2])<<"_"<<Output::f(G[4]);
            if(prompt!=0) namess<<"_"<<Output::f(prompt);
        }
        namess<<"_"<<Output::f(z0)<<"_"<<Output::f(D);
        if(!GEN) namess<<"_"<<Output::f(b0)<<ename;
        else namess<<"_go";
        namess<<"_"<<Output::f(mF);
        if(mF>0){
            if(nlow_!=ME) namess<<"_l"<<Output::f(nlow_);
            if(ebig_!=BIGENERGY) namess<<"_b"<<Output::f(ebig_);
        }
        if(sM>0) namess<<"_s"<<sM<<"-"<<Output::f(Efr);
        if(gE>0) namess<<"_g"<<Output::f(gE)<<"-"<<Output::f(gD);
        if(Output::raw) namess<<"_raw"; else namess<<"_ascii";
        namess<<".data";

        name+=namess.str();
        namess.clear();
        namess.str("");
        bool flag;
        do {
            if(Output::texi) return -1;
            flag=false;
            try{
                Output::open(name);

                if(!GEN) eM->interpolate();
                if(mF>0){
                    IntFlux::mF=mF;
                    IntFlux::eM=eM;
                    IntFlux::nT=nT;
                }
                if(sM>0){
                    IntFlux::sM=sM;
                    IntFlux::Efr=Efr;
                }
                if(gE>0){
                    IntFlux::gEcut=gE;
                    IntFlux::gD=gD;
                }
                if(!PROP && !(Amu==0 && Anm==0 && Ane==0)){
                    cerr<<"Parameterizing intflX ... "<<endl;
                    for(i=0; i<6; i++) iF.at(i)->interpolate(Ecut[i]);
                    cerr<<"done"<<endl;
                    interpolate();
                    cerr<<"Finished parameterizations"<<endl;
                }
                Output::close();
            }catch(int a){
                throw 0;
                flag=true;
                Output::Delete(name);
            }
        } while(flag);

        if(!PROP){
            for(i=0; i<6; i++){
                lept=i; tf[i]=A[i]*getTotFlux(1);
                af[i]=2*PI*a[i]*pow(Ecut[i<2?i+6:i], 1-g[i])/(g[i]-1);
                tF+=tf[i]+af[i]; sF[i]=tF;
            }
            for(i=0; i<6; i++) sF[i]/=tF;
            tT=nF/(tF*2.e4*PI*R*R);
        }

        if(f2k){
            cout<<"HI  "<<Output::version<<"\nHI  "<<param<<endl;
            if(!PROP){
                cout<<"HI  Total Flux is "<<Output::f(tF)<<" cm^-2s^-1 ( "<<Output::f(tf[0]+tf[1])<<" "<<Output::f(tf[2]+tf[3])<<" "<<
                                   Output::f(tf[4]+tf[5])<<" "<<Output::f(af[0]+af[1])<<" "<<Output::f(af[2]+af[3])<<" "<<Output::f(af[4]+af[5])<<" )"<<endl;
                cout<<"HI  corresponds to the energy cutoffs of Emu="<<
                                   Output::f(Ecut[0])<<" Enm="<<Output::f(Ecut[2])<<" Ene="<<Output::f(Ecut[4])<<" Ent="<<Output::f(Ecut[6])<<" GeV"<<endl;
                cout<<"HI  Amu="<<Output::f(A[0])<<" Anm="<<Output::f(A[2])<<" Ane="<<Output::f(A[4])
                                   <<"  Gmu="<<Output::f(G[0])<<" Gnm="<<Output::f(G[2])<<" Gne="<<Output::f(G[4])<<endl;
                cout<<"HI  ant="<<Output::f(a[0])<<" anm="<<Output::f(a[2])<<" ane="<<Output::f(a[4])
                                   <<"  gnt="<<Output::f(g[0])<<" gnm="<<Output::f(g[2])<<" gne="<<Output::f(g[4])<<endl;
                if(prompt>0) cout<<"HI  ratio of prompt to pion muons is "<<Output::f(prompt)<<endl;
                else if(prompt<0) cout<<"HI  only prompt muons will be generated: "<<Output::f(prompt)<<endl;
            }
            cout<<"HI  Detector with radius R="<<Output::f(R)<<" m is at depth="<<Output::f(D)<<" m"<<endl;
            cout<<"HI  Neutrino cross section multiplicative factors are nF="<<Output::f(nF)<<", mF="<<Output::f(mF)<<endl;
            if(sM>0) cout<<"HI  Solar modulation enabled with parameters: sM="<<sM<<", Efr="<<Output::f(Efr)<<" GeV"<<endl;
            if(gE>0) cout<<"HI  Geomagnetic cutoff enabled with parameters: gE="<<Output::f(gE)<<
                                        " GeV, gD="<<Output::f(gD)<<(sM>0?"":", Efr=")<<Output::f(Efr)<<" GeV"<<endl;
            if(!PROP) cout<<"HI  "<<Output::f(num)<<" showers correspond to the lifetime of "<<Output::f(tT*num)+" seconds"<<endl;
        }

        if(SEED){
            if(f2k) cout<<"HI  random number generator seed is set at "<<seed<<endl;
            MathModel::set_seed(seed);
        }

        return num;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * begins the f2k stream.
     */

    void AtmFlux::beginAtmFlux(){
        if(!GEN) if(f2k){
            mP->defUSER();
            tP->defUSER();
            cout<<"USER_DEF almc_e Name Eini Elost Theta"<<endl;

            mP->initBufs();
            tP->initBufs();
        }
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * concludes the f2k stream.
     */

    void AtmFlux::endAtmFlux(){
        if(GEN) cerr<<"generated "<<evt<<" events"<<endl;
        else{
            cerr<<"generated "<<evt<<" events, "<<mP->tracks<<" muon tracks, "<<tP->tracks<<" tau tracks"<<endl;
            cerr<<"average "<<(mP->vertices+tP->vertices)/evt<<" vertices, missed volume: "<<(mP->missed+tP->missed)<<endl;
        }
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * F2k stream lepton generator.
     */

//    static void AtmFlux::main(String[] args){
//        AtmFlux A = new AtmFlux();
//        A.dtct=".atmflux";
//        long num=A.initAtmFlux(args, true);

//        if(A.PROP){
//            A.mmcf2k();
//        }
//        else if(num>=0){
//            A.beginAtmFlux();
//            Output.out.println("TBEGIN ? ? ?");
//            for(long i=0; i<num; i++) A.findNext();
//            cout<<"TEND ? ? ?"<<endl;
//            cout<<"END"<<endl;
//            A.endAtmFlux();
//        }
//    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * I3m initializer.
     */

    void AtmFlux::setup(deque<string> args){
        dtct=".icecube";
        initAtmFlux(args, false);
       // I3hist = new Vector(Output::HISTSIZE);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * I3m initializer.
     */

    void AtmFlux::setup(string args){
        setup(Output::splitString(args," \t"));
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * I3m function - creates particle set of the next event.
     */

    vector<PROPOSALParticle*> AtmFlux::createNext(){

        vector<PROPOSALParticle*> I3p;
        if(f2k){
            I3p.resize(0);
             return I3p;
        }
        I3hist.clear();
        findNext();
        I3p.resize(I3hist.size());
        for(int i=0; i<I3hist.size(); i++) I3p.at(i)=I3hist.at(i);
        return I3p;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * creates the next event.
     */

    void AtmFlux::findNext(){
        evt++;
        if(!GEN) if(f2k){
            mP->iniUSER();
            tP->iniUSER();
        }

        double aux, r, q, fx, fy;
        double ct, st, c1, s1, sa, ca, sp, cp, b, c, tt, tts, ttn;
        double axx, axy, axz, ayx, ayy, ayz, azx, azy, azz, rx, ry, rz;
        double E, x, y, z, t;

        aux=MathModel::RandomDouble();
        lept=0; for(int i=0; i<6; i++) if(aux<sF[i]){ lept=i; break; }

        if(MathModel::RandomDouble()<af[lept]/(af[lept]+tf[lept])){
            ct=MathModel::RandomDouble();
            E=Ecut[lept<2?lept+6:lept]*pow(MathModel::RandomDouble(), 1/(1-g[lept]));
            if(lept<2) lept=lept+6;
        }
        else{
            ct=getTotFlux(1, MathModel::RandomDouble());
            E=iF.at(lept)->getE(Ecut[lept], ct, MathModel::RandomDouble());
        }
        th=acos(ct);
        phi=MathModel::RandomDouble()*2*PI;

        r=R*sqrt(MathModel::RandomDouble());
        q=MathModel::RandomDouble()*2*PI;
        x=r*cos(q);
        y=r*sin(q);

        st=sqrt(1-ct*ct);  // angle in the detector frame -> angle at the surface
        if(st==0){
            s1=0;
            c1=ct;
        }
        else{
            s1=st*(1-D/R0);
            c1=sqrt(1-s1*s1);
        }

        x/=c1;

        fx=x*cos(phi)-y*sin(phi);
        fy=x*sin(phi)+y*cos(phi);

        if(MathModel::RandomDouble()>0.5){ ct=-ct; th=PI-th; }

        sa=st*c1-ct*s1; if(sa>1) sa=1;  // angle difference between the two frames
        ca=ct*c1+st*s1; if(ca>1) ca=1; else if(ca<-1) ca=-1;

        cp=cos(phi);
        sp=sin(phi);

        axx=cp*cp*(ca-1)+1;   // transformation matrix computation
        axy=cp*sp*(ca-1);
        axz=cp*sa;
        ayx=axy;
        ayy=sp*sp*(ca-1)+1;
        ayz=sp*sa;
        azx=-axz;
        azy=-ayz;
        azz=ca;

        rx=R0*sa*cp;          // shift vector computation
        ry=R0*sa*sp;
        rz=R0*(ca-1)+D;       // z is additionally shifted by depth

        x=axx*fx+axy*fy+rx;   // coordinates on the tangent plane in the detector frame
        y=ayx*fx+ayy*fy+ry;
        z=azx*fx+azy*fy+rz;

        cp=cos(phi); sp=sin(phi);

        px=st*cp;             // direction in the CORSIKA frame
        py=st*sp;
        pz=ct;

        b=x*px+y*py+(z+R0-D)*pz;
        c=x*x+y*y+(z+R0-D)*(z+R0-D)-R0*R0;
        r=b-sqrt(b*b-c); // shift along (px,py,pz) to reach the surface from the tangent plane

        x-=px*r;              // final coordinates in the detector frame
        y-=py*r;
        z-=pz*r;
        t=r/(SPEED*1.e-11);

        th*=180/PI;
        phi*=180/PI;

        t0+=-log(MathModel::RandomDouble())*tT;
        tt=t0-(x*px+y*py+z*pz)/(SPEED*1.e-2);
        tts=floor(tt*1000)/1000;
        ttn=(tt-tts)*1.e9;
        if(f2k) cerr<<"EM "<<evt<<" 1 1970 0 "<<Output::f(tts)<<" "+Output::f(ttn)<<endl;

        Name=name[lept];
        eini=E;
        thini=th;

        if(GEN)	putOut(lept, 0, x, y, z, t, -1, E, NULL);
        else{
            eM->sett(-px, -py, -pz);

            if(mF>0 && (lept>=2 && lept<=5)){
                xn=eM->X(eM->ti, eM->tf)*mF;
                bool nu=(lept==2 || lept==4);
                double tI=(lept==5)?nT->dSdy(E):0;
                tI+=nT->dSdy(E, true, nu);
                tI+=nT->dSdy(E, false, nu);
                xn*=tI;
                if(xn>1) xn=1;
            }
            else xn=1;

            gens=0;
            elost=0;

            prop(lept, 0, E, x, y, z, t, eM->ti);
        }

        if(f2k){
            if(!GEN){
                mP->outUSER();
                tP->outUSER();
                cout<<"US almc_e "<<Name<<" "<<Output::f(eini)<<" "<<
                                   Output::f(elost)<<" "<<Output::f(thini)<<endl;
            }
            cout<<"EE"<<endl;
        }
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * general-purpose lepton propagator.
     */

    vector<PROPOSALParticle*> AtmFlux::propagate(PROPOSALParticle* p){

        vector<PROPOSALParticle*> I3p;
        if(f2k){
            I3p.resize(0);
            return I3p;
        }
        I3hist.clear();
        propagate(p->name, p->igen, p->gens, p->x*1.e-2, p->y*1.e-2, p->z*1.e-2, 180-p->theta, p->phi<180?p->phi+180:p->phi-180, p->r*1.e-2, p->e*1.e-3, p->t*1.e9);
        I3p.resize(I3hist.size());
        for(int i=0; i<I3hist.size(); i++) I3p[i]=(PROPOSALParticle *)I3hist.at(i);
        return I3p;
    }
    //----------------------------------------------------------------------------------------------------//

    /**
     * Main parser for the f2k file streams.
     */

//    void AtmFlux::mmcf2k(){
//        cerr<<" ---  *** Enter your input in F2000 format now ***  --- "<<endl;

//        try{
//            LineNumberReader file = new LineNumberReader(new InputStreamReader(Output.in));
//            StringTokenizer st;
//            String line, taux, prnt, type;
//            Vector buffer = new Vector(Output.HISTSIZE);
//            boolean tbegin=false;

//            int i, ipar, igen;
//            long events=0, tracks=0, vertices=0, missed;

//            file.readLine();

//            while((line=file.readLine())!=null){
//                st = new StringTokenizer(line);
//                if(!st.hasMoreTokens()) continue;
//                taux=st.nextToken();
//                if("TBEGIN".compare(taux)==0){
//                    beginAtmFlux();
//                    Output.out.println(line);
//                    tbegin=true;
//                    continue;
//                }
//                if(!tbegin){
//                    Output.out.println(line);
//                    continue;
//                }
//                buffer.addElement(line);
//                if("TR".compare(taux)==0){
//                    igen=Integer.parseInt(st.nextToken());
//                    if(igen>gens) gens=igen;
//                    prnt=st.nextToken();
//                    try{
//                        ipar=Integer.parseInt(prnt);
//                        if(ipar>gens) gens=ipar;
//                    }catch(NumberFormatException error){
//                        ipar=-1;
//                    }
//                }
//                else if("EE".compare(taux)==0 || "END".compare(taux)==0){
//                    mP.iniUSER();
//                    tP.iniUSER();

//                    for(i=0; i<buffer.size(); i++){
//                        line=(String)buffer.elementAt(i);
//                        st = new StringTokenizer(line);
//                        taux=st.nextToken();
//                        if("TR".compare(taux)==0){
//                            Output.out.println(line);

//                            igen=Integer.parseInt(st.nextToken());
//                            prnt=st.nextToken();
//                            type=st.nextToken();

//                            double x, y, z, th, phi, l, e, t;

//                            x=atof(st.nextToken());
//                            y=atof(st.nextToken());
//                            z=atof(st.nextToken());
//                            th=atof(st.nextToken());
//                            phi=atof(st.nextToken());
//                            try{
//                                l=atof(st.nextToken());
//                                if(l<0) l=0;
//                            }catch(NumberFormatException error){
//                                l=0;
//                            }
//                            e=atof(st.nextToken());
//                            t=atof(st.nextToken());

//                            Name=type;
//                            eini=e;
//                            thini=th;

//                            propagate(type, igen, gens, x, y, z, th, phi, l, e, t);
//                            tracks++;

//                        }
//                        else if("EM".compare(taux)==0){
//                            events++;
//                            Output.out.println(line);
//                        }
//                        else if("EE".compare(taux)==0){
//                            mP.outUSER();
//                            tP.outUSER();
//                            Output.out.println("US almc_e "+Name+" "+Output.f(eini)+" "+
//                                               Output.f(elost)+" "+Output.f(thini));
//                            Output.out.println(line);
//                        }
//                        else Output.out.println(line);
//                    }

//                    gens=0;
//                    buffer.clear();
//                }

//            }

//            if(events!=0) vertices=(mP.vertices+tP.vertices)/events;
//            missed=mP.missed+tP.missed;
//            Output.err.println("events read "+events+" tracks looped "+tracks+
//                               " average vertices "+vertices+" missed volume "+missed+
//                               " muon tracks "+mP.tracks+" tau tracks "+tP.tracks);
//        }catch(Exception error){
//            Output.err.println("Program finished with exception: "+error.toString());
//            throw new mmcException("input error");
//        }
//    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * general-purpose lepton propagator.
     */

     void AtmFlux::propagate(string type, int igen, int gens, double x, double y, double z, double th, double phi, double r, double E, double t){
        int i;
        double ct, st, cp, sp;
        this->gens=gens;
        this->th=th;
        this->phi=phi;

        for(i=0; i<13; i++)
        {
           // cout<<"In AtmFlux 964:\t"<<i<<"\t"<<type<<"\t"<<name[i]<<endl;
            if(type.compare(name[i])==0) break;
        } 
        if(i==13 || i==12 || i==8 || i==9) return;

        ct=cos(th*PI/180);
        st=sin(th*PI/180);
        cp=cos(phi*PI/180);
        sp=sin(phi*PI/180);

        px=st*cp;
        py=st*sp;
        pz=ct;

        x-=px*r;
        y-=py*r;
        z-=pz*r;
        t+=r/(SPEED*1.e-11);

        eM->sett(-px, -py, -pz);
        elost=0;

        xn=1;
        prop(i, igen, E, x, y, z, t, eM->r(x*1.e-3, y*1.e-3, z*1.e-3));
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * general-purpose lepton propagator.
     */

    void AtmFlux::prop(int lept, int igen, double E, double x, double y, double z, double t, double nf){

        if(nf>eM->tf){
            putOut(lept, igen, x, y, z, t, 0, E, NULL);
            return;
        }

        if(lept>=2 && lept<=7) while(true){
            if(E<nlow_*1.e-3){
                putOut(lept, igen, x, y, z, t, (eM->tf-nf)*1.e3, E, NULL);
                break;
            }
            int mlpt;
            bool nu, cc, el;
            double aux, l, tI, tIc, tIn, tIw=0, rnd, Xt, ml, ni;

            ni=nf;
            switch(lept){
            case 2:  nu=true;  el=false; mlpt=6;  break;
            case 3:  nu=false; el=false; mlpt=7;  break;
            case 4:  nu=true;  el=true;  mlpt=4;  break;
            case 5:  nu=false; el=true;  mlpt=5;  break;
            case 6:  nu=true;  el=false; mlpt=2;  break;
            case 7:  nu=false; el=false; mlpt=3;  break;
            default: nu=false; el=false; mlpt=12;
            }

            rnd=MathModel::RandomDouble();
            tIc=nT->dSdy(E, true, nu, rnd);
            tIn=nT->dSdy(E, false, nu, rnd);
            tIw=(lept==5)?nT->dSdy(E, rnd):0;
            tI=tIc+tIn+tIw;

            rnd=-xn*log(MathModel::RandomDouble())/(nF*tI);
            Xt=eM->X(max(ni, eM->ti), eM->tf, rnd);

            if(rnd>Xt) nf=eM->tf;
            else{
                nf=eM->t();
                if(abs(nf-eM->tf)<=abs(eM->tf-eM->ti)*COMPUTER_PRECISION) nf=eM->tf;  // computer precision control
            }
            l=(nf-ni)*1.e3;

            putOut(lept, igen, x, y, z, t, l, E, NULL);
            igen=gens;

            if(nf==eM->tf) break;

            x-=px*l;
            y-=py*l;
            z-=pz*l;
            t+=l/(SPEED*1.e-11);

            if(OSC) if(!el) if(MathModel::RandomDouble()<nT->Pmt(E, l)) lept=mlpt;

            switch(lept){
            case 2:  nu=true;  el=false; mlpt=0;  break;
            case 3:  nu=false; el=false; mlpt=1;  break;
            case 4:  nu=true;  el=true;  mlpt=8;  break;
            case 5:  nu=false; el=true;  mlpt=9;  break;
            case 6:  nu=true;  el=false; mlpt=10; break;
            case 7:  nu=false; el=false; mlpt=11; break;
            default: nu=false; el=false; mlpt=12;
            }
           // cout<<"ATMFLUX 1060: Neutrino!!!!!!!!!!!!!!!"<<endl;
            rnd=MathModel::RandomDouble()*tI;
            if(rnd<tIc){
                aux=nT->e(true, nu);
                E-=aux;
                putOut(12, igen, x, y, z, t, 0, aux, NULL);
                if(insDet(x, y, z)) elost+=aux;
                prop(mlpt, igen, E, x, y, z, t, nf);
                break;
            }
            else if(rnd<tIc+tIn){
                aux=nT->e(false, nu);
                putOut(12, igen, x, y, z, t, 0, aux, NULL);
                if(insDet(x, y, z)) elost+=aux;
                E-=aux;
            }
            else{
                aux=MathModel::RandomDouble();
                if(aux<1/9.){
                    aux=nT->e(ME);
                    putOut(8, igen, x, y, z, t, 0, aux, NULL);
                    if(insDet(x, y, z)) elost+=aux;
                    E-=aux;
                    lept=5;
                }
                if(aux<2/9.){
                    aux=nT->e(MMU);
                    prop(0, igen, aux, x, y, z, t, nf);
                    E-=aux;
                    lept=3;
                }
                if(aux<3/9.){
                    aux=nT->e(MTAU);
                    prop(10, igen, aux, x, y, z, t, nf);
                    E-=aux;
                    lept=7;
                }
                else{
                    putOut(12, igen, x, y, z, t, 0, E, NULL);
                    if(insDet(x, y, z)) elost+=aux;
                    break;
                }
            }
        }
        else if(lept==0 || lept==1 || lept==10 || lept==11){
            int nlpt=0;
            Amanda *aP;
            vector<PROPOSALParticle*> p;
            PROPOSALParticle *pi = new PROPOSALParticle(gens+1, gens+1, name[lept], x*1.e2, y*1.e2, z*1.e2, 180-th, phi<180?phi+180:phi-180, E*1.e3, t*1.e-9, 0);
            aP=lept<2?mP:tP;
            if(aP->medianum>1) aP->rho[aP->medianum-1]=eM->intRho(nf)/aP->srho[aP->medianum-1];
            p=aP->propagate(pi); if(f2k) aP->recUSER(); elost+=aP->elost;
            putOut(lept, igen, x, y, z, t, pi->r*1.e-2, E, pi);
            igen=gens;
            for(int i=0; i<p.size(); i++){
              
		bool flag=true;
                if(p.at(i)->name.compare("nu_e")==0) nlpt=4;
                else if(p.at(i)->name.compare("~nu_e")==0) nlpt=5;
                else if(p.at(i)->name.compare("nu_mu")==0) nlpt=2;
                else if(p.at(i)->name.compare("~nu_mu")==0) nlpt=3;
                else if(p.at(i)->name.compare("nu_tau")==0) nlpt=6;
                else if(p.at(i)->name.compare("~nu_tau")==0) nlpt=7;
                else if(p.at(i)->name.compare("mu")==0 || p.at(i)->name.compare("mu-")==0) nlpt=0;
                else if(p.at(i)->name.compare("mu+")==0) nlpt=1;
                else if(p.at(i)->name.compare("tau")==0 || p.at(i)->name.compare("tau-")==0) nlpt=10;
                else if(p.at(i)->name.compare("tau+")==0) nlpt=11;
                else flag=false;
                if(flag){
                    double ei=p.at(i)->e*1.e-3, xi=p.at(i)->x*1.e-2, yi=p.at(i)->y*1.e-2, zi=p.at(i)->z*1.e-2;
                    prop(nlpt, igen, ei, xi, yi, zi, p.at(i)->t*1.e9, nf+pi->r*1.e-5);
                    if(lept>=2 && lept<=7) if(insDet(xi, yi, zi)) elost-=ei;
                }
                else lptOut(p.at(i));
            }
        }
        else{
            putOut(lept, igen, x, y, z, t, 0, E, NULL);
            if(insDet(x, y, z)) elost+=E;
        }
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates if the coordinates are inside the detector
     */

    bool AtmFlux::insDet(double x, double y, double z){
        switch(mP->gdet){
        case 1:
            return abs(x)<mP->length/2 && abs(y)<mP->width/2 && abs(z)<mP->height/2;
        case 2:
            return x*x+y*y+z*z<mP->radius*mP->radius;
        case 0:
        default:
            return x*x+y*y<mP->radius*mP->radius && -mP->length/2<z && z<mP->length/2;
        }
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Constructs the lepton from the given information.
     */

    void AtmFlux::putOut(int lept, int igen, double x, double y, double z, double t, double l, double E, PROPOSALParticle* p){
        gens++;
        if(f2k){
            if(gens>0){
                cout<<"TR "<<gens<<" "<<igen<<" "<<name[lept]<<" "<<Output::f(x)<<" "<<Output::f(y)<<" "<<Output::f(z)
                              <<" "<<Output::f(th)<<" "<<Output::f(phi)<<" "<<(l>=0?Output::f(l):"?")<<" "<<Output::f(E)<<" "<<Output::f(t)<<endl;}
            else{
                cout<<"TR "<<gens<<" ? "<<name[lept]<<" "<<Output::f(x)<<" "<<Output::f(y)<<" "<<Output::f(z)
                    <<" "<<Output::f(th)<<" "<<Output::f(phi)<<" "<<(l>=0?Output::f(l):"?")<<" "<<Output::f(E)<<" "<<Output::f(t)<<endl;
            }
        }
        else I3hist.push_back(new PROPOSALParticle(igen, gens, name[lept], x*1.e2, y*1.e2, z*1.e2, 180-th, phi<180?phi+180:phi-180, E*1.e3, t*1.e-9, l*1.e2, p));
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Copies the lepton into the output stream.
     */

    void AtmFlux::lptOut(PROPOSALParticle *p){
        gens++;
        if(f2k) cout<<"TR "<<p->gens<<" "<<p->igen<<" "<<p->name<<" "<<
                                   Output::f(p->x*1.e-2)<<" "<<Output::f(p->y*1.e-2)<<" "<<Output::f(p->z*1.e-2)<<" "<<
                                   Output::f(180-p->theta)<<" "<<Output::f(p->phi<180?p->phi+180:p->phi-180)<<" "<<
                                   Output::f(p->r*1.e-2)<<" "<<Output::f(p->e*1.e-3)<<" "<<Output::f(p->t*1.e9)<<endl;
        else I3hist.push_back(p);
    }

    //----------------------------------------------------------------------------------------------------//

    double AtmFlux::getTotFlux(double x){
        return getTotFlux(x, -1);
    }

    //----------------------------------------------------------------------------------------------------//

    double AtmFlux::getTotFlux(double x, double rnd){
        if(rnd<0){
            if(jt) return J.at(lept)->interpolate(x);
            return I->integrateOpened(0, x, this);
        }
        else{
            if(jt) return min(max(J.at(lept)->findLimit(rnd*J.at(lept)->interpolate(x)), 0.), x);
            I->integrateOpened(0, x, this, rnd);
            return I->getUpperLimit();
        }
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Function for flux integration over zenith angles - interface to Integral.
     */

    double AtmFlux::function(double x){
        return 2*PI*iF.at(lept)->getIntFlux(Ecut[lept], x);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Parametrizes the integral of this class.
     */

    void AtmFlux::interpolate(){
        int g=4;  // do not change these settings
        jt=false;
        cerr<<"Parameterizing totflX ... "<<endl;
        J.resize(6);
        for(int i=0; i<6; i++){
            lept=i;
            J.at(i) = new Interpolate(NUM3, 0, 1, this, g, true, false, false, g, false, false, true);
        }
        jt=true;
        cerr<<"done"<<endl;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * 1d parametrization - interface to Interpolate
     */

    double AtmFlux::functionInt(double x){
        return getTotFlux(x);
    }


