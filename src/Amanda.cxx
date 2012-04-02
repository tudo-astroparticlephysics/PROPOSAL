#include "Amanda.h"
#include <cmath>
#include <algorithm>
#include "Output.h"
#include "methods.h"
#include "Bremsstrahlung.h"
#include "Photonuclear.h"
#include "Medium.h"



using namespace std;


bool Amanda::lfix =true;
double Amanda::LENGTH = 800;
double Amanda::RADIUS = 400;
double Amanda::WIDTH = 800;
double Amanda::HEIGHT = 800;
double Amanda::	BIG_ = 6.40e7;

//----------------------------------------------------------------------------------------------------//




Amanda::Amanda(){


    SURF=false;
    FACE=false;
    USER=false;
    USFI=false;
    SDEC=false;
    RECC=false;
    zset=0;


    mediamax=100;
    mediadef="";            // = NULL !!?!

    vcut[0] =  -1.0;
    vcut[1] =  -1.0;
    ecut[0] =  -1.0;
    ecut[1] =  -1.0;
    vaux=-1.0;
    eaux=-1.0;
    med="Ice";
    muta="mu";
    usna="mmc_en";
    elow=PhysicsModel::get_elow();
    ebig=PhysicsModel::get_ebig();

    timef=false;
    conti[0]=false;
    conti[1]=false;
    lpmef=false;
    scatt=false;
    frho=false;
    rfix=false;
    amasim=false;

    bspar=1;
    pncrs=1;
    pncbb=1;
    pncsh=1;

    crsci=1;
    crscb=1;
    crscp=1;
    crsce=1;
    crscd=1;
    drho=1;

    bspar=1;
    pncrs=1;
    pncbb=1;
    pncsh=1;

    romb=5;
    seed=0;
    SEED = false;
    raw=false;
    tdir = "";
    surf= 0;
    emax=0;

    gens = 0;
    imax = 0;

    param="Amanda";
    flag =0;
    dtct="";
    gdet=0;
    length=LENGTH, radius=RADIUS, width=WIDTH, height=HEIGHT;

        dw=false;
    rw=0, hw=0, fw=1;
    events=0, tracks=0, vertices=0, missed=0;

        medianum=1;
    DEBUG=false;
    z1=0;
    z2=0;
    h1=0;
    h2=0;
    nx=0;
    ny=0;
    nz=0;
    e1=0;
    e2=0;
    ec=0;
    elost=0;

    //Output::I3flag=false;

}
/**
 * Calculate points of intersection with the cylinder.
 */

    void Amanda::setcyl(double x, double y, double z, double cosph, double sinph, double costh,double sinth, double length, double radius){


        double aux, r1, r2;
        double b = x*cosph+y*sinph;
        double d = b*b+radius*radius-x*x-y*y;
        if(d>0){
            d=sqrt(d);
            if(costh!=0)
            {
                hx1=(z-length/2)/costh;
                hx2=(z+length/2)/costh;
                if(hx1>hx2)
                {
                    aux=hx1;
                    hx1=hx2;
                    hx2=aux;
                }
            }
            if(sinth!=0)
            {
                r1=(b-d)/sinth;
                r2=(b+d)/sinth;
                if(r1>r2)
                {
                    aux=r1;
                    r1=r2;
                    r2=aux;
                }
                if(costh==0)
                {
                    if(z>-length/2 && z<length/2)
                    {
                        hx1=r1;
                        hx2=r2;
                    }
                    else
                    {
                        hx1=0;
                        hx2=0;
                    }
                }
                else
                {
                    if(hx1>=r2 || hx2<=r1)
                    {
                        hx1=0;
                        hx2=0;
                    }
                    else
                    {
                        hx1=max(r1, hx1);
                        hx2=min(r2, hx2);
                    }
                }
            }
        }
        else{ hx1=0; hx2=0; }
        //cout<<"In Amanda::setcyl:\t hx1="<<hx1<<"\t hx2="<<hx2<<endl;
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculate points of intersection with the box.
     */

    void Amanda::setbox(double x, double y, double z, double dx, double dy, double dz, double length, double width, double height){

        bool top=false, bottom=false, tbflag=false;
        double hx;
        if(dx!=0)
        {
            if(!tbflag)
            {
                hx=-(x-length/2)/dx;
                if(fabs(y+hx*dy)<=width/2 && fabs(z+hx*dz)<=height/2)
                {
                    if(!top)
                    {
                        hx1=hx;
                        top=true;
                    }
                    else if(hx!=hx1 && !bottom)
                    {
                        hx2=hx;
                        bottom=true;
                    }
                    else tbflag=true;
                }
            }
            if(!tbflag)
            {
                hx=-(x+length/2)/dx;
                if(fabs(y+hx*dy)<=width/2 && fabs(z+hx*dz)<=height/2)
                {
                    if(!top)
                    {
                        hx1=hx;
                        top=true;
                    }
                    else if(hx!=hx1 && !bottom)
                    {
                        hx2=hx;
                        bottom=true;
                    }
                    else tbflag=true;
                }
            }
        }
        if(dy!=0)
        {
            if(!tbflag)
            {
                hx=-(y-width/2)/dy;
                if(fabs(x+hx*dx)<=length/2 && fabs(z+hx*dz)<=height/2)
                {
                    if(!top)
                    {
                        hx1=hx;
                        top=true;
                    }
                    else if(hx!=hx1 && !bottom)
                    {
                        hx2=hx;
                        bottom=true;
                    }
                    else tbflag=true;
                    }
            }
            if(!tbflag)
            {
                hx=-(y+width/2)/dy;
                if(fabs(x+hx*dx)<=length/2 && fabs(z+hx*dz)<=height/2)
                {
                    if(!top)
                    {
                        hx1=hx;
                        top=true;
                    }
                    else if(hx!=hx1 && !bottom)
                    {
                        hx2=hx;
                        bottom=true;
                    }
                    else tbflag=true;
                }
            }
        }
        if(dz!=0)
        {
            if(!tbflag)
            {
                hx=-(z-height/2)/dz;
                if(fabs(x+hx*dx)<=length/2 && fabs(y+hx*dy)<=width/2)
                {
                    if(!top)
                    {
                        hx1=hx;
                        top=true;
                    }
                    else if(hx!=hx1 && !bottom)
                    {
                        hx2=hx;
                        bottom=true;
                    }
                    else tbflag=true;
                }
            }
            if(!tbflag)
            {
                hx=-(z+height/2)/dz;
                if(fabs(x+hx*dx)<=length/2 && fabs(y+hx*dy)<=width/2)
                {
                    if(!top)
                    {
                        hx1=hx;
                        top=true;
                    }
                    else if(hx!=hx1 && !bottom)
                    {
                        hx2=hx;
                        bottom=true;
                    }
                    else tbflag=true;
                }
            }
        }
        if(!top)
        {
            hx1=0;
        }
        if(!bottom)
        {
            hx2=0;
        }
        if(hx1>hx2)
        {
            hx=hx1;
            hx1=hx2;
            hx2=hx;
        }
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculate points of intersection with the sphere.
     */

    void Amanda::setsph(double x, double y, double z, double cosph, double sinph, double costh, double sinth, double radius){
        double b = (x*cosph+y*sinph)*sinth+(z+radius)*costh;
        double d = b*b-(x*x+y*y+z*(z+2*radius));
        if(d>0)
        {
            d=sqrt(d);
            hx1=b-d;
            hx2=b+d;
        }
        else
        {
            hx1=0;
            hx2=0;
        }
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculate points of intersection with the box.
     */

    void Amanda::setpln(double x, double y, double z, double dx, double dy, double dz, double nx, double ny, double nz){
        double b = -(x*nx+y*ny+z*nz);
        double d = dx*nx+dy*ny+dz*nz;
        if(b>=0 || b*d>0)
        {
            if(b>=0)
            {
                hx1=0;
                        //hx2=d!=0?b/d:BIG_;
                if(d!=0)
                {
                    hx2=b/d;
                }
                else
                {
                    hx2=BIG_;
                }
            }
            else
            {
                       //hx1=d!=0?b/d:0;
                if(d!=0)
                {
                    hx1=b/d;
                }
                else
                {
                    hx1=0;
                }

                hx2=BIG_;
            }
        }
        else
        {
            hx1=0;
            hx2=0;
        }
    }




    //----------------------------------------------------------------------------------------------------//

    /**
     * checks the flag, corrects the units, initializes the particle and calls the propagator
     */

    double Amanda::propagateTo(double h, double e, double n, int i, Propagate *p, int igen, int gens,
                                      double time, double x, double y, double z, double theta, double phi){

        double result;
        switch(i)
        {
            case 1: p->get_output()->initDefault(igen, gens, type, time, x, y, z, theta, phi); break;
            case 2: p->get_output()->initF2000(igen, gens, type, time, x, y, z, theta, phi); break;
            case 3: p->get_output()->initDefault(igen, gens, type, time, x, y, z, theta, phi); break;
        }
        if(flag==i && n>=0 && n<h)
        {
            result = p->propagateTo(h*1.e2, e*1.e3, n*1.e2);
            ec = p->getPropEc();

            if(ec>0)
            {
                ec = ec*1.e-3;
            }
            else
            {
                ec = ec*1.e-2;
            }
            if(Output::I3flag)
            {
                pI->tc = p->getPropTc();
            }
        }
        else
        {
            result = p->propagateTo(h*1.e2, e*1.e3);
        }


                    //return result>0?result*1.e-3:result*1.e-2;
        if(result>0)
        {
            return result*1.e-3;
        }
        else
        {
            return result*1.e-2;
        }


    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Front-end for AMANDA, reads from and outputs to standard input/output in F2000 format.
     */
    // string[] -> string
     void Amanda::main(string args){
//        Amanda *A = new Amanda();
//        A->mmcf2k(args);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Initializes propagation for external applications, e.g.&nbsp;mmc-icetray (through jni).
     */

     void Amanda::setup(string args){
        dtct=".icecube";
        setup(Output::splitString(args, " \t"));
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Initializes propagation for external applications, e.g.&nbsp;mmc-icetray (through jni).
     */
// string[] -> string
     void Amanda::setup(deque<string> args){
        Output::I3flag=true;
        setp(args);
        I3hist.resize(Output::HISTSIZE);

    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Propagates particles for the external applications. Returns an array of secondaries.
     * If "-user" option is used, fills user-block variables.
     */
vector<PROPOSALParticle*> Amanda::propagate(PROPOSALParticle *p){


        //In Java NULL ist returned...
        vector<PROPOSALParticle*> I3p;

        if(!Output::I3flag)
        {
            I3p.resize(0);
            I3p.at(0) = NULL;
            return I3p;
        }

        pI=p;
        type=p->name;


        if(muta.compare(type)==0 || (muta+"-").compare(type)==0 || (muta+"+").compare(type)==0)
        {
            if(p->phi<180)
            {
                p->r=prop(p->igen, p->gens, p->x*1.e-2, p->y*1.e-2, p->z*1.e-2, 180-p->theta, p->phi+180, p->r*1.e-2, p->e*1.e-3, p->t*1.e9);
            }
            else
            {
                p->r=prop(p->igen, p->gens, p->x*1.e-2, p->y*1.e-2, p->z*1.e-2, 180-p->theta, p->phi-180, p->r*1.e-2, p->e*1.e-3, p->t*1.e9);

            }
            I3p.resize(I3hist.size());
            for(int i=0; i<(int)I3hist.size(); i++)
            {
                I3p[i]=(PROPOSALParticle*)I3hist.at(i);
            }
        }
        else{
            I3p.resize(1);
            I3p.at(0) = new PROPOSALParticle();
        }
        return I3p;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * This is the command-line option parser. Call with "-help" to list all options.
     */

     bool Amanda::setp(deque<string> args){
        string interpolString="";
        int i, j;

        int bnum=0;
        string pbad="";
        bool pflag;

        rho = (double *)calloc(mediamax,sizeof(double));
        srho = (double *)calloc(mediamax,sizeof(double));
        double vcut[mediamax][mediamax];// = {(double *)calloc(mediamax,sizeof(double)) , (double *)calloc(mediamax,sizeof(double))};
        double ecut[mediamax][mediamax];// = {(double *)calloc(mediamax,sizeof(double)) , (double *)calloc(mediamax,sizeof(double))};
        vector<string> med;
        med.resize(mediamax);
        bool conti[mediamax][mediamax];// = {(bool *)calloc(mediamax,sizeof(bool)) , (bool *)calloc(mediamax,sizeof(bool))};

        rho[0]=this->drho;
        vcut[0][0]=this->vcut[0];
        vcut[1][0]=this->vcut[1];
        ecut[0][0]=this->ecut[0];
        ecut[1][0]=this->ecut[1];
        med[0]=this->med;
        conti[0][0]=this->conti[0];
        conti[1][0]=this->conti[1];

        stringstream mediaopts;
        mediaopts.str("");
        mediaopts.clear();
        int mediaseg=1;
        medt = (int *)calloc(mediamax,sizeof(int));
        sphz = (double *)calloc(mediamax,sizeof(double));
        sphr = (double *)calloc(mediamax,sizeof(double));
        boxx = (double *)calloc(mediamax,sizeof(double));
        boxy = (double *)calloc(mediamax,sizeof(double));
        boxz = (double *)calloc(mediamax,sizeof(double));
        boxl = (double *)calloc(mediamax,sizeof(double));
        boxw = (double *)calloc(mediamax,sizeof(double));
        boxh = (double *)calloc(mediamax,sizeof(double));
        cylz = (double *)calloc(mediamax,sizeof(double));
        cylr = (double *)calloc(mediamax,sizeof(double));
        cyll = (double *)calloc(mediamax,sizeof(double));
        bool *medn =(bool *)calloc(mediamax,sizeof(bool));

        medianum=1;
        medt[0]=0;
        medn[0]=true;

        for(int n=0; n<(int)args.size(); n++){
            pflag=true;
            if(args.at(n).compare("-help")==0 || args.at(n).compare("-h")==0 || args.at(n).compare("--help")==0){
                cout<<"\n"<<
"This program propagates muons in "<<med[0]<<" to/through the detector\n"<<
"Available options are: -length=[LENGTH of the detector volume in meters]\n"<<
"                       -radius=[RADIUS of the detector volume in meters]\n"<<
"                       -width=[WIDTH of the detector volume in meters]\n"<<
"                       -height=[HEIGHT of the detector volume in meters]\n"<<
"                       -gdet=[0-2] detector is a cylinder/box/sphere\n"<<
"                       -vcut=[value of vcut used for the 1st region]\n"<<
"                       -ecut=[ecut in  MeV  used for the 2nd region]\n"<<
"                       -medi=[medium name]\n"<<
"                       -mediadef=[file with media definitions]\n"<<
"                       -tau  propagate taus instead of muons\n"<<
"                       -e    propagate electrons instead of muons\n"<<
"                       -monopole[=mass in GeV] propagate monopoles\n"<<
"                       -stau[=mass in GeV]     propagate staus\n"<<
"                       -sdec enable stopped muon decay treatment\n"<<
"                       -recc enable printout of continuous energy losses\n"<<
"                       -user     enable the mmc_en user line\n"<<
"                       -user=[z] same, but record energy at z, not CPD\n"<<
"                       -rdmc     enforce compliance with rdmc\n"<<
"                       -amasim   turn on workarounds for amasim\n"<<
"                       -time precise time of flight calculation\n"<<
"                       -cont enable continuous loss randomization\n"<<
"                       -scat enable Moliere scattering\n"<<
"                       -lpm  enable lpm treatment\n"<<
"                       -bs=[1-4]  bremsstrahlung: kkp, abb, ps, csc\n"<<
"                       -ph=[1-4]  photonuclear: bb, bb+bs, allm, bm\n"<<
"                       -bb=[bb/bs:1-4 3|4:bb|zeus, allm:1-2(91/7), bm:1]\n"<<
"                       -sh=[1-2]  nuclear structure function: dutt/butk\n"<<
"                       -c[i/b/p/e/d]=[cross section modifier, 0:disable]\n"<<
"                       -rho=[multiplicative factor] for medium density\n"<<
"                       -frho  enable smart density factor handling\n"<<
"                       -rw=[reweight cross sections by 1-x^+rw or x^-rw]\n"<<
"                       -elow=[muon energy in GeV below which it is lost]\n"<<
"                       -ebig=[upper bound in GeV of the paramet. tables]\n"<<
"                       -surf=[h in meters]  propagate to the plane z=[h]\n"<<
"                       -face  only if detector is on opposite side of it\n"<<
"                       -romb=[number of interpolation points]\n"<<
"                       -seed=[integer] sets random number generator seed\n"<<
"                       -raw  save tables in raw format\n"<<
"                       -tdir=[dir] specify directory for paramet. tables\n" <<
"                       -interpol=[string] default 'all'. Used option 'off'. More information in Propagate.h\n";
                return false;
            }
            else if(args.at(n).compare("-tau")==0){
                muta="tau";
                usna="mmc_et";
            }
            else if(args.at(n).compare("-e")==0){
                muta="e";
                usna="mmc_el";
            }
            else if(args.at(n).compare("-sdec")==0){
                SDEC=true;
            }
            else if(args.at(n).compare("-recc")==0){
                RECC=true;
            }
            else if(args.at(n).compare("-user")==0){
                USER=true;
            }
            else if(args.at(n).compare("-rdmc")==0){
                rfix=true;
            }
            else if(args.at(n).compare("-amasim")==0){
                amasim=true;
            }
            else if(args.at(n).compare("-time")==0){
                timef=true;
            }
            else if(args.at(n).compare("-cont")==0){
                conti[0][0]=true;
            }
            else if(args.at(n).compare("-scat")==0){
                scatt=true;
            }
            else if(args.at(n).compare("-lpm")==0){
                lpmef=true;
            }
            else if(args.at(n).compare("-frho")==0){
                frho=true;
            }
            else if(args.at(n).compare("-face")==0){
                FACE=true;
            }
            else if(args.at(n).compare("-raw")==0){
                raw=true;
            }
            else if(StartsWith(args.at(n),"-monopole")){
                muta=args.at(n).substr(1);
                usna="mmc_mn";
            }
            else if(StartsWith(args.at(n),"-stau")){
                muta=args.at(n).substr(1);
                usna="mmc_st";
            }
            else if(StartsWith(args.at(n),"-length="))
            {
                    length=atof(args.at(n).substr(8).c_str());
            }
            else if(StartsWith(args.at(n),"-interpol="))
            {
                    interpolString=args.at(n).substr(10).c_str();
                    while( (n+1)<(int)args.size() && args.at(++n).find("-")!=0){
                        interpolString += " ";
                        interpolString += args.at(n).c_str();
                    }

            }
            else if(StartsWith(args.at(n),"-radius="))
            {
                    radius=atof(args.at(n).substr(8).c_str());
            }
            else if(StartsWith(args.at(n),"-width="))
            {
                    width=atof(args.at(n).substr(7).c_str());
            }
            else if(StartsWith(args.at(n),"-height="))
            {
                    height=atof(args.at(n).substr(8).c_str());
            }
            else if(StartsWith(args.at(n),"-surf="))
            {
                    surf=atof(args.at(n).substr(6).c_str());
                    SURF=true;
            }
            else if(StartsWith(args.at(n),"-vcut="))
            {
                    double aux=atof(args.at(n).substr(6).c_str());
                    if(aux<0) ecut[0][0]=-aux;
                    else vcut[0][0]=aux;
            }
            else if(StartsWith(args.at(n),"-ecut="))
            {
                    double aux=atof(args.at(n).substr(6).c_str());
                    if(aux<0) vcut[1][0]=-aux;
                    else ecut[1][0]=aux;
            }
            else if(StartsWith(args.at(n),"-user="))
            {
                    USER=true;
                    USFI=true;
                    zset=atof(args.at(n).substr(6).c_str());
            }
            else if(StartsWith(args.at(n),"-gdet="))
            {
                    gdet=(int)atof(args.at(n).substr(6).c_str());
            }
            else if(StartsWith(args.at(n),"-bs="))
            {
                    bspar=(int)atof(args.at(n).substr(4).c_str());
            }
            else if(StartsWith(args.at(n),"-ph="))
            {
                    pncrs=(int)atof(args.at(n).substr(4).c_str());
            }
            else if(StartsWith(args.at(n),"-bb="))
            {
                    pncbb=(int)atof(args.at(n).substr(4).c_str());
            }
            else if(StartsWith(args.at(n),"-sh="))
            {
                    pncsh=(int)atof(args.at(n).substr(4).c_str());
            }
            else if(StartsWith(args.at(n),"-ci="))
            {
                    crsci=atof(args.at(n).substr(4).c_str());
            }
            else if(StartsWith(args.at(n),"-cb="))
            {
                    crscb=atof(args.at(n).substr(4).c_str());
            }
            else if(StartsWith(args.at(n),"-cp="))
            {
                    crscp=atof(args.at(n).substr(4).c_str());
            }
            else if(StartsWith(args.at(n),"-ce="))
            {
                    crsce=atof(args.at(n).substr(4).c_str());
            }
            else if(StartsWith(args.at(n),"-cd="))
            {
                    crscd=atof(args.at(n).substr(4).c_str());
            }
            else if(StartsWith(args.at(n),"-rho="))
            {
                    rho[0]=atof(args.at(n).substr(5).c_str());
                    if(rho[0]==0)
                    {
                        rho[0]=this->drho;
                    }

            }
            else if(StartsWith(args.at(n),"-rw="))
            {
                    rw=atof(args.at(n).substr(4).c_str());
            }
            else if(StartsWith(args.at(n),"-elow="))
            {
                    elow=atof(args.at(n).substr(6).c_str())*1.e3;
                    if(elow==0)
                    {
                        elow=PhysicsModel::get_elow();
                    }
            }
            else if(StartsWith(args.at(n),"-ebig=")){

                    ebig=atof(args.at(n).substr(6).c_str())*1.e3;
                    if(ebig==0)
                    {
                        ebig=PhysicsModel::get_ebig();
                    }
            }
            else if(StartsWith(args.at(n),"-romb="))
            {
                    romb=(int)atof(args.at(n).substr(6).c_str());
                    if(romb==0)
                    {
                        romb=5;
                    }
            }
            else if(StartsWith(args.at(n),"-seed="))
            {
                    seed=atol(args.at(n).substr(6).c_str());
                    SEED=true;
            }
            else if(StartsWith(args.at(n),"-medi="))
            {
                //--------------->  med[0]=args.at(n).substr(6).replace('-', ' ');
                med[0]=args.at(n).substr(6).c_str();
            }
            else if(StartsWith(args.at(n),"-mediadef="))
            {
                mediadef=args.at(n).substr(10).c_str();
            }
            else if(StartsWith(args.at(n),"-tdir="))
            {
                tdir=args.at(n).substr(6)+"/";
            }
            else
            {
                bnum++;
                if(bnum>1)
                {
                    pbad+=",";
                    pbad+=" \"";
                    pbad+=args.at(n);
                    pbad+="\"";
                }
                else
                {
                    pbad+=" ";
                    pbad+=" \"";
                    pbad+=args.at(n);
                    pbad+="\"";
                }

                pflag=false;
            }
            if(pflag)
            {
                param+=" "+args.at(n);
            }
        }

        cerr<<Output::version<<endl;
        cerr<<"Running \""<<param<<"\""<<endl;
        if(bnum>0)
        {
            param+=" (not used:"+pbad+")";
            pbad=bnum==1?pbad+" is":"s"+pbad+" are";
            cerr<<"Warning: Parameter"<<pbad<<" not recognized"<<endl;
        }

        if(gdet!=0 && gdet!=1 && gdet!=2){
            cerr<<"Warning: gdet is not a valid number"<<endl;
            gdet=0;
        }
        if(length<=0){
            cerr<<"Warning: length is not a valid number"<<endl;
            length=LENGTH;
        }
        if(radius<=0){
            cerr<<"Warning: radius is not a valid number"<<endl;
            radius=RADIUS;
        }
        if(width<=0){
            cerr<<"Warning: width is not a valid number"<<endl;
            width=WIDTH;
        }
        if(height<=0){
            cerr<<"Warning: height is not a valid number"<<endl;
            height=HEIGHT;
        }
        if(!SURF || surf==0) FACE=false;
        if(bspar<1 || bspar>4){
            cerr<<"Warning: bs is not a valid number"<<endl;
            bspar=1;
        }
        if(pncrs<1 || pncrs>4 || (pncrs==2 && !(muta.compare("mu")==0 || muta.compare("tau")==0))){
            cerr<<"Warning: ph is not a valid number"<<endl;
            pncrs=1;
        }
        if(((pncrs==1 || pncrs==2) && (pncbb<1 || pncbb>4)) ||
            (pncrs==3 && (pncbb<1 || pncbb>2)) || (pncrs==4 && pncbb!=1)){
            cerr<<"Warning: bb is not a valid number"<<endl;
            pncbb=1;
        }
        if(((pncrs==1 || pncrs==2) && (pncsh!=1)) || ((pncrs>2) && (pncsh<1 || pncsh>2))){
            cerr<<"Warning: sh is not a valid number"<<endl;
            pncsh=1;
        }
        if(romb<2 || romb>6){
            cerr<<"Warning: romb is not a valid number"<<endl;
            romb=5;
        }

        if(&mediadef!=NULL){
            if(Output::FileExist(mediadef))
            {
//                Reader ifile = Output.reader(mediadef);
//                LineNumberReader inf = new LineNumberReader(ifile);
                ifstream mediadef_file;
                mediadef_file.open(mediadef.c_str());
                int mediacur;
                char buf[256];
                string str_buf;
                deque<string> *token;
                mediaopts<<"HI  Attaching the media definition file:\n";

                while(mediadef_file.good()){
                    mediadef_file.getline(buf,256);
                    str_buf = string(buf);
                    mediaopts<<"HI  ! "<<str_buf<<"\n";
                    token =splitString(str_buf," \t");

                    if(token->empty()) continue;
                    string taux=nextToken(token);


                    if(taux.at(0)=='#')
                    {
                        continue;
                    }
                    if(medianum==mediamax){
                        cerr<<"Warning: Number of defined media will not exceed "<<mediamax<<endl;
                        break;
                    }
                    bool sflag=false;
                    if((toLowerCase(taux)).compare("all")==0)
                    {
                        mediacur=0;
                        medt[mediacur]=0;
                    }
                    else if((toLowerCase(taux)).compare("sphere")==0)
                    {
                        mediacur=medianum;
                        medt[mediacur]=1;
                        sphz[mediacur]=atof(nextToken(token).c_str());
                        sphr[mediacur]=atof(nextToken(token).c_str());
                        if(sphr[mediacur]<0)
                        {
                            sphr[mediacur]*=-1;
                            sflag=true;
                        }
                        if(sflag)
                        {
                            medt[mediacur]*=-1;
                        }
                        for(i=0; i<medianum; i++)
                        {
                            if(medt[i]==medt[mediacur])
                            {
                                if(sphz[i]==sphz[mediacur] && sphr[i]==sphr[mediacur])
                                {
                                    break;
                                }
                            }
                        }
                        if(i==medianum)
                        {
                            medianum++;
                            if(sflag)
                            {
                                mediaseg+=4;
                            }
                            else
                            {
                                mediaseg+=2;
                            }
                        }
                    }
                    else if((toLowerCase(taux)).compare("box")==0){
                        mediacur=medianum;
                        medt[mediacur]=2;
                        boxx[mediacur]=atof(nextToken(token).c_str());
                        boxy[mediacur]=atof(nextToken(token).c_str());
                        boxz[mediacur]=atof(nextToken(token).c_str());
                        boxl[mediacur]=atof(nextToken(token).c_str());
                        if(boxl[mediacur]<0)
                        {
                            boxl[mediacur]*=-1;
                            sflag=true;
                        }
                        boxw[mediacur]=atof(nextToken(token).c_str());
                        if(boxw[mediacur]<0)
                        {
                            boxw[mediacur]*=-1;
                            sflag=true;
                        }
                        boxh[mediacur]=atof(nextToken(token).c_str());
                        if(boxh[mediacur]<0)
                        {
                            boxh[mediacur]*=-1;
                            sflag=true;
                        }
                        if(sflag)
                        {
                            medt[mediacur]*=-1;
                        }
                        for(i=0; i<medianum; i++)
                        {
                            if(medt[i]==medt[mediacur])
                            {
                                if(boxx[i]==boxx[mediacur] && boxy[i]==boxy[mediacur] && boxz[i]==boxz[mediacur] &&
                                   boxl[i]==boxl[mediacur] && boxw[i]==boxw[mediacur] && boxh[i]==boxh[mediacur])
                                {
                                    break;
                                }
                            }
                        }
                        if(i==medianum)
                        {
                            medianum++;
                            if(sflag)
                            {
                                mediaseg+=4;
                            }
                            else
                            {
                                mediaseg+=2;
                            }

                        }
                    }
                    else if((toLowerCase(taux)).compare("cyl")==0)
                    {
                        mediacur=medianum;
                        medt[mediacur]=3;
                        cylz[mediacur]=atof(nextToken(token).c_str());
                        cylr[mediacur]=atof(nextToken(token).c_str());
                        if(cylr[mediacur]<0)
                        {
                            cylr[mediacur]*=-1;
                            sflag=true;
                        }
                        cyll[mediacur]=atof(nextToken(token).c_str());
                        if(cyll[mediacur]<0)
                        {
                            cyll[mediacur]*=-1;
                            sflag=true;
                        }
                        if(sflag)
                        {
                            medt[mediacur]*=-1;
                        }
                        for(i=0; i<medianum; i++) if(medt[i]==medt[mediacur])
                        {
                            if(cylz[i]==cylz[mediacur] && cylr[i]==cylr[mediacur] && cyll[i]==cyll[mediacur])
                            {
                                break;
                            }
                        }
                        if(i==medianum)
                        {
                            medianum++;
                            if(sflag)
                            {
                                mediaseg+=4;
                            }
                            else
                            {
                                mediaseg+=2;
                            }

                        }
                    }
                    else if((toLowerCase(taux)).compare("plane")==0)
                    {
                        mediacur=medianum;
                        medt[mediacur]=4;
                        boxx[mediacur]=atof(nextToken(token).c_str());
                        boxy[mediacur]=atof(nextToken(token).c_str());
                        boxz[mediacur]=atof(nextToken(token).c_str());
                        boxl[mediacur]=atof(nextToken(token).c_str());
                        boxw[mediacur]=atof(nextToken(token).c_str());
                        boxh[mediacur]=atof(nextToken(token).c_str());
                        for(i=0; i<medianum; i++)
                        {
                            if(medt[i]==4)
                            {
                                if(boxx[i]==boxx[mediacur] && boxy[i]==boxy[mediacur] && boxz[i]==boxz[mediacur] &&
                                boxl[i]==boxl[mediacur] && boxw[i]==boxw[mediacur] && boxh[i]==boxh[mediacur])
                                {
                                    break;
                                }
                            }
                        }
                        if(i==medianum)
                        {
                            medianum++;
                            mediaseg+=2;
                        }
                    }
                    else
                    {
                        continue;
                    }
;


                    if(atoi(nextToken(token).c_str())>0)
                    {
                        medn[mediacur]=true;
                    }
                    else
                    {
                        medn[mediacur]=false;
                    }

                    vcut[0][mediacur]=atof(nextToken(token).c_str());
                    ecut[0][mediacur]=atof(nextToken(token).c_str());

                    if(atoi(nextToken(token).c_str())>0)
                    {
                        conti[0][mediacur]=true;
                    }
                    else
                    {
                        conti[0][mediacur]=false;
                    }
                    vcut[1][mediacur]=atof(nextToken(token).c_str());
                    ecut[1][mediacur]=atof(nextToken(token).c_str());

                    if(atoi(nextToken(token).c_str())>0)
                    {
                        conti[1][mediacur]=true;
                    }
                    else
                    {
                        conti[1][mediacur]=false;
                    }

                    rho[mediacur]=atof(nextToken(token).c_str());
                    taux=nextToken(token);
                    while(!token->empty())
                    {
                        taux=taux+" "+nextToken(token);
                    }
                    med[mediacur]=taux;
                }
                mediadef_file.close();

            }
            else
            {
                cerr<<"Error: File "<<mediadef<<" cannot be found"<<endl;
                //cerr<<"Error: File "<<mediadef<<" is corrupt: "<<error.tostring()<<endl;
                medianum=1;
                mediaseg=1;
                medn[0]=true;
            }

            mediaopts<<"HI  total number of defined media is ";
            mediaopts<<medianum<<"\n";

        }

        fnm = (int *)calloc(mediaseg,sizeof(int));
        fnx = (double *)calloc(mediaseg,sizeof(double));

        if(gdet==2)
        {
            options<<"HI  radius = "<<Output::f(radius)<<" m  medium = \""<<med[0]<<"\"\n";

        }
        else
        {
            options<<"HI  length = "<<Output::f(length);
            if(gdet==1)
            {
                options<<" m  width = "<<Output::f(width)<<" m  height = "<<Output::f(height)<<" m  medium = \""<<med[0]<<"\"\n";
            }
            else
            {
                options<<" m  radius = "<<Output::f(radius)<<" m  medium = \""<<med[0]+"\"\n";
            }
        }
        //options=(gdet==2?"HI  radius = "+Output::f(radius):"HI  length = "+Output::f(length)+(gdet==1?" m  width = "+Output::f(width)+" m  height = "+Output::f(height):" m  radius = "+Output::f(radius)))+" m  medium = \""+med[0]+"\"\n";

        if(SURF){
            options<<"HI  propagating to surface z = "<<Output::f(surf)<<" m";
            if(FACE)
            {
                if(surf<0)
                {
                    options<<" (up only)";
                }
                else if(surf>0)
                {
                    options<<" (down only)";
                }
            }
            options<<"\n";
        }
        for(i=0; i<medianum; i++){
            options<<"HI  vcut = "<<Output::f(vcut[0][i])<<"/"<<Output::f(vcut[1][i])<<" ecut = "<<Output::f(ecut[0][i])<<"/"<<Output::f(ecut[1][i])<<" MeV";
            if(i==0){
                if(SDEC) options<<"  sdec";
                if(RECC) options<<"  recc";
                if(timef) options<<"  time";
            }
            if(conti[0][i]) options<<"  cont=1";
            if(conti[1][i]) options<<"  cont=2";

            if(i==0){
                options<<"  bs="<<bspar;
                options<<"  ph="<<pncrs;
                options<<"  bb="<<pncbb;
                options<<"  sh="<<pncsh;
                if(lpmef) options<<"  lpm";
                if(scatt) options<<"  scat";
                options<<"  romb="<<romb;
            }
            else options<<"  medium = \""<<med[i]<<"\"";
            if(rho[i]!=1) options<<"  rho="<<Output::f(rho[i]);
            options<<"\n";
        }
        options<<mediaopts.str();
        if(rfix) options<<"HI  Using smart density factor handling, some loss of precision may occur\n";

        Output::AMASIM=amasim;
        PhysicsModel::set_elow(elow);
        PhysicsModel::set_ebig(ebig);
        Propagate::g=romb;
        if(SEED){
            options<<"HI  random number generator seed is set at ";
            options<<seed<<"\n";
            MathModel::set_seed(seed);
        }
        Output::raw=raw;
        cerr<<options.str();

        p1.resize(medianum);
        p2.resize(medianum);
        p3.resize(medianum);
        stringstream setting_stream;
        vector<string> settings;
        settings.resize(medianum);

        for(i=0; i<medianum; i++){

            setting_stream<<Output::f(vcut[0][i])<<" "<<Output::f(ecut[0][i])<<" "<<conti[0][i]<<" ";
            setting_stream<<Output::f(vcut[1][i])<<" "<<Output::f(ecut[1][i])<<" "<<conti[1][i];

            if(frho)
            {
                setting_stream<<" "<<toLowerCase(med[i]);
            }
            else
            {
                setting_stream<<" "<<Output::f(rho[i])<<" "<<toLowerCase(med[i]);

            }
            settings.at(i)=setting_stream.str();

            for(j=0; j<i; j++) if(settings.at(j).compare(settings.at(i))==0) break;
            if(i==j){
//                double vaux=this->vaux, eaux=this->eaux;
//                if(vcut[i]<0)
//                {
//                    eaux=-vcut[i];
//                    vcut[i]=-1;
//                }
//                if(ecut[i]<0)
//                {
//                    vaux=-ecut[i];
//                    ecut[i]=-1;
//                }

                p1.at(i) = new Propagate(med[i], ecut[0][i], vcut[0][i], muta, frho?1:rho[i]);
                p2.at(i) = new Propagate(med[i], ecut[1][i], vcut[1][i], muta, frho?1:rho[i]);
                p3.at(i) = new Propagate(med[i], -1.0, -1.0, muta, frho?1:rho[i]);

                // To enable stopped muon decay
                if(Output::RecDec)
                {
                    p1.at(i)->sdec=SDEC;
                }
                // To enable stopped muon decay
                if(Output::RecDec)
                {
                    p3.at(i)->sdec=SDEC;
                }
                p1.at(i)->contiCorr=conti[0][i]; // To randomize the continuous energy losses, use only for small vcut
                p1.at(i)->exactTime=timef;    // To compute local time of the particle exactly
                p1.at(i)->get_cros()->set_lpm(lpmef);        // Enable lpm and dielectric suppression effects
                p1.at(i)->get_cros()->get_bremsstrahlung()->set_form(bspar);     // Choose parametrization of the bremsstrahlung cross section
                p1.at(i)->get_cros()->get_photonuclear()->set_form(pncrs);     // Choose parametrization of the photon-nucleon cross section
                p1.at(i)->get_cros()->get_photonuclear()->set_bb(pncbb);       // Choose parametrization of the photon-nucleon cross section
                p1.at(i)->get_cros()->get_photonuclear()->set_shadow(pncsh);   // Choose parametrization of the photon-nucleon cross section
                p1.at(i)->molieScat=scatt;    // To enable Moliere scattering of the muon
                p1.at(i)->get_cros()->set_ci(crsci);         // Cross section multiplicative modifier
                p1.at(i)->get_cros()->set_cb(crscb);         // Cross section multiplicative modifier
                p1.at(i)->get_cros()->set_cp(crscp);         // Cross section multiplicative modifier
                p1.at(i)->get_cros()->set_ce(crsce);         // Cross section multiplicative modifier
                p1.at(i)->get_cros()->set_cd(crscd);         // Cross section multiplicative modifier
                p2.at(i)->contiCorr=conti[1][i];  // To randomize the continuous energy losses, use only for small vcut
                p2.at(i)->sdec=SDEC;          // To enable stopped muon decay
                p2.at(i)->recc=RECC;          // To enable printout of continuous energy losses
                p2.at(i)->exactTime=timef;    // To compute local time of the particle exactly
                p2.at(i)->get_cros()->set_lpm(lpmef);        // Enable lpm and dielectric suppression effects
                p2.at(i)->get_cros()->get_bremsstrahlung()->set_form(bspar);     // Choose parametrization of the bremsstrahlung cross section
                p2.at(i)->get_cros()->get_photonuclear()->set_form(pncrs);     // Choose parametrization of the photon-nucleon cross section
                p2.at(i)->get_cros()->get_photonuclear()->set_bb(pncbb);       // Choose parametrization of the photon-nucleon cross section
                p2.at(i)->get_cros()->get_photonuclear()->set_shadow(pncsh);   // Choose parametrization of the photon-nucleon cross section
                p2.at(i)->molieScat=scatt;    // To enable Moliere scattering of the muon
                p2.at(i)->get_cros()->set_ci(crsci);         // Cross section multiplicative modifier
                p2.at(i)->get_cros()->set_cb(crscb);         // Cross section multiplicative modifier
                p2.at(i)->get_cros()->set_cp(crscp);         // Cross section multiplicative modifier
                p2.at(i)->get_cros()->set_ce(crsce);         // Cross section multiplicative modifier
                p2.at(i)->get_cros()->set_cd(crscd);         // Cross section multiplicative modifier
                p3.at(i)->get_cros()->set_lpm(lpmef);        // Enable lpm and dielectric suppression effects
                p3.at(i)->get_cros()->get_bremsstrahlung()->set_form(bspar);     // Choose parametrization of the bremsstrahlung cross section
                p3.at(i)->get_cros()->get_photonuclear()->set_form(pncrs);     // Choose parametrization of the photon-nucleon cross section
                p3.at(i)->get_cros()->get_photonuclear()->set_bb(pncbb);       // Choose parametrization of the photon-nucleon cross section
                p3.at(i)->get_cros()->get_photonuclear()->set_shadow(pncsh);   // Choose parametrization of the photon-nucleon cross section
                p3.at(i)->get_cros()->set_ci(crsci);         // Cross section multiplicative modifier
                p3.at(i)->get_cros()->set_cb(crscb);         // Cross section multiplicative modifier
                p3.at(i)->get_cros()->set_cp(crscp);         // Cross section multiplicative modifier
                p3.at(i)->get_cros()->set_ce(crsce);         // Cross section multiplicative modifier
                p3.at(i)->get_cros()->set_cd(crscd);         // Cross section multiplicative modifier

                ////////// Setting interpolation options

                if(interpolString.compare("")==false)interpolString="all";

                if(interpolString.find("off")==interpolString.npos){

                    p1.at(i)->interpolate(interpolString, tdir+dtct);
                    if(medn[i])
                    {
                        p2.at(i)->interpolate(interpolString, tdir+dtct);
                    }
                    p3.at(i)->interpolate(interpolString, tdir+dtct);

                }
            }
            else{
                p1.at(i)=p1.at(j);
                p2.at(i)=p2.at(j);
                p3.at(i)=p3.at(j);
                if(medn[i])
                {
                    if(!p2.at(i)->jt)
                    {
                        ////////// Setting interpolation options
                        if(interpolString.find("off")==interpolString.npos)p2.at(i)->interpolate(interpolString, tdir+dtct);
                    }
                }
            }

            if(!frho)
            {
                rho[i]=1;
            }
        }

        return true;
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * History lines for the f2k file streams.
     */

     void Amanda::history(){
//        Output.out.println("HI  "+Output.version+"\nHI  "+param);
//        Output.out.print(options);
//        Output.out.println("HI  All "+muta+"'s must have energies from "+
//                           Output.f(p1[0].p.low*1.e-3)+" GeV to "+Output.f(PhysicsModel.ebig*1.e-3)+" GeV");
//        if(USFI) options=" for fixed z = "+Output.f(zset)+" m"; else options="";
//        if(USER) Output.out.println("HI  User line "+usna+" will be enabled"+options);
//        else{ Output.out.println("HI  Use \"-user\" option to enable "+usna+" user line"); rfix=false; }
//        if(rfix) Output.out.println("HI  enforcing compliance with rdmc: exactly one "+usna+" user line per event");
//        if(amasim) Output.out.println("HI  some particles will be substituted with AMASIM-compatible types");
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Define user line block for the f2k file streams.
     */

     void Amanda::defUSER(){
//        if(USER){
//            Output.out.println("USER_DEF "+usna+" NR E_INI E_CPD E_IN E_OUT CDP_X CDP_Y CDP_Z Z_IN Z_OUT");
//            if(rw!=0) Output.out.println("USER_DEF ev_wght WEIGHT DIST_R");
//        }
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize user line block for the f2k file streams.
     */

     void Amanda::iniUSER(){
//        if(USER){ userbf.clear(); imax=0; emax=0; }
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Record user line block for the f2k file streams.
     */

     void Amanda::recUSER(){
//        userbf.addElement(igen+" "+Output.f(e)+" "+(ec>0?Output.f(ec):"?")+" "+
//                          (e1!=0?Output.f(e1):"?")+" "+(e2!=0?Output.f(e2):"?")+" "+
//                          Output.f(nx)+" "+Output.f(ny)+" "+Output.f(nz)+" "+
//                          (h1>0?Output.f(z1):"?")+" "+(h2>0?Output.f(z2):"?"));
//        if(rfix) if(e>emax){ imax=userbf.size()-1; emax=e; }
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Output user line block for the f2k file streams.
     */

     void Amanda::outUSER(){
//        if(USER){
//            if(rfix){
//                if(userbf.size()==0) Output.out.println("US "+usna+" ? ? ? ? ? ? ? ? ? ?");
//                else Output.out.println("US "+usna+" "+userbf.elementAt(imax));
//            }
//            else for(int i=0; i<userbf.size(); i++) Output.out.println("US "+usna+" "+userbf.elementAt(i));
//            if(rw!=0) Output.out.println("US ev_whgt "+fw+" "+hw);
//        }
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize particle buffers for the f2k file streams.
     */

     void Amanda::initBufs(){
//        hist = new stringBuffer(Output::HISTSIZE);
//        userbf = new Vector(Output::HISTSIZE);
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Main parser for the f2k file streams.
     */

     void Amanda::mmcf2k(vector<string> args){
//        Output::I3flag=false;
//        dtct=".amanda";
//        if(!setp(args)) return;

//        //cerr<<" ---  *** Enter your input in F2000 format now ***  --- ");
//        cerr<<" ---  *** Enter your input in F2000 format now ***  --- "<<endl;

//        try{
//            LineNumberReader file = new LineNumberReader(new InputStreamReader(Output.in));
//            stringTokenizer st;
//            string line, taux;
//            Vector buffer = new Vector(Output.HISTSIZE);
//            initBufs();
//            bool tbegin=false;

//            int i, ipar;

//            if((line=file.readLine())!=null){
//                Output.out.println(line);
//                history();
//            }

//            while((line=file.readLine())!=null){
//                st = new stringTokenizer(line);
//                if(!st.hasMoreTokens()) continue;
//                taux=st.nextToken();
//                if("TBEGIN".equals(taux)){
//                    defUSER();
//                    Output.out.println(line);
//                    tbegin=true;
//                    continue;
//                }
//                if(!tbegin){
//                    Output.out.println(line);
//                    continue;
//                }
//                buffer.addElement(line);
//                if("TR".equals(taux)){
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
//                else if("EE".equals(taux) || "END".equals(taux)){
//                    iniUSER();
//                    if(rw!=0) dw=true;

//                    for(i=0; i<buffer.size(); i++){
//                        line=(string)buffer.elementAt(i);
//                        st = new stringTokenizer(line);
//                        taux=st.nextToken();
//                        if("TR".equals(taux)){

//                            igen=Integer.parseInt(st.nextToken());
//                            prnt=st.nextToken();
//                            type=st.nextToken();
//                            if(muta.equals(type) || (muta+"-").equals(type) || (muta+"+").equals(type)){

//                                double x, y, z, th, phi, l, e, t;

//                                x=Double.parseDouble(st.nextToken());
//                                y=Double.parseDouble(st.nextToken());
//                                z=Double.parseDouble(st.nextToken());
//                                th=Double.parseDouble(st.nextToken());
//                                phi=Double.parseDouble(st.nextToken());
//                                if(lfix){
//                                    st.nextToken();
//                                    l=0;
//                                }
//                                else try{
//                                    l=Double.parseDouble(st.nextToken());
//                                    if(l<0) l=0;
//                                    else l*=1.e2;
//                                }catch(NumberFormatException error){
//                                    l=0;
//                                }
//                                e=Double.parseDouble(st.nextToken());
//                                t=Double.parseDouble(st.nextToken());

//                                prop(x, y, z, th, phi, l, e, t);

//                            }
//                            else{
//                                type="TR "+igen+" "+prnt+" "+type+st.nextToken("\n");
//                                Output.out.println(type);
//                            }
//                        }
//                        else if("EM".equals(taux)){
//                            events++;
//                            Output.out.println(line);
//                        }
//                        else if("EE".equals(taux)){
//                            outUSER();
//                            Output.out.println(line);
//                        }
//                        else Output.out.println(line);
//                    }

//                    gens=0;
//                    buffer.clear();
//                    if(rw!=0) dw=false;
//                }

//            }

//            if(events!=0) vertices/=events;
//            cerr<<"events read "+events+" tracks looped "+tracks+
//                               " average vertices "+vertices+" missed volume "+missed);
//        }catch(Exception error){
//            cerr<<"Program finished with exception: "+error.tostring());
//            throw new mmcException("input error");
//        }
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Main propagator routine.
     */

     double Amanda::prop(int igen, int gens, double x, double y, double z, double th, double phi, double l, double e, double t){
        this->igen=igen;
        this->gens=gens;
        return prop(x, y, z, th, phi, l, e, t);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Main propagator routine.
     */

     double Amanda::prop(double x, double y, double z, double th, double phi, double l, double e, double t)
     {
        this->e=e;
        tracks++;
        int j, m, mi, mf, md, mj;

        Propagate *pp=NULL;
        double r1=0, r2=0, r3, rr, h3, n2=0, aux;
        double dx, dy, dz;
        double result;

        double sinth, sinph, costh, cosph;

        sinph=sin(phi*PI/180);
        cosph=cos(phi*PI/180);
        sinth=sin(th*PI/180);
        costh=cos(th*PI/180);

        dx=-sinth*cosph;
        dy=-sinth*sinph;
        dz=-costh;


        switch(gdet)
       {

            case 1:
                setbox(x, y, z, dx, dy, dz, length, width, height);
                break;
            case 2:
                setsph(x, y, z-radius, cosph, sinph, costh, sinth, radius);
                break;
            case 0:
            default:

                setcyl(x, y, z, cosph, sinph, costh, sinth, length, radius);

        }
//cout<<"hx1="<<hx1<<"\t hx2="<<hx2<<endl;
        h1=hx1;
        h2=hx2;

        if(h1<0)
        {
            h1=0;
        }
        if(h2<0)
        {
            h2=0;
        }
        if(h1==h2)
        {
            missed++;
        }
        h3=BIG_;
        if(dw)
        {
            hw=h1+(h2-h1)*MathModel::RandomDouble();
            fw=1;
        }


        if(SURF){
            if(costh!=0) aux=(z-surf)/costh;
            else aux=BIG_;
            if(aux<0) aux=BIG_;
            if(FACE) if((z-surf)*surf<0) aux=0;
            if(h1>aux) h1=aux;
            if(h2>aux) h2=aux;
            if(h3>aux) h3=aux;
        }
//cout<<"h1="<<h1<<"\t h2="<<h2<<"\t h3="<<h3<<endl;
        if(USER){
            if(USFI){
                if(costh!=0) n2=(z-zset)/costh;
                else n2=0;
            }
            else n2=x*sinth*cosph+y*sinth*sinph+z*costh;
            if(n2<h1) flag=1;
            else if(n2<h2) flag=2;
            else if(n2<h3) flag=3;
            else flag=0;
        }

        fnu=1;
        fnm[0]=0;
        fnx[0]=BIG_;
        for(m=1; m<medianum; m++)
        {
//            cout<<"p1.at("<<m-1<<")->jt = "<<p1.at(m-1)->jt<<endl;
//            cout<<"p2.at("<<m-1<<")->jt = "<<p2.at(m-1)->jt<<endl;
//            cout<<"p3.at("<<m-1<<")->jt = "<<p3.at(m-1)->jt<<endl;
//cout<<"m="<<m<<"\t medianum="<<medianum<<"\t medt[m]="<<medt[m]<<"\t abs(medt[m])="<<abs(medt[m])<<endl;
            switch(abs(medt[m]))
            {
                case 1:
                    setsph(x, y, z-sphz[m], cosph, sinph, costh, sinth, sphr[m]);
                    break;
                case 2:
                    setbox(x-boxx[m], y-boxy[m], z-boxz[m], dx, dy, dz, boxl[m], boxw[m], boxh[m]);
                    break;
                case 3:
                    setcyl(x, y, z-cylz[m], cosph, sinph, costh, sinth, cyll[m], cylr[m]);
                    break;
                case 4:
                    setpln(x-boxx[m], y-boxy[m], z-boxz[m], dx, dy, dz, boxl[m], boxw[m], boxh[m]);
                    break;
                default: hx1=0; hx2=0;
            }
//cout<<"Line 1659: hx1="<<hx1<<"\t hx2="<<hx2<<endl;
            for(int ij=0; ij<(medt[m]>0?1:2); ij++)
            {
                if(medt[m]>0)
                {
                    r1=hx1;
                    r2=hx2;
                }
                else
                {
                    if(ij==0)
                    {
                        r1=0;
                        r2=hx1;
                    }
                    else
                    {
                        r1=hx2;
                        r2=BIG_;
                    }
                }
                if(r2>r1)
                {
                    md=0;
                    mi=0;
                    if(r1>0)
                    {
                        md++;
                        for(mi=0; mi<fnu; mi++)
                        {
                            if(fnx[mi]>r1)
                            {
                                break;
                            }
                        }
                    }
                    mf=0;
                    if(r2>0)
                    {
                        md++;
                        for(mf=mi; mf<fnu; mf++)
                        {
                            if(fnx[mf]>r2)
                            {
                                break;
                            }
                        }
                    md+=mi-mf;
                    }
                    if(md>0)
                    {
                        for(mj=fnu-1; mj>=mf; mj--)
                        {
                            fnx[mj+md]=fnx[mj];
                            fnm[mj+md]=fnm[mj];
                        }
                    }
                    else if(md<0)
                    {
                        for(mj=mf; mj<fnu; mj++)
                        {
                            fnx[mj+md]=fnx[mj];
                            fnm[mj+md]=fnm[mj];
                        }
                    }
                    fnu+=md;
                    if(r1>0)
                    {
                        fnx[mi]=r1;
                        mi++;
                    }
                    if(r2>0)
                    {
                        fnx[mi]=r2;
                        fnm[mi]=m;
                    }
                }
            }
        }
//cout<<"r1="<<r1<<"\t r2="<<r2<<endl;
        if(DEBUG){
            for(j=0; j<fnu; j++)
            {
                //Output.out.println(j+"\t"+Output.f(fnx[j])+"\t"+fnm[j]+"\t"+p1[fnm[j]].m.name+"\t"+Output.f(rho[fnm[j]]));
            }
           // Output.out.println();
        }

        if(Output::I3flag)
        {
            I3hist.clear();
        }
        else
        {
            hist.str("");
            hist.clear();
        }
        if(USER)
        {
            e1=0;
            e2=0;
            ec=0;
            elost=0;
        }
        pp=NULL;
        result=e;
        rr=0;

        for(j=0; j<fnu; j++)
        {
            m=fnm[j];
            aux=fnx[j];
            r1=min(h1, aux);
            r2=min(h2, aux);
            r3=min(h3, aux);
            p1.at(m)->set_rho(rho[m]);
            p2.at(m)->set_rho(rho[m]);
            p3.at(m)->set_rho(rho[m]);
//cout<<"r1="<<r1<<"\t r2="<<r2<<"\t r3="<<r3<<endl;
            if(pp==NULL)
            {
//cout<<"r1="<<r1<<endl;

                if(phi<180)
                {
                    result=propagateTo(r1, e, n2, 1, p1.at(m), igen, gens,t*1.e-9, x*1.e2, y*1.e2, z*1.e2, 180-th, phi+180);
                }
                else
                {
                    result=propagateTo(r1, e, n2, 1, p1.at(m), igen, gens,t*1.e-9, x*1.e2, y*1.e2, z*1.e2, 180-th, phi-180);
                }


            }
            else
            {
//cout<<"r1-rr="<<r1-rr<<endl;

                result=propagateTo(r1-rr, result, n2-rr, 1, p1.at(m), 0, 0, pp->get_particle()->t, pp->get_particle()->x, pp->get_particle()->y, pp->get_particle()->z, pp->get_particle()->theta, pp->get_particle()->phi);
            }            
            pp=p1[m];
            l+=pp->get_particle()->r;
//cout<<"l="<<l<<"\t result="<<result<<endl;
            if(USER)
            {
                if(r1>=rr && r1==h1)
                {
                    e1=result;
                    if(Output::I3flag)
                    {
                        pI->ti=pp->get_particle()->t;
                    }
                }
            }

            if(Output::RecDec)
            {
                if(Output::I3flag)
                {

                    I3hist.insert(I3hist.begin(), pp->get_output()->I3hist.begin(), pp->get_output()->I3hist.end());
                }
                else
                {
                    hist<<pp->get_output()->history.str();
                    gens=pp->get_output()->gens;
                }
            }
//cout<<"e2="<<e2<<"\t r2="<<r2<<"\t rr="<<rr<<"\t h2="<<h2 <<endl;
//cout<<"r1="<<r1<<"\t r2="<<r2<<"\t rr="<<rr<<"\t result="<<result<<endl;
            if(result>0)
            {

                if(r1>rr)
                {
                    rr=r1;
                    if(USER)
                    {
                        elost+=result;
                    }
                }
                if(dw)
                {
                    p2.at(m)->dw=true;
                    p2.at(m)->rw=rw;
                    p2.at(m)->hw=(hw-rr)*1.e2;
                }
//cout<<"m=========================================================== " <<m<<endl;
//cout<<"p2.at(m)->jt = "<<p2.at(m)->jt<<endl;
                result=propagateTo(r2-rr, result, n2-rr, 2, p2.at(m), igen, gens, pp->get_particle()->t, pp->get_particle()->x, pp->get_particle()->y, pp->get_particle()->z, pp->get_particle()->theta, pp->get_particle()->phi);
                pp=p2.at(m);
                l+=pp->get_particle()->r;
                if(USER)
                {
                    if(r2>=rr && r2==h2)
                    {
                        e2=result;
                        if(Output::I3flag)
                        {
                            pI->tf=pp->get_particle()->t;
                        }
                    }
                }
                if(dw)
                {
                    if(!p2.at(m)->dw)
                    {
                        dw=false;
                        fw=p2.at(m)->rw;
                        hw=p2.at(m)->hw*1.e-2+rr;
                    }
                    else
                    {
                        p2.at(m)->dw=false;
                    }
                    p2.at(m)->rw=0;
                    p2.at(m)->hw=0;
                }
                if(Output::I3flag)
                {
                    I3hist.insert(I3hist.begin(), pp->get_output()->I3hist.begin(), pp->get_output()->I3hist.end());
                }
                else
                {
                    hist<<pp->get_output()->history.str();
                }
                vertices+=pp->get_output()->gens-gens;
                gens=pp->get_output()->gens;
                if(result>0)
                {
                    if(r2>rr)
                    {
                        rr=r2;
                    }
                    if(USER)
                    {
                        elost-=result;


                        if(elost<0)
                        {
                            elost=0;
                        }
                    }
//cout<<"r3-rr="<<r3-rr<<endl;
                    result=propagateTo(r3-rr, result, n2-rr, 3, p3.at(m), igen, gens, pp->get_particle()->t, pp->get_particle()->x, pp->get_particle()->y, pp->get_particle()->z, pp->get_particle()->theta, pp->get_particle()->phi);
                    pp=p3.at(m);
                    l+=pp->get_particle()->r;
                    if(Output::RecDec)
                    {
                        if(Output::I3flag)
                        {
                            I3hist.insert(I3hist.begin(), pp->get_output()->I3hist.begin(), pp->get_output()->I3hist.end());
                        }
                        else
                        {
                            hist<<pp->get_output()->history.str();
                        }
                        gens=pp->get_output()->gens;
                    }
                    if(result>0)
                    {
                        if(r3>rr)
                        {
                            rr=r3;
                        }
                    }
                }
            }
            if(result<=0)
            {
                break;
            }
        }

        if(USER){
            if(h1==0)
            {
                e1=0;
            }
            if(h2==0)
            {
                e2=0;
            }
            if(e1<0)
            {
                e1-=rr;
            }
            if(e2<0)
            {
                e2-=rr;
            }
            z1=z-costh*h1;
            z2=z-costh*h2;
            nx=x-sinth*cosph*n2;
            ny=y-sinth*sinph*n2;
            nz=z-costh*n2;

            if(!Output::I3flag)
            {
                recUSER();
            }
            else{
                pI->xi=x-sinth*cosph*h1;
                pI->yi=y-sinth*sinph*h1;
                pI->zi=z-costh*h1;
                pI->xf=x-sinth*cosph*h2;
                pI->yf=y-sinth*sinph*h2;
                pI->zf=z-costh*h2;
                pI->xc=nx;
                pI->yc=ny;
                pI->zc=nz;
                pI->Ei=e1;
                pI->Ef=e2;
                if(ec>0)
                {
                    pI->Ec=ec;
                }
                else
                {
                    pI->Ec=0;
                }

                pI->Elost=elost;
 //cout<<"pI->Elost="<<pI->Elost<<"\t elost"<<elost<<endl;
//  cout<<"pI->yi="<<pI->yi<<endl;
//  cout<<"pI->zi="<<pI->zi<<endl;
//  cout<<"pI->xf="<<pI->xf<<endl;
//  cout<<"pI->yf="<<pI->yf<<endl;
//  cout<<"pI->zf="<<pI->zf<<endl;
//  cout<<"pI->xc="<<pI->xc<<endl;
//  cout<<"pI->yc="<<pI->yc<<endl;
//  cout<<"pI->zc="<<pI->zc<<endl;
//  cout<<"pI->Ei="<<pI->Ei<<endl;
//  cout<<"pI->Ef="<<pI->Ef<<endl;
//  cout<<"pI->Ec="<<pI->Ec<<endl;
//  cout<<"l="<<l<<"\t SURF="<<SURF<<endl;
            }
        }

        if(SURF)
        {
            x=pp->get_particle()->x*1.e-2;
            y=pp->get_particle()->y*1.e-2;
            z=pp->get_particle()->z*1.e-2;
            th=180-pp->get_particle()->theta;
            phi=pp->get_particle()->phi;
            if(phi<180)
            {
                phi=phi+180;
            }
            else
            {
                phi=phi-180;
            }
            if(result>0)
            {
                e=result;
                l=-BIG_;
            }
            else
            {
                e=0;
                l=0;
            }
            t=pp->get_particle()->t*1.e9;
        }
        if(!Output::I3flag)
        {
//            Output.out.println("TR "+igen+" "+prnt+" "+type+" "+Output.f(x)+" "+Output.f(y)+" "+Output.f(z)+" "+
//                               Output.f(th)+" "+Output.f(phi)+" "+(l>=0?Output.f(l*1.e-2):"?")+" "+
//                               Output.f(e)+" "+Output.f(t));
//            Output.out.print(hist);
        }

        return l;
    }

