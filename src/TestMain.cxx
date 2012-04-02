#include "MMCexception.h"
#include "MathModel.h"
//#include "typeinfo"
#include "Output.h"
#include <cmath>
#include "Mpairproduction.h"
#include "MpairContinuous.h"
#include "PhysicsModel.h"
#include "Amanda.h"
#include "Photonuclear.h"
#include "BremsContinuous.h"
#include "EpairContinuous.h"
#include "PhotoContinuous.h"
#include "IonizContinuous.h"

#include "Photonuclear.h"
#include "Ionizationloss.h"

#include "EpairStochastic.h"
#include "IonizStochastic.h"
#include "BremsStochastic.h"
#include "PhotoStochastic.h"
#include "Mpairproduction.h"
#include "MpairContinuous.h"
#include "PhysicsModel.h"
#include "Amanda.h"
#include <sstream>
#include <fstream>

using namespace std;




std::ifstream randomIn;
//double PhysicsModel::ebig_ =BIGENERGY;

int main(){


//        Propagate *pr1 = new Propagate("ice" ,-1,0.005,"mu");
//        Propagate *pr2 = new Propagate("frejus rock" ,-1,0.005,"mu");
//        pr1->get_cros()->get_bremsstrahlung()->set_form(1);
//        pr2->get_cros()->get_bremsstrahlung()->set_form(1);
//        pr1->get_cros()->get_photonuclear()->set_form(3);
//        pr2->get_cros()->get_photonuclear()->set_form(3);
//        pr1->get_cros()->get_photonuclear()->set_bb(2);
//        pr2->get_cros()->get_photonuclear()->set_bb(2);
//        pr1->get_cros()->get_photonuclear()->set_shadow(2);
//        pr2->get_cros()->get_photonuclear()->set_shadow(2);


//        pr1->set_contiCorr(true);
//        pr2->set_contiCorr(true);
//        pr1->interpolate("all");
//        pr2->interpolate("all");

//        fstream out1;
//        fstream out2;


//        double energy;

//        string filename1="/home/koehne/Desktop/Range/Ice/Range_IcecubeStandard1e";
//        string filename2="/home/koehne/Desktop/Range/Frock/Range_IcecubeStandard1e";
//        string extension=".txt";

//        stringstream ss1;
//        stringstream ss2;

//        for(double i = 2;i<=14;){

//            energy =pow(10,i);
//            ss1<<filename1<<i<<extension;
//            ss2<<filename2<<i<<extension;

//            out1.precision(10);
//            out2.precision(10);

//            out1.open(ss1.str().c_str(),ios::out);
//            out2.open(ss2.str().c_str(),ios::out);

//            cout<<i<<endl;

//            for(int j =0; j<1e3;j++){

//                out1<<pr1->propagateTo(1e20,energy)<<endl;
//                out2<<pr2->propagateTo(1e20,energy)<<endl;
//            }
//            i=i+0.1;

//            ss1.str("");
//            ss2.str("");
//            ss1.clear();
//            ss2.clear();

//            out1.close();
//            out2.close();

//        }

//    Output::I3flag=true;
//    vector<PROPOSALParticle*> I3p;

    Propagate *pr1 = new Propagate("water" ,-1,0.2,"mu");
    Propagate *pr2 = new Propagate("water" ,-1,0.2,"mu");
    Propagate *pr3 = new Propagate("water" ,-1,0.05,"mu");
    Propagate *pr4 = new Propagate("water" ,-1,0.05,"mu");
    Propagate *pr5 = new Propagate("water" ,-1,0.01,"mu");
    Propagate *pr6 = new Propagate("water" ,-1,0.01,"mu");
    Propagate *pr7 = new Propagate("water" ,-1,0.001,"mu");
    Propagate *pr8 = new Propagate("water" ,-1,0.001,"mu");

//    pr1->get_cros()->get_bremsstrahlung()->set_form(1);
//    pr2->get_cros()->get_bremsstrahlung()->set_form(1);
//    pr3->get_cros()->get_bremsstrahlung()->set_form(1);
//    pr4->get_cros()->get_bremsstrahlung()->set_form(1);
//    pr5->get_cros()->get_bremsstrahlung()->set_form(1);
//    pr6->get_cros()->get_bremsstrahlung()->set_form(1);
//    pr7->get_cros()->get_bremsstrahlung()->set_form(1);
//    pr8->get_cros()->get_bremsstrahlung()->set_form(1);

    pr1->get_cros()->get_photonuclear()->set_form(2);
    pr2->get_cros()->get_photonuclear()->set_form(2);
    pr3->get_cros()->get_photonuclear()->set_form(2);
    pr4->get_cros()->get_photonuclear()->set_form(2);
    pr5->get_cros()->get_photonuclear()->set_form(2);
    pr6->get_cros()->get_photonuclear()->set_form(2);
    pr7->get_cros()->get_photonuclear()->set_form(2);
    pr8->get_cros()->get_photonuclear()->set_form(2);

    pr1->get_cros()->get_photonuclear()->set_bb(4);
    pr2->get_cros()->get_photonuclear()->set_bb(4);
    pr3->get_cros()->get_photonuclear()->set_bb(4);
    pr4->get_cros()->get_photonuclear()->set_bb(4);
    pr5->get_cros()->get_photonuclear()->set_bb(4);
    pr6->get_cros()->get_photonuclear()->set_bb(4);
    pr7->get_cros()->get_photonuclear()->set_bb(4);
    pr8->get_cros()->get_photonuclear()->set_bb(4);

//    pr1->get_cros()->get_photonuclear()->set_shadow(2);
//    pr2->get_cros()->get_photonuclear()->set_shadow(2);
//    pr3->get_cros()->get_photonuclear()->set_shadow(2);
//    pr4->get_cros()->get_photonuclear()->set_shadow(2);
//    pr5->get_cros()->get_photonuclear()->set_shadow(2);
//    pr6->get_cros()->get_photonuclear()->set_shadow(2);
//    pr7->get_cros()->get_photonuclear()->set_shadow(2);
//    pr8->get_cros()->get_photonuclear()->set_shadow(2);


    pr2->set_contiCorr(true);
    pr4->set_contiCorr(true);
    pr6->set_contiCorr(true);
    pr8->set_contiCorr(true);

    pr1->interpolate("all");
    pr2->interpolate("all");
    pr3->interpolate("all");
    pr4->interpolate("all");
    pr5->interpolate("all");
    pr6->interpolate("all");
    pr7->interpolate("all");
    pr8->interpolate("all");


    fstream out;
    fstream out1;
    fstream out2;
    fstream out3;
    fstream out4;
    fstream out5;
    fstream out6;
    fstream out7;
    fstream out8;

    out.open("/home/koehne/Desktop/FinalEnergyPaper/Info5.txt", ios::out | ios::app);
    out1.open("/home/koehne/Desktop/FinalEnergyPaper/std1_water_10e6TeV_40km_vcut02_nocont.txt", ios::out);
    out2.open("/home/koehne/Desktop/FinalEnergyPaper/std1_water_10e6TeV_40km_vcut02_cont.txt", ios::out);
    out3.open("/home/koehne/Desktop/FinalEnergyPaper/std1_water_10e6TeV_40km_vcut005_nocont.txt", ios::out);
    out4.open("/home/koehne/Desktop/FinalEnergyPaper/std1_water_10e6TeV_40km_vcut005_cont.txt", ios::out);
    out5.open("/home/koehne/Desktop/FinalEnergyPaper/std1_water_10e6TeV_40km_vcut001_nocont.txt", ios::out);
    out6.open("/home/koehne/Desktop/FinalEnergyPaper/std1_water_10e6TeV_40km_vcut001_cont.txt", ios::out);
    out7.open("/home/koehne/Desktop/FinalEnergyPaper/std1_water_10e6TeV_40km_vcut0001_nocont.txt", ios::out);
    out8.open("/home/koehne/Desktop/FinalEnergyPaper/std1_water_10e6TeV_40km_vcut0001_cont.txt", ios::out);

    out1.precision(10);
    out2.precision(10);
    out3.precision(10);
    out4.precision(10);
    out5.precision(10);
    out6.precision(10);
    out7.precision(10);
    out8.precision(10);

    out<<"File: \t";
    out<<"std1_water_10e6TeV_40km_vcut02_nocont.txt \t";
    out<<"brems=def \t";
    out<<"photo=2 \t";
    out<<"bb=4 \t";
    out<<"shadow=def"<<endl;

    out<<"File: \t";
    out<<"std1_water_10e6TeV_40km_vcut02_cont.txt \t";
    out<<"brems=def \t";
    out<<"photo=2 \t";
    out<<"bb=4 \t";
    out<<"shadow=def"<<endl;

    out<<"File: \t";
    out<<"std1_water_10e6TeV_40km_vcut005_nocont.txt \t";
    out<<"brems=def \t";
    out<<"photo=2 \t";
    out<<"bb=4 \t";
    out<<"shadow=def"<<endl;

    out<<"File: \t";
    out<<"std1_water_10e6TeV_40km_vcut005_cont.txt \t";
    out<<"brems=def \t";
    out<<"photo=2 \t";
    out<<"bb=4 \t";
    out<<"shadow=def"<<endl;

    out<<"File: \t";
    out<<"std1_water_10e6TeV_40km_vcut001_nocont.txt \t";
    out<<"brems=def \t";
    out<<"photo=2 \t";
    out<<"bb=4 \t";
    out<<"shadow=def"<<endl;

    out<<"File: \t";
    out<<"std1_water_10e6TeV_40km_vcut001_cont.txt \t";
    out<<"brems=def \t";
    out<<"photo=2 \t";
    out<<"bb=4 \t";
    out<<"shadow=def"<<endl;

    out<<"File: \t";
    out<<"std1_water_10e6TeV_40km_vcut0001_nocont.txt \t";
    out<<"brems=def \t";
    out<<"photo=2 \t";
    out<<"bb=4 \t";
    out<<"shadow=def"<<endl;

    out<<"File: \t";
    out<<"std1_water_10e6TeV_40km_vcut0001_cont.txt \t";
    out<<"brems=def \t";
    out<<"photo=2 \t";
    out<<"bb=4 \t";
    out<<"shadow=def"<<endl;


    for(int i =0;i<1e6;i++){
        if(i%10000==0)cout<<i<<endl;
        pr1->get_particle()->setEnergy(1e12);
        pr2->get_particle()->setEnergy(1e12);
        pr3->get_particle()->setEnergy(1e12);
        pr4->get_particle()->setEnergy(1e12);
        pr5->get_particle()->setEnergy(1e12);
        pr6->get_particle()->setEnergy(1e12);
        pr7->get_particle()->setEnergy(1e12);
        pr8->get_particle()->setEnergy(1e12);

        out1<<pr1->propagateTo(4e6,1e12)<<endl;
        out2<<pr2->propagateTo(4e6,1e12)<<endl;
        out3<<pr3->propagateTo(4e6,1e12)<<endl;
        out4<<pr4->propagateTo(4e6,1e12)<<endl;
        out5<<pr5->propagateTo(4e6,1e12)<<endl;
        out6<<pr6->propagateTo(4e6,1e12)<<endl;
        out7<<pr7->propagateTo(4e6,1e12)<<endl;
        out8<<pr8->propagateTo(4e6,1e12)<<endl;


    }
//    out.close();
//    cout.precision(16);
//    randomIn.open("randomNumbers.txt");

string filename = "/home/koehne/IcecubeMMC++/Daten/spectrum_test2.dat";
string opts;

    opts="-mediadef=mediadef -cont -lpm -raw -romb=5 -scat -user -sdec -time -lpm -bs=1 -ph=3 -bb=2 -sh=2";
//double energy;
double theta;
double phi;
double x;
double y;
double z;
double r;
double t;
Amanda *a = new Amanda();
a->setup(opts);
PROPOSALParticle *p;
bool raw=true;

double exp= 2.0;
double par=1e6;

//
//for(int j =0;j<10000;j++){
//    //cout<<j<<endl;

//    energy=par*pow((1-MathModel::RandomDouble()),(1/(-exp+1)));
//    if(j%1000==0)cout<<j/1000+1<<endl;
//    phi=0;
//    theta=0;
//    x=0;
//    y=0;
//    z=-50000;
//    r=0;
//    t=1e-6;
//    p = new PROPOSALParticle(0,0,"mu",x,y,z,theta,phi,1e7,t,r);

//    p->setEnergy(energy);

//    vector<PROPOSALParticle*> I3p=a->propagate(p);
//    if(raw){
//        Output::particleToFileRaw(p,filename,true);
//    }
//    else{
//        Output::particleToFileAscii(p,filename,true);
//    }

//        for(int i =0; i<I3p.size();i++){
//            if(raw){
//                Output::particleToFileRaw(I3p.at(i),filename,false);
//            }
//            else
//            {
//                Output::particleToFileAscii(I3p.at(i),filename,false);

//            }
//        }

//    delete p;

//        for(int i =0; i<I3p.size();i++){
//            delete I3p.at(i);
//        }
//}


};










