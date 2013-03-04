#include "config.h"
#include <iostream>
#include "PROPOSAL/Propagate.h"
#include <sstream>
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/BremsContinuous.h"
#include "PROPOSAL/Epairproduction.h"
#include "PROPOSAL/EpairContinuous.h"
#include "PROPOSAL/Photonuclear.h"
#include "PROPOSAL/PhotoContinuous.h"
#include "PROPOSAL/Ionizationloss.h"
#include "PROPOSAL/IonizContinuous.h"


#ifndef ROOT_SUPPORT
    #define ROOT_SUPPORT 0
#endif
//#if ROOT_SUPPORT
//    #include "TH1F.h"
//    #include "TTree.h"
//#endif

using namespace std;

int main(){


//    double min=4.;
//    double max=3.;
//    double aux;

//    for(long int i=0;i<1e10;i++){

////        aux =   min;
////        min =   max;
////        max =   aux;
//        SWAP(min,max,double);

//    }

//    Propagate *pr1 = new Propagate("ice",500.,0.05,"mu");
//    pr1->interpolate("all");
//    double energy;
//    ofstream out;
//    out.open("/home/koehne/Daten/referenz2.txt");
//    out.precision(14);

//    for(int i =0;i<1e4;i++){
//        if(i%100==0)cout<<i/100+1<<endl;
//        energy=i*1e4;
//        out<<pr1->propagateTo(1e8,energy)<<endl;

//    }

//    out.close();
//    return 0;


//    Propagate *pr1 = new Propagate("ice",500.,0.05,"mu");
//    Propagate *pr2 = new Propagate("ice",500.,0.05,"mu");
//    pr1->setCShigh();
//    pr2->setCSlow();
//    pr1->set_contiCorr(true);
//    pr2->set_contiCorr(true);
//    pr1->interpolate("all","");
//    pr2->interpolate("all","");

//    double energy;
//    ofstream out1;
//    out1.open("/home/koehne/Daten/high_r1e6_e1e8.txt");
//    out1.precision(14);

//    ofstream out2;
//    out2.open("/home/koehne/Daten/low_r1e6_e1e8.txt");
//    out2.precision(14);

//    for(int i =0;i<1e5;i++){
//        if(i%1000==0)cout<<i/1000+1<<endl;
//        energy=1e8;
//        out1<<pr1->propagateTo(1e6,energy)<<endl;
//        out2<<pr2->propagateTo(1e6,energy)<<endl;


//    }

//    out1.close();
//    out2.close();

    stringstream ss;
    ofstream out;
    string muta="mu";

    double energy=1;

    Propagate *pr1 = new Propagate("ice",-1,-1,"mu");
    //Output::raw=true;
    for(int bspar=1;bspar<5;bspar++)
    {
        for(int pncrs=1;pncrs<5;pncrs++)
        {
            for(int pncbb=1;pncbb<5;pncbb++)
            {
                for(int pncsh=1;pncsh<3;pncsh++)
                {

//                    if(bspar<1 || bspar>4)
//                    {
//                        cout<<"Warning: bs is not a valid number"<<endl;
//                        bspar   =   1;
//                    }

//                    if(pncrs<1 || pncrs>4 || (pncrs==2 && !(muta.compare("mu")==0 || muta.compare("tau")==0)))
//                    {
//                        cout<<"Warning: ph is not a valid number"<<endl;
//                        pncrs   =   1;
//                    }

//                    if(((pncrs==1 || pncrs==2) && (pncbb<1 || pncbb>4)) ||
//                            (pncrs==3 && (pncbb<1 || pncbb>2)) || (pncrs==4 && pncbb!=1))
//                    {
//                        cout<<"Warning: bb is not a valid number"<<endl;
//                        pncbb   =   1;
//                    }

//                    if(((pncrs==1 || pncrs==2) && (pncsh!=1)) || ((pncrs>2) && (pncsh<1 || pncsh>2)))
//                    {
//                        cout<<"Warning: sh is not a valid number"<<endl;
//                        pncsh   =   1;
//                    }

                    pr1->setCSform(bspar,pncrs,pncbb,pncsh,"mu");
                    //pr1->set_contiCorr(true);
                    pr1->interpolate("all","");
                    ss<<"/home/koehne/Daten/dEdx_"<<bspar<<pncrs<<pncbb<<pncsh<<".txt";
                    out.open(ss.str().c_str());
                    cout<<"-----"<<ss.str().c_str()<<"-----"<<endl;
                    ss.clear();
                    ss.str("");

                    for(double i=1;i<12;i=i+0.1)
                    {
                        energy=100*pow(10,i);

                        pr1->get_particle()->setEnergy(energy);
                        pr1->get_cros()->get_ionization()->get_Continuous()->get_particle()->setEnergy(energy);
                        pr1->get_cros()->get_bremsstrahlung()->get_Continuous()->get_particle()->setEnergy(energy);
                        pr1->get_cros()->get_epairproduction()->get_Continuous()->get_particle()->setEnergy(energy);
                        pr1->get_cros()->get_photonuclear()->get_Continuous()->get_particle()->setEnergy(energy);
                        out<<energy<<"\t";
                        out<<pr1->get_cros()->get_ionization()->get_Continuous()->dEdx()<<"\t";
                        out<<pr1->get_cros()->get_bremsstrahlung()->get_Continuous()->dEdx()<<"\t";
                        out<<pr1->get_cros()->get_epairproduction()->get_Continuous()->dEdx()<<"\t";
                        out<<pr1->get_cros()->get_photonuclear()->get_Continuous()->dEdx()<<endl;

                    }

                    //pr1 = new Propagate("ice",500.,0.05,"mu");
                    out.close();

                }


            }


        }


    }

    return 0;

}


