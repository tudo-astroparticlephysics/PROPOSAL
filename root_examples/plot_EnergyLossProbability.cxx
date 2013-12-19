#include <iostream>
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/Ionization.h"
#include "PROPOSAL/Epairproduction.h"
#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/Output.h"
#include "time.h"
#include <fstream>
#include "TCanvas.h"
#include "TH2D.h"
#include "TTree.h"

using namespace std;


int main(int argc, char** argv)
{
    cout<<"\n------------------------------------------------------------------\n"
        <<"This is an example to plot the probabilities for are energy losses.\n"
        <<"A ROOT-file named EnergyLosses.root will be generated,\n"
        <<"which contains TCanvas (stochastic) with an overview of all losses.\n"
        <<"You have to call the app. with 4 parameters: \n"
        <<"ecut seed vcut #particles \n"
        <<"e.g. 500 42 -1 1000\n"
        <<"------------------------------------------------------------------\n"
        <<endl;

    if(argc<5)
    {
        cerr<<"Not enough parameters(" << argc-1 << ") given! Expected 4 parameters! \n Aborting application!";
        exit(1);
    }
    stringstream ss;

    double ecut = atof(argv[1]);
    int seed = atoi(argv[2]);
    double vcut = atof(argv[3]);
    int statistic = atoi(argv[4]);
    string vcut_o;
    string ecut_o;

    ss<<"EnergyLosses.root";
    Medium *med = new Medium("ice",1.);
    EnergyCutSettings* cuts = new EnergyCutSettings(ecut,vcut);

    TFile *file     =   new TFile(ss.str().c_str(),"RECREATE");


    Propagator *pr = new Propagator(med,cuts,"mu","resources/tables",false,false,false,false,1);


    TH2D *ioniz = new TH2D("ioniz","Ionization;primary log10(energy/MeV);secondary log10(energy/MeV)",512,2,14,512,-4,14);
    TH2D *brems = new TH2D("brems","Bremsstrahlung;primary log10(energy/MeV);secondary log10(energy/MeV)",512,2,14,512,-4,14);
    TH2D *photo = new TH2D("photo","Photonuclear Effect;primary log10(energy/MeV);secondary log10(energy/MeV)",512,2,14,512,-4,14);
    TH2D *epair = new TH2D("epair","e-Pair Production;primary log10(energy/MeV);secondary log10(energy/MeV)",512,2,14,512,-4,14);

    double sec_energy;
    int ww;
    double energy;

    TTree *tree = new TTree("secondarys","secondarys");
    tree->SetDirectory(0);
    tree->Branch("sec_energy",&sec_energy,"sec_energy/D");
    tree->Branch("energy",&energy,"energy/D");
    tree->Branch("ww",&ww,"ww/I");


    pr->set_seed(seed);


    for(int i =0;i<statistic;i++)
    {
        pr->GetParticle()->SetEnergy(pow(10,14));
        pr->GetParticle()->SetPropagatedDistance(0);
        pr->Propagate(1e20);
        cout<<i<<endl;
        vector<Particle*> bla (Output::getInstance().GetSecondarys());
        for(unsigned int k =0 ; k<bla.size();k++)
        {
            sec_energy  =   bla.at(k)->GetEnergy();
            energy  =   bla.at(k)->GetParentParticleEnergy();

            if(bla.at(k)->GetType()==-1002)
            {
                ww=1;
                ioniz->Fill(log10(energy),log10(sec_energy));

            }
            else if(bla.at(k)->GetType()==-1001)
            {
                ww=2;
                brems->Fill(log10(energy),log10(sec_energy));

            }
            else if(bla.at(k)->GetType()==-1004)
            {
                ww=3;
                photo->Fill(log10(energy),log10(sec_energy));

            }
            else if(bla.at(k)->GetType()==-1003)
            {
                ww=4;
                epair->Fill(log10(energy),log10(sec_energy));

            }

            tree->Fill();

        }
        Output::getInstance().ClearSecondaryVector();



    }

    TCanvas *can = new TCanvas("stochatic","stochastic",1400,1050);
    can->Divide(2,2);
    can->cd(1);
    can->SetGridx();
    can->SetGridy();
    can->SetTickx();
    can->SetTicky();
    gPad->SetLogz(1);
    ioniz->Draw("colz");

    can->cd(2);
    can->SetGridx();
    can->SetGridy();
    can->SetTickx();
    can->SetTicky();
    gPad->SetLogz(1);
    brems->Draw("colz");

    can->cd(3);
    can->SetGridx();
    can->SetGridy();
    can->SetTickx();
    can->SetTicky();
    gPad->SetLogz(1);
    photo->Draw("colz");

    can->cd(4);
    can->SetGridx();
    can->SetGridy();
    can->SetTickx();
    can->SetTicky();
    gPad->SetLogz(1);
    epair->Draw("colz");

    can->Write();
    ioniz->Write();
    epair->Write();
    brems->Write();
    photo->Write();

    //tree->Write();

    file->Close();
}
