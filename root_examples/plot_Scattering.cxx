#include "TH1D.h"
#include "TCanvas.h"
#include "PROPOSAL/Propagator.h"
#include "TFile.h"
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include "TGraphErrors.h"


using namespace std;

int main()
{

    cout<<"This program will propagte 1e6 muons with an energy of 9 TeV."<<endl;
    Propagator *prpgtr = new Propagator("resources/configuration_IceOnly");
    prpgtr->set_seed(1234);

    stringstream hist_name;
    stringstream hist_name2;
    stringstream hist_title;
    stringstream hist_title2;


    vector<TH1D*> tmp_hist_dev_vec;
    vector<TH1D*> tmp_hist_angle_vec;

    vector<double> e;

    int number_of_particles =   3e5;

    double distance;
    double deviation;
    double angle;

    double Emin=2e2;
    double Emax=1e12;
    double multiplier = 1;


    for(double energy = 0; energy <= Emax ; energy *= multiplier)
    {        
        if(energy<1e3)
        {
            multiplier=1;
            energy = energy + Emin;
        }

        else if(energy<1e4)
        {
            multiplier=1;
            energy = energy + 5*Emin;
        }
        else if(energy <1e5)
        {
            energy  =   10*energy;
            multiplier =10;
        }
        else
        {
            multiplier =10;
        }


        hist_name.str("");
        hist_name.clear();
        hist_name2.str("");
        hist_name2.clear();
        hist_title.str("");
        hist_title.clear();
        hist_title2.str("");
        hist_title2.clear();

        hist_name   << "Scattering_dev_"    <<energy;
        hist_name2  << "Scattering_angle_"  <<energy;
        hist_title  << "dev. of mouns with E="<<energy<< " and range(";
        hist_title2 << "angle of mouns with E="<<energy<< " and range(";

        distance    = 0;
        deviation   = 0;
        angle       = 0;

        TH1D *tmp_hist_dev = new TH1D(hist_name.str().c_str(),"dd",80,-2,6);
        TH1D *tmp_hist_angle = new TH1D(hist_name2.str().c_str(),"ss",50,-5,0);

        cout << "Energy: " << energy << endl;
        bool cross = false;
        for(int i =0 ;i< number_of_particles; i++)
        {
            if(i%(number_of_particles/10)==0)
            {
                if(cross)cout<< "+";
                if(!cross)cout<< "#";
                cross = !cross;
                fflush(stdout);
            }
            Particle * prtcl = new Particle("mu",0,0,0,0,0,energy,0);

            prpgtr->Propagate(prtcl);

            distance    += prtcl->GetPropagatedDistance();
            deviation   = sqrt(prtcl->GetX()*prtcl->GetX()+prtcl->GetY()*prtcl->GetY());
            angle       = prtcl->GetTheta();

            tmp_hist_dev->Fill(log10(deviation));
            tmp_hist_angle->Fill(log10(angle));
        }
        hist_title  <<   distance/number_of_particles   <<")";
        hist_title2 <<   distance/number_of_particles   <<")";

        tmp_hist_dev->GetXaxis()->SetTitle("log10( range / [cm] )");
        tmp_hist_dev->GetYaxis()->SetTitle("#");
        tmp_hist_dev->SetTitle(hist_title.str().c_str());
        tmp_hist_dev_vec.push_back(tmp_hist_dev);
        e.push_back(energy);

        tmp_hist_angle->GetXaxis()->SetTitle("log( angle [deg] )");
        tmp_hist_angle->GetYaxis()->SetTitle("#");
        tmp_hist_angle->SetTitle(hist_title2.str().c_str());
        tmp_hist_angle_vec.push_back(tmp_hist_angle);

        cout << endl;
        if(log10(energy)>=4 && (number_of_particles>=1e3))
        {
            number_of_particles/=10;
        }
    }

    TFile *file = new TFile("Scattering_distribution.root","RECREATE");

    TGraphErrors *angle_vs_energy = new TGraphErrors();
    angle_vs_energy->SetName("angle_vs_energy");
    angle_vs_energy->SetTitle("angle vs energy");


    for(unsigned int i =0 ;i< tmp_hist_dev_vec.size();i++)
    {
        angle_vs_energy->SetPoint(i,e.at(i),pow(10 ,tmp_hist_angle_vec.at(i)->GetMean()));
        angle_vs_energy->SetPointError(i,0,pow(10,tmp_hist_angle_vec.at(i)->GetRMS()));
        tmp_hist_dev_vec.at(i)->Write();
        tmp_hist_angle_vec.at(i)->Write();
    }


    TCanvas *can = new TCanvas("can","can",1400,1050);
    can->cd();
    can->SetLogx();
    can->SetGridx();
    can->SetGridy();
    angle_vs_energy->SetMarkerStyle(3);
    angle_vs_energy->Draw("AP");
    angle_vs_energy->GetXaxis()->SetTitle("energy [deg]");
    angle_vs_energy->GetYaxis()->SetTitle("angle [deg]");

    can->Write();
    file->Close();

	return 0;
}
