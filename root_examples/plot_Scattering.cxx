#include "TH1D.h"
#include "TCanvas.h"
#include "PROPOSAL/Propagator.h"
#include "TFile.h"
#include <iostream>
#include <sstream>

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


    int number_of_particles =   1e3;

    double distance;
    double deviation;
    double angle;

    double Emin=1e3;
    double Emax=1e12;
    double multiplier = 10;

    for(double energy = Emin; energy <= Emax ; energy *= multiplier)
    {
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
        hist_title  << "dev. of mouns with E="<<energy<< "and range(";
        hist_title2 << "angle of mouns with E="<<energy<< "and range(";

        distance    = 0;
        deviation   = 0;
        angle       = 0;

        TH1D *tmp_hist_dev = new TH1D(hist_name.str().c_str(),"dd",64,1,1000);
        TH1D *tmp_hist_angle = new TH1D(hist_name2.str().c_str(),"ss",64,1e-8,1);
        cout<<tmp_hist_dev<<endl;
        cout<<tmp_hist_angle<<endl;
        for(int i =0 ;i< number_of_particles; i++)
        {
            if(i%(number_of_particles/10)==0)cout<<"Progress: " <<100.*i/number_of_particles<<"%"<<endl;

            Particle * prtcl = new Particle("mu",0,0,0,0,0,energy,0);

            prpgtr->Propagate(prtcl);

            distance    += prtcl->GetPropagatedDistance();
            deviation   = sqrt(prtcl->GetX()*prtcl->GetX()+prtcl->GetY()*prtcl->GetY());
            angle       = prtcl->GetTheta();

            tmp_hist_dev->Fill(deviation);
            tmp_hist_angle->Fill(angle);

            //delete prtcl;
        }
        hist_title  <<   distance/number_of_particles   <<")";
        hist_title2 <<   distance/number_of_particles   <<")";

        tmp_hist_dev->GetXaxis()->SetTitle("range [cm]");
        tmp_hist_dev->GetYaxis()->SetTitle("#");
        tmp_hist_dev->SetTitle(hist_title.str().c_str());
        tmp_hist_dev_vec.push_back(tmp_hist_dev);

        //delete tmp_hist_dev;

        tmp_hist_angle->GetXaxis()->SetTitle("angle [deg]");
        tmp_hist_angle->GetYaxis()->SetTitle("#");
        tmp_hist_angle->SetTitle(hist_title2.str().c_str());
        tmp_hist_angle_vec.push_back(tmp_hist_angle);
        //delete tmp_hist_angle;

    }

    TFile *file = new TFile("Scattering_distribution.root","RECREATE");


    for(unsigned int i =0 ;i< tmp_hist_dev_vec.size();i++)
    {
        tmp_hist_dev_vec.at(i)->Write();
        tmp_hist_angle_vec.at(i)->Write();
    }
    file->Close();

	return 0;
}
