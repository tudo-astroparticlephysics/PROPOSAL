#include "TH1D.h"
#include "TCanvas.h"
#include "TBrowser.h"
#include "TTree.h"
#include "PROPOSAL/Propagator.h"
#include "TFile.h"
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TLegend.h"

#include <boost/math/special_functions/erf.hpp>
#define erfInv(x)   boost::math::erf_inv(x)

using namespace std;

int main()
{
    cout<<"\n----------------------------------------------------------------\n"
        <<"This is an example comparison plot of the implemented class\n"
        <<"StandardNormal and the boost function erf_inv(x).\n"
        <<"A file StandardNormal_Comparison.root will be created.\n"
        <<"The index i means that the resulting value is derived by \n"
        <<"using interpolation of derived sampling points.\n"
        <<"----------------------------------------------------------------\n"
        <<endl;

    double sigma        =2e-3;

    double Nsigma = 3.;
    double theta_max    =Nsigma    *sigma;

    double Xlow = -Nsigma*sigma;
    double Xup  =  Nsigma*sigma;

    int Nbins   =   100/2;


    TRandom3 *Gen;
    StandardNormal* PropStd     = new StandardNormal();
    StandardNormal* PropStd_i   = new StandardNormal();
    PropStd_i->EnableInterpolation("",false);

    TFile *file = new TFile("StandardNormal_Comparison.root","RECREATE");

    TH1D* Hist_PropStd          = new TH1D("Hist_PropStd","",Nbins,Xlow,Xup);
    TH1D* Hist_PropStd_i        = new TH1D("Hist_PropStd_i","",Nbins,Xlow,Xup);
    TH1D* Hist_BoostStd         = new TH1D("Hist_BoostStd","",Nbins,Xlow,Xup);


    double rnd;
    double rnd_Boost,rnd_PropStd,rnd_PropStd_i;

    int NSamples = (int)1e6;

    cout << "Calculating " << NSamples << " points...\n";

    Gen = new TRandom3(1234);
    cout << "StandardNormal...\n";
    for(int i=0; i<NSamples ; i++)
    {
        rnd = Gen->Rndm();
        rnd_PropStd     = PropStd->StandardNormalRandomNumber(rnd,0,sigma,Xlow,Xup,false);
        Hist_PropStd->Fill(rnd_PropStd);
    }

    Gen = new TRandom3(1234);
    cout << "StandardNormal interpolation...\n";
    for(int i=0; i<NSamples ; i++)
    {
        rnd = Gen->Rndm();
        rnd_PropStd_i   = PropStd_i->StandardNormalRandomNumber(rnd,0,sigma,Xlow,Xup,false);
        Hist_PropStd_i->Fill(rnd_PropStd_i);
    }

    Gen = new TRandom3(1234);
    cout << "Boost erf_inv...\n";
    for(int i=0; i<NSamples ; i++)
    {
        rnd = Gen->Rndm();
        rnd_Boost       = SQRT2*sigma*erfInv( 2*(rnd-0.5) );
        Hist_BoostStd->Fill(rnd_Boost);
    }

    cout << "Done!\n";
    Hist_PropStd->GetXaxis()->SetTitle("sigma");
    Hist_PropStd->GetYaxis()->SetTitle("#");
    Hist_PropStd->SetLineColor(kBlue);
    Hist_PropStd->SetLineWidth(2);
    Hist_PropStd->Write();

    Hist_PropStd_i->GetXaxis()->SetTitle("sigma");
    Hist_PropStd_i->GetYaxis()->SetTitle("#");
    Hist_PropStd_i->SetLineColor(kBlue);
    Hist_PropStd_i->SetLineWidth(2);
    Hist_PropStd_i->Write();

    Hist_BoostStd->GetXaxis()->SetTitle("sigma");
    Hist_BoostStd->GetYaxis()->SetTitle("#");
    Hist_BoostStd->SetLineColor(kRed);
    Hist_BoostStd->SetLineWidth(2);
    Hist_BoostStd->Write();


    TCanvas* Can_PropStd = new TCanvas("Can_PropStd","Proposal Std Normal",1024,768);
    Can_PropStd->cd();

    Hist_PropStd->Draw("E");
    Hist_BoostStd->Draw("E Same");

    TLegend* Leg_PropStd = new TLegend(0.1,0.7,0.3,0.9);
    Leg_PropStd->SetFillColor(0);
    Leg_PropStd->AddEntry(Hist_PropStd,"PROPOSAL","l");
    Leg_PropStd->AddEntry(Hist_BoostStd,"Boost","l");

    Leg_PropStd->Draw();

    Can_PropStd->Write();


    TCanvas* Can_PropStd_i = new TCanvas("Can_PropStd_i","Proposal Std Normal",1024,768);
    Can_PropStd_i->cd();

    Hist_PropStd_i->Draw("E");
    Hist_BoostStd->Draw("E Same");

    TLegend* Leg_PropStd_i = new TLegend(0.1,0.7,0.3,0.9);
    Leg_PropStd_i->SetFillColor(0);
    Leg_PropStd_i->AddEntry(Hist_PropStd_i,"PROPOSAL","l");
    Leg_PropStd_i->AddEntry(Hist_BoostStd,"Boost","l");

    Leg_PropStd_i->Draw();

    Can_PropStd_i->Write();

    file->Write();
    file->Close();

    cout << "StandardNormal_Comparison.root created" << endl;
	return 0;
}
