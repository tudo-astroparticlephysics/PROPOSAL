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
        break;
    }

    Gen = new TRandom3(1234);
    cout << "StandardNormal interpolation...\n";
    for(int i=0; i<NSamples ; i++)
    {
        rnd = Gen->Rndm();
        rnd_PropStd_i   = PropStd_i->StandardNormalRandomNumber(rnd,0,sigma,Xlow,Xup,false);
        Hist_PropStd_i->Fill(rnd_PropStd_i);
                break;
    }

    Gen = new TRandom3(1234);
    cout << "Boost erf_inv...\n";
    for(int i=0; i<NSamples ; i++)
    {
        rnd = Gen->Rndm();
        rnd_Boost       = SQRT2*sigma*erfInv( 2*(rnd-0.5) );
        Hist_BoostStd->Fill(rnd_Boost);
                break;
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

    ////////////////////////////////////
    TCanvas* Can_PropStd_i_VGL = new TCanvas("Can_PropStd_i_VGL","Comp",1024,768);
    Can_PropStd_i_VGL->cd();

    TGraph* Graph_PropStd_i = new TGraph();
    TGraph* Graph_BoostStd  = new TGraph();
    cout << "Comparison interpolation...\n";
    int ctr=0;
    double avrg = 841764.74;
    sigma = 3060939.1;
    Xlow = 105.65839;
    Xup = 10000000;
    double xhi,xlo;
    double rndtmp;
    for(double rnd = 1e-5; rnd<1 ; rnd += 1e-5)
    {
        rnd_PropStd_i   = PropStd_i->StandardNormalRandomNumber(rnd,avrg,sigma,Xlow,Xup,false);
        xhi =  0.5+boost::math::erf((Xup-avrg)/(SQRT2*sigma))/2;
        xlo =  0.5+boost::math::erf((Xlow-avrg)/(SQRT2*sigma))/2;
        rndtmp =  xlo + (xhi-xlo)*rnd;
        rnd_Boost       = SQRT2*sigma*erfInv( 2*(rndtmp-0.5) )+avrg;
        Graph_PropStd_i->SetPoint(ctr,rnd,rnd_PropStd_i);
        Graph_BoostStd->SetPoint(ctr,rnd,rnd_Boost);
        ctr++;
    }


    Graph_PropStd_i->GetXaxis()->SetTitle("sigma");
    Graph_PropStd_i->GetYaxis()->SetTitle("#");
    Graph_PropStd_i->SetMarkerColor(kBlue);
    Graph_PropStd_i->SetMarkerStyle(2);

    Graph_BoostStd->GetXaxis()->SetTitle("sigma");
    Graph_BoostStd->GetYaxis()->SetTitle("#");
    Graph_BoostStd->SetLineColor(kRed);
    Graph_BoostStd->SetLineStyle(1);
    Graph_BoostStd->SetLineWidth(3);


    Graph_PropStd_i->Draw("AP");
    Graph_BoostStd->Draw("sameL");

    TLegend* Leg_Graph_Comp = new TLegend(0.1,0.7,0.3,0.9);
    Leg_Graph_Comp->SetFillColor(0);
    Leg_Graph_Comp->AddEntry(Graph_PropStd_i,"PROPOSAL","p");
    Leg_Graph_Comp->AddEntry(Graph_BoostStd,"Boost","l");

    Leg_Graph_Comp->Draw();

    Can_PropStd_i_VGL->Write();


    file->Write();
    file->Close();

    cout << "StandardNormal_Comparison.root created" << endl;
	return 0;
}
