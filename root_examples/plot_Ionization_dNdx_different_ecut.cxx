#include "TGraph.h"
#include "TCanvas.h"
#include "PROPOSAL/Ionization.h"
#include "PROPOSAL/Particle.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/Medium.h"
#include "TFile.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include <cmath>
#include <sstream>

using namespace std;

int main()
{
    TFile *file     =   new TFile("Ionization_dNdx_different_ecut.root","RECREATE");

    Particle *mu    =   new Particle("mu");
    Particle *e     =   new Particle("e");
    Particle *tau   =   new Particle("tau");

    Medium  *med    =   new Medium("water",1.);

    EnergyCutSettings*  cuts1    =   new EnergyCutSettings(e->GetMass(),-1);
    EnergyCutSettings*  cuts2    =   new EnergyCutSettings(1,-1);
    EnergyCutSettings*  cuts3    =   new EnergyCutSettings(10,-1);
    EnergyCutSettings*  cuts4    =   new EnergyCutSettings(100,-1);
    EnergyCutSettings*  cuts5    =   new EnergyCutSettings(500,-1);
    EnergyCutSettings*  cuts6    =   new EnergyCutSettings(1000,-1);
    EnergyCutSettings*  cuts7    =   new EnergyCutSettings(10000,-1);

    vector<Ionization*> ioniz;

    ioniz.push_back(new Ionization(mu,med,cuts1));
    ioniz.push_back(new Ionization(mu,med,cuts2));
    ioniz.push_back(new Ionization(mu,med,cuts3));
    ioniz.push_back(new Ionization(mu,med,cuts4));
    ioniz.push_back(new Ionization(mu,med,cuts5));
    ioniz.push_back(new Ionization(mu,med,cuts6));
    ioniz.push_back(new Ionization(mu,med,cuts7));

    ioniz.push_back(new Ionization(e,med,cuts1));
    ioniz.push_back(new Ionization(e,med,cuts2));
    ioniz.push_back(new Ionization(e,med,cuts3));
    ioniz.push_back(new Ionization(e,med,cuts4));
    ioniz.push_back(new Ionization(e,med,cuts5));
    ioniz.push_back(new Ionization(e,med,cuts6));
    ioniz.push_back(new Ionization(e,med,cuts7));

    ioniz.push_back(new Ionization(tau,med,cuts1));
    ioniz.push_back(new Ionization(tau,med,cuts2));
    ioniz.push_back(new Ionization(tau,med,cuts3));
    ioniz.push_back(new Ionization(tau,med,cuts4));
    ioniz.push_back(new Ionization(tau,med,cuts5));
    ioniz.push_back(new Ionization(tau,med,cuts6));
    ioniz.push_back(new Ionization(tau,med,cuts7));

    vector<TGraph*> graphs;
    graphs.resize(ioniz.size());


    vector<double> dNdx;
    dNdx.resize(ioniz.size());

    double energy;
    vector<int> counter(graphs.size(),0);

    stringstream graph_name;
    stringstream graph_title;

    vector<TCanvas*> cans;
    cans.resize(graphs.size());

    for(unsigned int i = 0 ; i < ioniz.size() ; i++)
    {
        ioniz.at(i)->EnableDNdxInterpolation();
        graphs.at(i)    =   new TGraph();

        graph_name<<ioniz.at(i)->GetParticle()->GetName()<<"_"<<ioniz.at(i)->GetMedium()->GetName();
        graph_name<<"_"<<ioniz.at(i)->GetEnergyCutSettings()->GetEcut();
        graph_title<<ioniz.at(i)->GetParticle()->GetName()<<" in "<<ioniz.at(i)->GetMedium()->GetName();
        graph_title<<" with ecut = "<<ioniz.at(i)->GetEnergyCutSettings()->GetEcut();

        graphs.at(i)->SetName(graph_name.str().c_str());
        graphs.at(i)->SetTitle(graph_title.str().c_str());

        cans.at(i)      =    new TCanvas(graph_name.str().c_str(),graph_name.str().c_str(),1024,768);

        graph_name.str("");
        graph_name.clear();
        graph_title.str("");
        graph_title.clear();

    }

    for(double log_energy = 0 ;  log_energy < 12; log_energy +=0.2)
    {
        for(unsigned int i = 0 ; i < ioniz.size() ; i++)
        {
            energy  =   ioniz.at(i)->GetParticle()->GetMass()   +   pow(10,log_energy);

            ioniz.at(i)->GetParticle()->SetEnergy(energy);

            dNdx.at(i)  =   ioniz.at(i)->CalculatedNdx();

            if(dNdx.at(i) !=0 )
            {
                graphs.at(i)->SetPoint(counter.at(i),energy,dNdx.at(i));
                counter.at(i)++;
            }
        }

    }

    for(unsigned int i = 0 ; i < graphs.size() ; i++)
    {
        cans.at(i)->cd();
        cans.at(i)->SetLogx();
        cans.at(i)->SetLogy();
        cans.at(i)->SetGridx();
        cans.at(i)->SetGridy();
        graphs.at(i)->SetMarkerStyle(3);
        graphs.at(i)->GetXaxis()->SetTitle("energy [MeV]");
        graphs.at(i)->GetXaxis()->SetTitleOffset(1.1);
        graphs.at(i)->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
        graphs.at(i)->GetYaxis()->SetTitleOffset(1.1);
        graphs.at(i)->Draw("AP");
        cans.at(i)->Write();

    }

    TCanvas* muons      =   new TCanvas("muons","muons",1024,768);
    TCanvas* taus       =   new TCanvas("taus","taus",1024,768);
    TCanvas* electrons  =   new TCanvas("electrons","electrons",1024,768);

    TMultiGraph* muons_gr      =   new TMultiGraph("muons","muons");
    TMultiGraph* taus_gr       =   new TMultiGraph("taus","taus");
    TMultiGraph* electrons_gr  =   new TMultiGraph("electrons","electrons");

    TLegend* muons_leg = new TLegend(0.7,0.1,0.9,0.3);
    muons_leg->SetFillColor(0);

    TLegend* taus_leg = new TLegend(0.7,0.1,0.9,0.3);
    taus_leg->SetFillColor(0);

    TLegend* electrons_leg = new TLegend(0.7,0.1,0.9,0.3);
    electrons_leg->SetFillColor(0);

    stringstream entry_name;

    for(unsigned int i = 0 ; i < graphs.size() ; i++)
    {
        if(ioniz.at(i)->GetParticle()->GetName().compare("mu")==0)
        {
            muons_gr->Add(graphs.at(i),"P");

            entry_name<<"ecut = "<<ioniz.at(i)->GetEnergyCutSettings()->GetEcut();
            muons_leg->AddEntry(graphs.at(i),entry_name.str().c_str(),"p");

            if(muons_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(muons_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(muons_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(muons_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }
            if(muons_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kMagenta);
            }
            if(muons_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kCyan+2);
            }
            if(muons_gr->GetListOfGraphs()->Capacity()==7)
            {
                graphs.at(i)->SetMarkerColor(kOrange);
            }

            entry_name.str("");
            entry_name.clear();
        }
        if(ioniz.at(i)->GetParticle()->GetName().compare("tau")==0)
        {
            taus_gr->Add(graphs.at(i),"P");

            entry_name<<"ecut = "<<ioniz.at(i)->GetEnergyCutSettings()->GetEcut();
            taus_leg->AddEntry(graphs.at(i),entry_name.str().c_str(),"p");

            if(taus_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(taus_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(taus_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(taus_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }
            if(taus_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kMagenta);
            }
            if(taus_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kCyan+2);
            }
            if(taus_gr->GetListOfGraphs()->Capacity()==7)
            {
                graphs.at(i)->SetMarkerColor(kOrange);
            }

            entry_name.str("");
            entry_name.clear();
        }
        if(ioniz.at(i)->GetParticle()->GetName().compare("e")==0)
        {
            electrons_gr->Add(graphs.at(i),"P");

            entry_name<<"ecut = "<<ioniz.at(i)->GetEnergyCutSettings()->GetEcut();
            electrons_leg->AddEntry(graphs.at(i),entry_name.str().c_str(),"p");

            if(electrons_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(electrons_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(electrons_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(electrons_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }
            if(electrons_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kMagenta);
            }
            if(electrons_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kCyan+2);
            }
            if(electrons_gr->GetListOfGraphs()->Capacity()==7)
            {
                graphs.at(i)->SetMarkerColor(kOrange);
            }

            entry_name.str("");
            entry_name.clear();
        }
    }

    muons->cd();
    muons->SetLogx();
    muons->SetLogy();
    muons->SetGridx();
    muons->SetGridy();
    muons_gr->Draw("A");
    muons_gr->GetXaxis()->SetTitle("energy [MeV]");
    muons_gr->GetXaxis()->SetTitleOffset(1.1);
    muons_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    muons_gr->GetYaxis()->SetTitleOffset(1.1);
    muons_leg->Draw("Same");
    muons->Write();

    taus->cd();
    taus->SetLogx();
    taus->SetLogy();
    taus->SetGridx();
    taus->SetGridy();
    taus_gr->Draw("A");
    taus_gr->GetXaxis()->SetTitle("energy [MeV]");
    taus_gr->GetXaxis()->SetTitleOffset(1.1);
    taus_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    taus_gr->GetYaxis()->SetTitleOffset(1.1);
    taus_leg->Draw("Same");
    taus->Write();

    electrons->cd();
    electrons->SetGridx();
    electrons->SetGridy();
    electrons->SetLogx();
    electrons->SetLogy();
    electrons_gr->Draw("A");
    electrons_gr->GetXaxis()->SetTitle("energy [MeV]");
    electrons_gr->GetXaxis()->SetTitleOffset(1.1);
    electrons_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    electrons_gr->GetYaxis()->SetTitleOffset(1.1);
    electrons_leg->Draw("Same");
    electrons->Write();


    file->Close();

    return 0;
}
















