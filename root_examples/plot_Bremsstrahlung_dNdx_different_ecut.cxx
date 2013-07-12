#include "TGraph.h"
#include "TCanvas.h"
#include "PROPOSAL/Bremsstrahlung.h"
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
    TFile *file     =   new TFile("Bremsstrahlung_dNdx_different_ecut.root","RECREATE");

    Particle *mu    =   new Particle("mu");
    Particle *tau   =   new Particle("tau");
    Particle *e     =   new Particle("e");

    Medium  *med   =   new Medium("water",1.);

    EnergyCutSettings*  cuts1    =   new EnergyCutSettings(0.1,-1);
    EnergyCutSettings*  cuts2    =   new EnergyCutSettings(1,-1);
    EnergyCutSettings*  cuts3    =   new EnergyCutSettings(10,-1);
    EnergyCutSettings*  cuts4    =   new EnergyCutSettings(100,-1);
    EnergyCutSettings*  cuts5    =   new EnergyCutSettings(500,-1);
    EnergyCutSettings*  cuts6    =   new EnergyCutSettings(1000,-1);
    EnergyCutSettings*  cuts7    =   new EnergyCutSettings(10000,-1);

    vector<Bremsstrahlung*> brems;

    brems.push_back(new Bremsstrahlung(mu,med,cuts1));
    brems.push_back(new Bremsstrahlung(tau,med,cuts1));
    brems.push_back(new Bremsstrahlung(e,med,cuts1));

    brems.push_back(new Bremsstrahlung(mu,med,cuts2));
    brems.push_back(new Bremsstrahlung(tau,med,cuts2));
    brems.push_back(new Bremsstrahlung(e,med,cuts2));

    brems.push_back(new Bremsstrahlung(mu,med,cuts3));
    brems.push_back(new Bremsstrahlung(tau,med,cuts3));
    brems.push_back(new Bremsstrahlung(e,med,cuts3));

    brems.push_back(new Bremsstrahlung(mu,med,cuts4));
    brems.push_back(new Bremsstrahlung(tau,med,cuts4));
    brems.push_back(new Bremsstrahlung(e,med,cuts4));

    brems.push_back(new Bremsstrahlung(mu,med,cuts5));
    brems.push_back(new Bremsstrahlung(tau,med,cuts5));
    brems.push_back(new Bremsstrahlung(e,med,cuts5));

    brems.push_back(new Bremsstrahlung(mu,med,cuts6));
    brems.push_back(new Bremsstrahlung(tau,med,cuts6));
    brems.push_back(new Bremsstrahlung(e,med,cuts6));

    brems.push_back(new Bremsstrahlung(mu,med,cuts7));
    brems.push_back(new Bremsstrahlung(tau,med,cuts7));
    brems.push_back(new Bremsstrahlung(e,med,cuts7));


    brems.push_back(new Bremsstrahlung(mu,med,cuts1));
    brems.push_back(new Bremsstrahlung(tau,med,cuts1));
    brems.push_back(new Bremsstrahlung(e,med,cuts1));

    brems.push_back(new Bremsstrahlung(mu,med,cuts2));
    brems.push_back(new Bremsstrahlung(tau,med,cuts2));
    brems.push_back(new Bremsstrahlung(e,med,cuts2));

    brems.push_back(new Bremsstrahlung(mu,med,cuts3));
    brems.push_back(new Bremsstrahlung(tau,med,cuts3));
    brems.push_back(new Bremsstrahlung(e,med,cuts3));

    brems.push_back(new Bremsstrahlung(mu,med,cuts4));
    brems.push_back(new Bremsstrahlung(tau,med,cuts4));
    brems.push_back(new Bremsstrahlung(e,med,cuts4));

    brems.push_back(new Bremsstrahlung(mu,med,cuts5));
    brems.push_back(new Bremsstrahlung(tau,med,cuts5));
    brems.push_back(new Bremsstrahlung(e,med,cuts5));

    brems.push_back(new Bremsstrahlung(mu,med,cuts6));
    brems.push_back(new Bremsstrahlung(tau,med,cuts6));
    brems.push_back(new Bremsstrahlung(e,med,cuts6));

    brems.push_back(new Bremsstrahlung(mu,med,cuts7));
    brems.push_back(new Bremsstrahlung(tau,med,cuts7));
    brems.push_back(new Bremsstrahlung(e,med,cuts7));


    brems.push_back(new Bremsstrahlung(mu,med,cuts1));
    brems.push_back(new Bremsstrahlung(tau,med,cuts1));
    brems.push_back(new Bremsstrahlung(e,med,cuts1));

    brems.push_back(new Bremsstrahlung(mu,med,cuts2));
    brems.push_back(new Bremsstrahlung(tau,med,cuts2));
    brems.push_back(new Bremsstrahlung(e,med,cuts2));

    brems.push_back(new Bremsstrahlung(mu,med,cuts3));
    brems.push_back(new Bremsstrahlung(tau,med,cuts3));
    brems.push_back(new Bremsstrahlung(e,med,cuts3));

    brems.push_back(new Bremsstrahlung(mu,med,cuts4));
    brems.push_back(new Bremsstrahlung(tau,med,cuts4));
    brems.push_back(new Bremsstrahlung(e,med,cuts4));

    brems.push_back(new Bremsstrahlung(mu,med,cuts5));
    brems.push_back(new Bremsstrahlung(tau,med,cuts5));
    brems.push_back(new Bremsstrahlung(e,med,cuts5));

    brems.push_back(new Bremsstrahlung(mu,med,cuts6));
    brems.push_back(new Bremsstrahlung(tau,med,cuts6));
    brems.push_back(new Bremsstrahlung(e,med,cuts6));

    brems.push_back(new Bremsstrahlung(mu,med,cuts7));
    brems.push_back(new Bremsstrahlung(tau,med,cuts7));
    brems.push_back(new Bremsstrahlung(e,med,cuts7));


    brems.push_back(new Bremsstrahlung(mu,med,cuts1));
    brems.push_back(new Bremsstrahlung(tau,med,cuts1));
    brems.push_back(new Bremsstrahlung(e,med,cuts1));

    brems.push_back(new Bremsstrahlung(mu,med,cuts2));
    brems.push_back(new Bremsstrahlung(tau,med,cuts2));
    brems.push_back(new Bremsstrahlung(e,med,cuts2));

    brems.push_back(new Bremsstrahlung(mu,med,cuts3));
    brems.push_back(new Bremsstrahlung(tau,med,cuts3));
    brems.push_back(new Bremsstrahlung(e,med,cuts3));

    brems.push_back(new Bremsstrahlung(mu,med,cuts4));
    brems.push_back(new Bremsstrahlung(tau,med,cuts4));
    brems.push_back(new Bremsstrahlung(e,med,cuts4));

    brems.push_back(new Bremsstrahlung(mu,med,cuts5));
    brems.push_back(new Bremsstrahlung(tau,med,cuts5));
    brems.push_back(new Bremsstrahlung(e,med,cuts5));

    brems.push_back(new Bremsstrahlung(mu,med,cuts6));
    brems.push_back(new Bremsstrahlung(tau,med,cuts6));
    brems.push_back(new Bremsstrahlung(e,med,cuts6));

    brems.push_back(new Bremsstrahlung(mu,med,cuts7));
    brems.push_back(new Bremsstrahlung(tau,med,cuts7));
    brems.push_back(new Bremsstrahlung(e,med,cuts7));

    for(unsigned int i = 0; i<brems.size() ; i++)
    {
        if(i<21) brems.at(i)->SetParametrization(1);
        else if(i<42) brems.at(i)->SetParametrization(2);
        else if(i<63) brems.at(i)->SetParametrization(3);
        else if(i<84) brems.at(i)->SetParametrization(4);

    }

    vector<TGraph*> graphs;
    graphs.resize(brems.size());


    vector<double> dNdx;
    dNdx.resize(brems.size());

    double energy;

    vector<int> counter(graphs.size(),0);

    stringstream graph_name;
    stringstream graph_title;

    vector<TCanvas*> cans;
    cans.resize(graphs.size());

    for(unsigned int i = 0 ; i < brems.size() ; i++)
    {
        brems.at(i)->EnableDNdxInterpolation();
        graphs.at(i)    =   new TGraph();

        graph_name<<brems.at(i)->GetParticle()->GetName()<<"_"<<brems.at(i)->GetMedium()->GetName();
        graph_name<<"_"<<brems.at(i)->GetParametrization();
        graph_name<<"_"<<brems.at(i)->GetEnergyCutSettings()->GetEcut();

        graph_title<<brems.at(i)->GetParticle()->GetName()<<" in "<<brems.at(i)->GetMedium()->GetName();
        graph_title<<" with ecut = "<<brems.at(i)->GetEnergyCutSettings()->GetEcut();

        if(brems.at(i)->GetParametrization()==1)
        {
            graph_title<<", parametrization: Kelner-Kakoulin-Petrukhin";
            graphs.at(i)->SetMarkerStyle(4);

        }
        if(brems.at(i)->GetParametrization()==2)
        {
            graph_title<<", parametrization: Andreev-Bezrukov-Bugaev";
            graphs.at(i)->SetMarkerStyle(3);

        }
        if(brems.at(i)->GetParametrization()==3)
        {
            graph_title<<", parametrization: Petrukhin-Shestakov";
            graphs.at(i)->SetMarkerStyle(2);

        }
        if(brems.at(i)->GetParametrization()==4)
        {
            graph_title<<", parametrization: Complete screening case";
            graphs.at(i)->SetMarkerStyle(28);

        }


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
        for(unsigned int i = 0 ; i < brems.size() ; i++)
        {
            energy  =   brems.at(i)->GetParticle()->GetMass()   +   pow(10,log_energy);

            brems.at(i)->GetParticle()->SetEnergy(energy);

            dNdx.at(i)  =   brems.at(i)->CalculatedNdx();

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
        graphs.at(i)->GetXaxis()->SetTitle("energy [MeV]");
        graphs.at(i)->GetXaxis()->SetTitleOffset(1.1);
        graphs.at(i)->GetYaxis()->SetTitle("dN/dx [ cm^{-1}]");
        graphs.at(i)->Draw("AP");
        cans.at(i)->Write();

    }

    TCanvas* muons_para1      =   new TCanvas("muons_para1","muons in water, parametrization = KKP",1024,768);
    TCanvas* taus_para1       =   new TCanvas("taus_para1","taus in water, parametrization = KKP",1024,768);
    TCanvas* electrons_para1  =   new TCanvas("electrons_para1","electrons in water, parametrization = KKP",1024,768);

    TMultiGraph* muons_para1_gr      =   new TMultiGraph("muons_para1","muons in water, parametrization = KKP");
    TMultiGraph* taus_para1_gr       =   new TMultiGraph("taus_para1","taus in water, parametrization = KKP");
    TMultiGraph* electrons_para1_gr  =   new TMultiGraph("electrons_para1","electrons in water, parametrization = KKP");


    TCanvas* muons_para2      =   new TCanvas("muons_para2","muons in water, parametrization = ABB",1024,768);
    TCanvas* taus_para2       =   new TCanvas("taus_para2","taus in water, parametrization = ABB",1024,768);
    TCanvas* electrons_para2  =   new TCanvas("electrons_para2","electrons in water, parametrization = ABB",1024,768);

    TMultiGraph* muons_para2_gr      =   new TMultiGraph("muons_para2","muons in water, parametrization = ABB");
    TMultiGraph* taus_para2_gr       =   new TMultiGraph("taus_para2","taus in water, parametrization = ABB");
    TMultiGraph* electrons_para2_gr  =   new TMultiGraph("electrons_para2","electrons in water, parametrization = ABB");


    TCanvas* muons_para3      =   new TCanvas("muons_para3","muons in water, parametrization = PS",1024,768);
    TCanvas* taus_para3       =   new TCanvas("taus_para3","taus in water, parametrization = PS",1024,768);
    TCanvas* electrons_para3  =   new TCanvas("electrons_para3","electrons in water, parametrization = PS",1024,768);

    TMultiGraph* muons_para3_gr      =   new TMultiGraph("muons_para3","muons in water, parametrization = PS");
    TMultiGraph* taus_para3_gr       =   new TMultiGraph("taus_para3","taus in water, parametrization = PS");
    TMultiGraph* electrons_para3_gr  =   new TMultiGraph("electrons_para3","electrons in water, parametrization = PS");


    TCanvas* muons_para4      =   new TCanvas("muons_para4","muons in water, parametrization = CSC",1024,768);
    TCanvas* taus_para4       =   new TCanvas("taus_para4","taus in water, parametrization = CSC",1024,768);
    TCanvas* electrons_para4  =   new TCanvas("electrons_para4","electrons in water, parametrization = CSC",1024,768);

    TMultiGraph* muons_para4_gr      =   new TMultiGraph("muons_para4","muons in water, parametrization = CSC");
    TMultiGraph* taus_para4_gr       =   new TMultiGraph("taus_para4","taus in water, parametrization = CSC");
    TMultiGraph* electrons_para4_gr  =   new TMultiGraph("electrons_para4","electrons in water, parametrization = CSC");

    TLegend* muons_leg_para1 = new TLegend(0.7,0.2,0.9,0.4);
    muons_leg_para1->SetFillColor(0);

    TLegend* muons_leg_para2 = new TLegend(0.7,0.2,0.9,0.4);
    muons_leg_para2->SetFillColor(0);

    TLegend* muons_leg_para3 = new TLegend(0.7,0.2,0.9,0.4);
    muons_leg_para3->SetFillColor(0);

    TLegend* muons_leg_para4 = new TLegend(0.7,0.2,0.9,0.4);
    muons_leg_para4->SetFillColor(0);


    TLegend* taus_leg_para1 = new TLegend(0.7,0.2,0.9,0.4);
    taus_leg_para1->SetFillColor(0);

    TLegend* taus_leg_para2 = new TLegend(0.7,0.2,0.9,0.4);
    taus_leg_para2->SetFillColor(0);

    TLegend* taus_leg_para3 = new TLegend(0.7,0.2,0.9,0.4);
    taus_leg_para3->SetFillColor(0);

    TLegend* taus_leg_para4 = new TLegend(0.7,0.2,0.9,0.4);
    taus_leg_para4->SetFillColor(0);


    TLegend* electrons_leg_para1 = new TLegend(0.7,0.2,0.9,0.4);
    electrons_leg_para1->SetFillColor(0);

    TLegend* electrons_leg_para2 = new TLegend(0.7,0.2,0.9,0.4);
    electrons_leg_para2->SetFillColor(0);

    TLegend* electrons_leg_para3 = new TLegend(0.7,0.2,0.9,0.4);
    electrons_leg_para3->SetFillColor(0);

    TLegend* electrons_leg_para4 = new TLegend(0.7,0.2,0.9,0.4);
    electrons_leg_para4->SetFillColor(0);


    stringstream leg_entry;

    for(unsigned int i = 0 ; i < graphs.size() ; i++)
    {
        if(brems.at(i)->GetParticle()->GetName().compare("mu")==0 && brems.at(i)->GetParametrization()==1)
        {
            muons_para1_gr->Add(graphs.at(i),"P");

            leg_entry<<"ecut = "<<brems.at(i)->GetEnergyCutSettings()->GetEcut();
            muons_leg_para1->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(muons_para1_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(muons_para1_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(muons_para1_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(muons_para1_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }
            if(muons_para1_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kMagenta);
            }
            if(muons_para1_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kCyan+2);
            }
            if(muons_para1_gr->GetListOfGraphs()->Capacity()==7)
            {
                graphs.at(i)->SetMarkerColor(kOrange);
            }

            leg_entry.str("");
            leg_entry.clear();

        }
        if(brems.at(i)->GetParticle()->GetName().compare("mu")==0 && brems.at(i)->GetParametrization()==2)
        {
            muons_para2_gr->Add(graphs.at(i),"P");

            leg_entry<<"ecut = "<<brems.at(i)->GetEnergyCutSettings()->GetEcut();
            muons_leg_para2->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(muons_para2_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(muons_para2_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(muons_para2_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(muons_para2_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }
            if(muons_para2_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kMagenta);
            }
            if(muons_para2_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kCyan+2);
            }
            if(muons_para2_gr->GetListOfGraphs()->Capacity()==7)
            {
                graphs.at(i)->SetMarkerColor(kOrange);
            }

            leg_entry.str("");
            leg_entry.clear();

        }
        if(brems.at(i)->GetParticle()->GetName().compare("mu")==0 && brems.at(i)->GetParametrization()==3)
        {
            muons_para3_gr->Add(graphs.at(i),"P");

            leg_entry<<"ecut = "<<brems.at(i)->GetEnergyCutSettings()->GetEcut();
            muons_leg_para3->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(muons_para3_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(muons_para3_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(muons_para3_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(muons_para3_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }
            if(muons_para3_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kMagenta);
            }
            if(muons_para3_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kCyan+2);
            }
            if(muons_para3_gr->GetListOfGraphs()->Capacity()==7)
            {
                graphs.at(i)->SetMarkerColor(kOrange);
            }

            leg_entry.str("");
            leg_entry.clear();

        }
        if(brems.at(i)->GetParticle()->GetName().compare("mu")==0 && brems.at(i)->GetParametrization()==4)
        {
            muons_para4_gr->Add(graphs.at(i),"P");

            leg_entry<<"ecut = "<<brems.at(i)->GetEnergyCutSettings()->GetEcut();
            muons_leg_para4->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(muons_para4_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(muons_para4_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(muons_para4_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(muons_para4_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }
            if(muons_para4_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kMagenta);
            }
            if(muons_para4_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kCyan+2);
            }
            if(muons_para4_gr->GetListOfGraphs()->Capacity()==7)
            {
                graphs.at(i)->SetMarkerColor(kOrange);
            }

            leg_entry.str("");
            leg_entry.clear();

        }



        if(brems.at(i)->GetParticle()->GetName().compare("tau")==0 && brems.at(i)->GetParametrization()==1)
        {
            taus_para1_gr->Add(graphs.at(i),"P");

            leg_entry<<"ecut = "<<brems.at(i)->GetEnergyCutSettings()->GetEcut();
            taus_leg_para1->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(taus_para1_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(taus_para1_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(taus_para1_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(taus_para1_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }
            if(taus_para1_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kMagenta);
            }
            if(taus_para1_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kCyan+2);
            }
            if(taus_para1_gr->GetListOfGraphs()->Capacity()==7)
            {
                graphs.at(i)->SetMarkerColor(kOrange);
            }

            leg_entry.str("");
            leg_entry.clear();

        }
        if(brems.at(i)->GetParticle()->GetName().compare("tau")==0 && brems.at(i)->GetParametrization()==2)
        {
            taus_para2_gr->Add(graphs.at(i),"P");

            leg_entry<<"ecut = "<<brems.at(i)->GetEnergyCutSettings()->GetEcut();
            taus_leg_para2->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(taus_para2_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(taus_para2_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(taus_para2_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(taus_para2_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }
            if(taus_para2_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kMagenta);
            }
            if(taus_para2_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kCyan+2);
            }
            if(taus_para2_gr->GetListOfGraphs()->Capacity()==7)
            {
                graphs.at(i)->SetMarkerColor(kOrange);
            }

            leg_entry.str("");
            leg_entry.clear();

        }
        if(brems.at(i)->GetParticle()->GetName().compare("tau")==0 && brems.at(i)->GetParametrization()==3)
        {
            taus_para3_gr->Add(graphs.at(i),"P");

            leg_entry<<"ecut = "<<brems.at(i)->GetEnergyCutSettings()->GetEcut();
            taus_leg_para3->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(taus_para3_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(taus_para3_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(taus_para3_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(taus_para3_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }
            if(taus_para3_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kMagenta);
            }
            if(taus_para3_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kCyan+2);
            }
            if(taus_para3_gr->GetListOfGraphs()->Capacity()==7)
            {
                graphs.at(i)->SetMarkerColor(kOrange);
            }

            leg_entry.str("");
            leg_entry.clear();

        }
        if(brems.at(i)->GetParticle()->GetName().compare("tau")==0 && brems.at(i)->GetParametrization()==4)
        {
            taus_para4_gr->Add(graphs.at(i),"P");

            leg_entry<<"ecut = "<<brems.at(i)->GetEnergyCutSettings()->GetEcut();
            taus_leg_para4->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(taus_para4_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(taus_para4_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(taus_para4_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(taus_para4_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }
            if(taus_para4_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kMagenta);
            }
            if(taus_para4_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kCyan+2);
            }
            if(taus_para4_gr->GetListOfGraphs()->Capacity()==7)
            {
                graphs.at(i)->SetMarkerColor(kOrange);
            }

            leg_entry.str("");
            leg_entry.clear();

        }


        if(brems.at(i)->GetParticle()->GetName().compare("e")==0 && brems.at(i)->GetParametrization()==1)
        {
            electrons_para1_gr->Add(graphs.at(i),"P");

            leg_entry<<"ecut = "<<brems.at(i)->GetEnergyCutSettings()->GetEcut();
            electrons_leg_para1->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(electrons_para1_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(electrons_para1_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(electrons_para1_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(electrons_para1_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }
            if(electrons_para1_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kMagenta);
            }
            if(electrons_para1_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kCyan+2);
            }
            if(electrons_para1_gr->GetListOfGraphs()->Capacity()==7)
            {
                graphs.at(i)->SetMarkerColor(kOrange);
            }

            leg_entry.str("");
            leg_entry.clear();

        }
        if(brems.at(i)->GetParticle()->GetName().compare("e")==0 && brems.at(i)->GetParametrization()==2)
        {
            electrons_para2_gr->Add(graphs.at(i),"P");

            leg_entry<<"ecut = "<<brems.at(i)->GetEnergyCutSettings()->GetEcut();
            electrons_leg_para2->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(electrons_para2_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(electrons_para2_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(electrons_para2_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(electrons_para2_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }
            if(electrons_para2_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kMagenta);
            }
            if(electrons_para2_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kCyan+2);
            }
            if(electrons_para2_gr->GetListOfGraphs()->Capacity()==7)
            {
                graphs.at(i)->SetMarkerColor(kOrange);
            }

            leg_entry.str("");
            leg_entry.clear();

        }
        if(brems.at(i)->GetParticle()->GetName().compare("e")==0 && brems.at(i)->GetParametrization()==3)
        {
            electrons_para3_gr->Add(graphs.at(i),"P");

            leg_entry<<"ecut = "<<brems.at(i)->GetEnergyCutSettings()->GetEcut();
            electrons_leg_para3->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(electrons_para3_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(electrons_para3_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(electrons_para3_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(electrons_para3_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }
            if(electrons_para3_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kMagenta);
            }
            if(electrons_para3_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kCyan+2);
            }
            if(electrons_para3_gr->GetListOfGraphs()->Capacity()==7)
            {
                graphs.at(i)->SetMarkerColor(kOrange);
            }

            leg_entry.str("");
            leg_entry.clear();

        }
        if(brems.at(i)->GetParticle()->GetName().compare("e")==0 && brems.at(i)->GetParametrization()==4)
        {
            electrons_para4_gr->Add(graphs.at(i),"P");

            leg_entry<<"ecut = "<<brems.at(i)->GetEnergyCutSettings()->GetEcut();
            electrons_leg_para4->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(electrons_para4_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(electrons_para4_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(electrons_para4_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(electrons_para4_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }
            if(electrons_para4_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kMagenta);
            }
            if(electrons_para4_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kCyan+2);
            }
            if(electrons_para4_gr->GetListOfGraphs()->Capacity()==7)
            {
                graphs.at(i)->SetMarkerColor(kOrange);
            }

            leg_entry.str("");
            leg_entry.clear();

        }

    }

    muons_para1->cd();
    muons_para1->SetLogx();
    muons_para1->SetLogy();
    muons_para1->SetGridx();
    muons_para1->SetGridy();
    muons_para1_gr->Draw("A");
    muons_para1_gr->GetXaxis()->SetTitle("energy [MeV]");
    muons_para1_gr->GetXaxis()->SetTitleOffset(1.1);
    muons_para1_gr->GetYaxis()->SetTitle("dN/dx [ cm^{-1}]");
    muons_leg_para1->Draw("Same");
    muons_para1->Write();


    muons_para2->cd();
    muons_para2->SetLogx();
    muons_para2->SetLogy();
    muons_para2->SetGridx();
    muons_para2->SetGridy();
    muons_para2_gr->Draw("A");
    muons_para2_gr->GetXaxis()->SetTitle("energy [MeV]");
    muons_para2_gr->GetXaxis()->SetTitleOffset(1.1);
    muons_para2_gr->GetYaxis()->SetTitle("dN/dx [ cm^{-1}]");
    muons_leg_para2->Draw("Same");
    muons_para2->Write();


    muons_para3->cd();
    muons_para3->SetLogx();
    muons_para3->SetLogy();
    muons_para3->SetGridx();
    muons_para3->SetGridy();
    muons_para3_gr->Draw("A");
    muons_para3_gr->GetXaxis()->SetTitle("energy [MeV]");
    muons_para3_gr->GetXaxis()->SetTitleOffset(1.1);
    muons_para3_gr->GetYaxis()->SetTitle("dN/dx [ cm^{-1}]");
    muons_leg_para3->Draw("Same");
    muons_para3->Write();


    muons_para4->cd();
    muons_para4->SetLogx();
    muons_para4->SetLogy();
    muons_para4->SetGridx();
    muons_para4->SetGridy();
    muons_para4_gr->Draw("A");
    muons_para4_gr->GetXaxis()->SetTitle("energy [MeV]");
    muons_para4_gr->GetXaxis()->SetTitleOffset(1.1);
    muons_para4_gr->GetYaxis()->SetTitle("dN/dx [ cm^{-1}]");
    muons_leg_para4->Draw("Same");
    muons_para4->Write();


    taus_para1->cd();
    taus_para1->SetLogx();
    taus_para1->SetLogy();
    taus_para1->SetGridx();
    taus_para1->SetGridy();
    taus_para1_gr->Draw("A");
    taus_para1_gr->GetXaxis()->SetTitle("energy [MeV]");
    taus_para1_gr->GetXaxis()->SetTitleOffset(1.1);
    taus_para1_gr->GetYaxis()->SetTitle("dN/dx [ cm^{-1}]");
    taus_leg_para1->Draw("Same");
    taus_para1->Write();


    taus_para2->cd();
    taus_para2->SetLogx();
    taus_para2->SetLogy();
    taus_para2->SetGridx();
    taus_para2->SetGridy();
    taus_para2_gr->Draw("A");
    taus_para2_gr->GetXaxis()->SetTitle("energy [MeV]");
    taus_para2_gr->GetXaxis()->SetTitleOffset(1.1);
    taus_para2_gr->GetYaxis()->SetTitle("dN/dx [ cm^{-1}]");
    taus_leg_para2->Draw("Same");
    taus_para2->Write();


    taus_para3->cd();
    taus_para3->SetLogx();
    taus_para3->SetLogy();
    taus_para3->SetGridx();
    taus_para3->SetGridy();
    taus_para3_gr->Draw("A");
    taus_para3_gr->GetXaxis()->SetTitle("energy [MeV]");
    taus_para3_gr->GetXaxis()->SetTitleOffset(1.1);
    taus_para3_gr->GetYaxis()->SetTitle("dN/dx [ cm^{-1}]");
    taus_leg_para3->Draw("Same");
    taus_para3->Write();


    taus_para4->cd();
    taus_para4->SetLogx();
    taus_para4->SetLogy();
    taus_para4->SetGridx();
    taus_para4->SetGridy();
    taus_para4_gr->Draw("A");
    taus_para4_gr->GetXaxis()->SetTitle("energy [MeV]");
    taus_para4_gr->GetXaxis()->SetTitleOffset(1.1);
    taus_para4_gr->GetYaxis()->SetTitle("dN/dx [ cm^{-1}]");
    taus_leg_para4->Draw("Same");
    taus_para4->Write();


    electrons_para1->cd();
    electrons_para1->SetLogx();
    electrons_para1->SetLogy();
    electrons_para1->SetGridx();
    electrons_para1->SetGridy();
    electrons_para1_gr->Draw("A");
    electrons_para1_gr->GetXaxis()->SetTitle("energy [MeV]");
    electrons_para1_gr->GetXaxis()->SetTitleOffset(1.1);
    electrons_para1_gr->GetYaxis()->SetTitle("dN/dx [ cm^{-1}]");
    electrons_leg_para1->Draw("Same");
    electrons_para1->Write();


    electrons_para2->cd();
    electrons_para2->SetLogx();
    electrons_para2->SetLogy();
    electrons_para2->SetGridx();
    electrons_para2->SetGridy();
    electrons_para2_gr->Draw("A");
    electrons_para2_gr->GetXaxis()->SetTitle("energy [MeV]");
    electrons_para2_gr->GetXaxis()->SetTitleOffset(1.1);
    electrons_para2_gr->GetYaxis()->SetTitle("dN/dx [ cm^{-1}]");
    electrons_leg_para2->Draw("Same");
    electrons_para2->Write();


    electrons_para3->cd();
    electrons_para3->SetLogx();
    electrons_para3->SetLogy();
    electrons_para3->SetGridx();
    electrons_para3->SetGridy();
    electrons_para3_gr->Draw("A");
    electrons_para3_gr->GetXaxis()->SetTitle("energy [MeV]");
    electrons_para3_gr->GetXaxis()->SetTitleOffset(1.1);
    electrons_para3_gr->GetYaxis()->SetTitle("dN/dx [ cm^{-1}]");
    electrons_leg_para3->Draw("Same");
    electrons_para3->Write();


    electrons_para4->cd();
    electrons_para4->SetLogx();
    electrons_para4->SetLogy();
    electrons_para4->SetGridx();
    electrons_para4->SetGridy();
    electrons_para4_gr->Draw("A");
    electrons_para4_gr->GetXaxis()->SetTitle("energy [MeV]");
    electrons_para4_gr->GetXaxis()->SetTitleOffset(1.1);
    electrons_para4_gr->GetYaxis()->SetTitle("dN/dx [ cm^{-1}]");
    electrons_leg_para4->Draw("Same");
    electrons_para4->Write();


    file->Close();

    return 0;
}
















