/** \file
 * This is an example to plot dN/dx for Epairproduction.\n
 * If you have ROOT installed this file will be compilied when running make\n
 * and you can find an executable in your build directory in root_examples:\n
 * -->   Plot_Epairproduction_dNdx   <--\n
 * A ROOT-file named Epairproduction_dNdx.root will be generated,
 * which contains TCanvas's with TGraphs dN/dx vs. energy.\n
 * These TGraphs are generated for muons, taus and electrons in \n
 * water, hydrogen and uranium \n
 * (with and without Landau-Pomeranchuk-Migdal supression)\n.
 * You will also find summary plots for every particle and medium.
 * @brief Example to plot dN/dx for Epairproduction
 * @author Jan-Hendrik Köhne
 */

#include "TGraph.h"
#include "TCanvas.h"
#include "PROPOSAL/Epairproduction.h"
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

    cout<<"\n----------------------------------------------------------------\n"
        <<"This is an example to plot dN/dx for Epairproduction.\n"
        <<"A ROOT-file named Epairproduction_dNdx.root will be generated,\n"
        <<"which contains TCanvas's with TGraphs dN/dx vs. energy.\n"
        <<"These TGraphs are generated for muons, taus and electrons in \n"
        <<"water, hydrogen and uranium \n"
        <<"(with and without Landau-Pomeranchuk-Migdal supression).\n"
        <<"You will also find summary plots for every particle and medium.\n"
        <<"----------------------------------------------------------------\n"
        <<endl;

    TFile *file     =   new TFile("Epairproduction_dNdx.root","RECREATE");

    Particle *mu    =   new Particle("mu");
    Particle *tau   =   new Particle("tau");
    Particle *e     =   new Particle("e");

    Medium  *med1   =   new Medium("hydrogen",1.);
    Medium  *med2   =   new Medium("water",1.);
    Medium  *med3   =   new Medium("uranium",1.);

    EnergyCutSettings*  cuts    =   new EnergyCutSettings(2*e->GetMass(),-1);

    vector<Epairproduction*> epair;

    epair.push_back(new Epairproduction(mu,med1,cuts));
    epair.push_back(new Epairproduction(mu,med2,cuts));
    epair.push_back(new Epairproduction(mu,med3,cuts));
    epair.push_back(new Epairproduction(tau,med1,cuts));
    epair.push_back(new Epairproduction(tau,med2,cuts));
    epair.push_back(new Epairproduction(tau,med3,cuts));
    epair.push_back(new Epairproduction(e,med1,cuts));
    epair.push_back(new Epairproduction(e,med2,cuts));
    epair.push_back(new Epairproduction(e,med3,cuts));

    epair.push_back(new Epairproduction(mu,med1,cuts));
    epair.at(epair.size()-1)->EnableLpmEffect(true);
    epair.push_back(new Epairproduction(mu,med2,cuts));
    epair.at(epair.size()-1)->EnableLpmEffect(true);
    epair.push_back(new Epairproduction(mu,med3,cuts));
    epair.at(epair.size()-1)->EnableLpmEffect(true);
    epair.push_back(new Epairproduction(tau,med1,cuts));
    epair.at(epair.size()-1)->EnableLpmEffect(true);
    epair.push_back(new Epairproduction(tau,med2,cuts));
    epair.at(epair.size()-1)->EnableLpmEffect(true);
    epair.push_back(new Epairproduction(tau,med3,cuts));
    epair.at(epair.size()-1)->EnableLpmEffect(true);
    epair.push_back(new Epairproduction(e,med1,cuts));
    epair.at(epair.size()-1)->EnableLpmEffect(true);
    epair.push_back(new Epairproduction(e,med2,cuts));
    epair.at(epair.size()-1)->EnableLpmEffect(true);
    epair.push_back(new Epairproduction(e,med3,cuts));
    epair.at(epair.size()-1)->EnableLpmEffect(true);

    vector<TGraph*> graphs;
    graphs.resize(epair.size());


    vector<double> dNdx;
    dNdx.resize(epair.size());

    double energy;

    vector<int> counter(graphs.size(),0);
    stringstream graph_name;
    stringstream graph_title;

    vector<TCanvas*> cans;
    cans.resize(graphs.size());

    vector<TCanvas*> cans_with_and_without;
    cans_with_and_without.resize(graphs.size()/2);

    for(unsigned int i = 0 ; i < epair.size() ; i++)
    {
        epair.at(i)->EnableDNdxInterpolation();
        graphs.at(i)    =   new TGraph();

        graph_name<<epair.at(i)->GetParticle()->GetName()<<"_"<<epair.at(i)->GetMedium()->GetName();
        graph_title<<epair.at(i)->GetParticle()->GetName()<<" in "<<epair.at(i)->GetMedium()->GetName();

        if(i<9)
            cans_with_and_without.at(i)      =    new TCanvas(graph_name.str().c_str(),graph_name.str().c_str(),1024,768);

        if(epair.at(i)->GetLpmEffectEnabled())
        {
            graph_name<<"_lpm";
            graph_title<<" with lpm effect";
            graphs.at(i)->SetMarkerStyle(4);

        }
        else
        {
            graph_name<<"_nolpm";
            graph_title<<" without lpm effect";
            graphs.at(i)->SetMarkerStyle(3);

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
        for(unsigned int i = 0 ; i < epair.size() ; i++)
        {
            energy  =   epair.at(i)->GetParticle()->GetMass()   +   pow(10,log_energy);

            epair.at(i)->GetParticle()->SetEnergy(energy);

            dNdx.at(i)  =   epair.at(i)->CalculatedNdx();

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
        graphs.at(i)->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
        graphs.at(i)->GetYaxis()->SetTitleOffset(1.1);
        graphs.at(i)->Draw("AP");
        cans.at(i)->Write();

    }

    TLegend* leg;

    for(unsigned int i = 0 ; i < graphs.size()/2 ; i++)
    {
        leg = new TLegend(0.7,0.2,0.9,0.4);
        leg->SetFillColor(0);
        cans_with_and_without.at(i)->cd();
        cans_with_and_without.at(i)->SetLogx();
        cans_with_and_without.at(i)->SetLogy();
        cans_with_and_without.at(i)->SetGridx();
        cans_with_and_without.at(i)->SetGridy();
        graphs.at(i)->GetXaxis()->SetTitle("energy [MeV]");
        graphs.at(i)->GetXaxis()->SetTitleOffset(1.1);
        graphs.at(i)->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
        graphs.at(i)->GetYaxis()->SetTitleOffset(1.1);
        graphs.at(i)->Draw("AP");
        leg->AddEntry(graphs.at(i),"no Lpm","P");
        graphs.at(i+9)->Draw("P same");
        leg->AddEntry(graphs.at(i+9),"Lpm","P");
        leg->Draw("same");
        cans_with_and_without.at(i)->Write();

    }

    TCanvas* muons      =   new TCanvas("muons","muons",1024,768);
    TCanvas* taus       =   new TCanvas("taus","taus",1024,768);
    TCanvas* electrons  =   new TCanvas("electrons","electrons",1024,768);

    TCanvas* hydrogen   =   new TCanvas("hydrogen","hydrogen",1024,768);
    TCanvas* water      =   new TCanvas("water","water",1024,768);
    TCanvas* uranium    =   new TCanvas("uranium","uranium",1024,768);

    TMultiGraph* muons_gr      =   new TMultiGraph("muons","muons");
    TMultiGraph* taus_gr       =   new TMultiGraph("taus","taus");
    TMultiGraph* electrons_gr  =   new TMultiGraph("electrons","electrons");

    TMultiGraph* hydrogen_gr   =   new TMultiGraph("hydrogen","hydrogen");
    TMultiGraph* water_gr      =   new TMultiGraph("water","water");
    TMultiGraph* uranium_gr    =   new TMultiGraph("uranium","uranium");

    TLegend* muons_leg = new TLegend(0.7,0.2,0.9,0.4);
    muons_leg->SetFillColor(0);

    TLegend* taus_leg = new TLegend(0.7,0.2,0.9,0.4);
    taus_leg->SetFillColor(0);

    TLegend* electrons_leg = new TLegend(0.7,0.2,0.9,0.4);
    electrons_leg->SetFillColor(0);

    TLegend* hydrogen_leg = new TLegend(0.7,0.2,0.9,0.4);
    hydrogen_leg->SetFillColor(0);

    TLegend* water_leg = new TLegend(0.7,0.2,0.9,0.4);
    water_leg->SetFillColor(0);

    TLegend* uranium_leg = new TLegend(0.7,0.2,0.9,0.4);
    uranium_leg->SetFillColor(0);

    stringstream leg_entry;

    for(unsigned int i = 0 ; i < graphs.size() ; i++)
    {
        if(epair.at(i)->GetParticle()->GetName().compare("mu")==0)
        {
            muons_gr->Add(graphs.at(i),"P");
            leg_entry<<epair.at(i)->GetMedium()->GetName().c_str();
            if(epair.at(i)->GetLpmEffectEnabled())
            {
                leg_entry<<" lpm";
            }

            muons_leg->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(muons_gr->GetListOfGraphs()->Capacity()==1||muons_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(muons_gr->GetListOfGraphs()->Capacity()==2||muons_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(muons_gr->GetListOfGraphs()->Capacity()==3||muons_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }

            leg_entry.str("");
            leg_entry.clear();

        }
        if(epair.at(i)->GetParticle()->GetName().compare("tau")==0)
        {
            taus_gr->Add(graphs.at(i),"P");

            leg_entry<<epair.at(i)->GetMedium()->GetName().c_str();
            if(epair.at(i)->GetLpmEffectEnabled())
            {
                leg_entry<<" lpm";
            }

            taus_leg->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(taus_gr->GetListOfGraphs()->Capacity()==1||taus_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(taus_gr->GetListOfGraphs()->Capacity()==2||taus_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(taus_gr->GetListOfGraphs()->Capacity()==3||taus_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }

            leg_entry.str("");
            leg_entry.clear();

        }
        if(epair.at(i)->GetParticle()->GetName().compare("e")==0)
        {
            electrons_gr->Add(graphs.at(i),"P");

            leg_entry<<epair.at(i)->GetMedium()->GetName().c_str();
            if(epair.at(i)->GetLpmEffectEnabled())
            {
                leg_entry<<" lpm";
            }

            electrons_leg->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(electrons_gr->GetListOfGraphs()->Capacity()==1||electrons_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(electrons_gr->GetListOfGraphs()->Capacity()==2||electrons_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(electrons_gr->GetListOfGraphs()->Capacity()==3||electrons_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }

            leg_entry.str("");
            leg_entry.clear();

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

    for(unsigned int i = 0 ; i < graphs.size() ; i++)
    {
        if(epair.at(i)->GetMedium()->GetName().compare("hydrogen")==0)
        {
            hydrogen_gr->Add(graphs.at(i),"P");

            leg_entry<<epair.at(i)->GetParticle()->GetName().c_str();
            if(epair.at(i)->GetLpmEffectEnabled())
            {
                leg_entry<<" lpm";
            }

            hydrogen_leg->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(hydrogen_gr->GetListOfGraphs()->Capacity()==1||hydrogen_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(hydrogen_gr->GetListOfGraphs()->Capacity()==2||hydrogen_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(hydrogen_gr->GetListOfGraphs()->Capacity()==3||hydrogen_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }

            leg_entry.str("");
            leg_entry.clear();

        }
        if(epair.at(i)->GetMedium()->GetName().compare("water")==0)
        {
            water_gr->Add(graphs.at(i),"P");

            leg_entry<<epair.at(i)->GetParticle()->GetName().c_str();
            if(epair.at(i)->GetLpmEffectEnabled())
            {
                leg_entry<<" lpm";
            }

            water_leg->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(water_gr->GetListOfGraphs()->Capacity()==1||water_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(water_gr->GetListOfGraphs()->Capacity()==2||water_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(water_gr->GetListOfGraphs()->Capacity()==3||water_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }

            leg_entry.str("");
            leg_entry.clear();

        }
        if(epair.at(i)->GetMedium()->GetName().compare("uranium")==0)
        {
            uranium_gr->Add(graphs.at(i),"P");

            leg_entry<<epair.at(i)->GetParticle()->GetName().c_str();
            if(epair.at(i)->GetLpmEffectEnabled())
            {
                leg_entry<<" lpm";
            }

            uranium_leg->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(uranium_gr->GetListOfGraphs()->Capacity()==1||uranium_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(uranium_gr->GetListOfGraphs()->Capacity()==2||uranium_gr->GetListOfGraphs()->Capacity()==5)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(uranium_gr->GetListOfGraphs()->Capacity()==3||uranium_gr->GetListOfGraphs()->Capacity()==6)
            {
                graphs.at(i)->SetMarkerColor(kGreen);

            }

            leg_entry.str("");
            leg_entry.clear();

        }
    }


    hydrogen->cd();
    hydrogen->SetLogx();
    hydrogen->SetLogy();
    hydrogen->SetGridx();
    hydrogen->SetGridy();
    hydrogen_gr->Draw("A");
    hydrogen_gr->GetXaxis()->SetTitle("energy [MeV]");
    hydrogen_gr->GetXaxis()->SetTitleOffset(1.1);
    hydrogen_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    hydrogen_gr->GetYaxis()->SetTitleOffset(1.1);
    hydrogen_leg->Draw("Same");
    hydrogen->Write();

    water->cd();
    water->SetLogx();
    water->SetLogy();
    water->SetGridx();
    water->SetGridy();
    water_gr->Draw("A");
    water_gr->GetXaxis()->SetTitle("energy [MeV]");
    water_gr->GetXaxis()->SetTitleOffset(1.1);
    water_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    water_gr->GetYaxis()->SetTitleOffset(1.1);
    water_leg->Draw("Same");
    water->Write();

    uranium->cd();
    uranium->SetLogx();
    uranium->SetLogy();
    uranium->SetGridx();
    uranium->SetGridy();
    uranium_gr->Draw("A");
    uranium_gr->GetXaxis()->SetTitle("energy [MeV]");
    uranium_gr->GetXaxis()->SetTitleOffset(1.1);
    uranium_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    uranium_gr->GetYaxis()->SetTitleOffset(1.1);
    uranium_leg->Draw("Same");
    uranium->Write();


    file->Close();

    return 0;
}
















