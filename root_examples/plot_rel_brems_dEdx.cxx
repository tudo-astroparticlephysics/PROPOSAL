/** \file
 * This is an example to plot dE/dx per energy for Bremsstrahlung.\n
 * If you have ROOT installed this file will be compilied when running make\n
 * and you can find an executable in your build directory in root_examples:\n
 * -->   Plot_rel_Bremsstrahlung_dEdx   <--\n
 * A ROOT-file named Bremsstrahlung_rel_dEdx.root will be generated,
 * which contains TCanvas's with TGraphs dE/dx per energy vs. energy.\n
 * These TGraphs are generated for muons, taus and electrons in
 * water, hydrogen and uranium. All four parametrizations are shown:\n
 * Kelner-Kakoulin-Petrukhin\n
 * Andreev-Bezrukov-Bugaev\n
 * Petrukhin-Shestakov\n
 * Complete screening case\n
 * You will also find summary plots for every particle and medium.
 * @brief Example to plot dE/dx for Bremsstrahlung
 * @author Jan-Hendrik Köhne
 */


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
    cout<<"\n------------------------------------------------------------------\n"
        <<"This is an example to plot dE/dx for Bremsstrahlung.\n"
        <<"A ROOT-file named Bremsstrahlung_dEdx.root will be generated,\n"
        <<"which contains TCanvas's with TGraphs dE/dx vs. energy.\n"
        <<"These TGraphs are generated for muons, taus and electrons in\n"
        <<"water, hydrogen and uranium. All four parametrizations are shown:\n"
        <<"Kelner-Kakoulin-Petrukhin\n"
        <<"Andreev-Bezrukov-Bugaev\n"
        <<"Petrukhin-Shestakov\n"
        <<"Complete screening case\n"
        <<"You will also find summary plots for every particle and medium.\n"
        <<"------------------------------------------------------------------\n"
        <<endl;


    TFile *file     =   new TFile("Bremsstrahlung_rel_dEdx.root","RECREATE");

    Particle *mu    =   new Particle("mu");
    Particle *tau   =   new Particle("tau");
    Particle *e     =   new Particle("e");

    Medium  *med1   =   new Medium("hydrogen",1.);
    Medium  *med2   =   new Medium("water",1.);
    Medium  *med3   =   new Medium("uranium",1.);

    EnergyCutSettings*  cuts    =   new EnergyCutSettings(-1,-1);

    vector<Bremsstrahlung*> brems;

    brems.push_back(new Bremsstrahlung(mu,med1,cuts));
    brems.push_back(new Bremsstrahlung(mu,med2,cuts));
    brems.push_back(new Bremsstrahlung(mu,med3,cuts));
    brems.push_back(new Bremsstrahlung(tau,med1,cuts));
    brems.push_back(new Bremsstrahlung(tau,med2,cuts));
    brems.push_back(new Bremsstrahlung(tau,med3,cuts));
    brems.push_back(new Bremsstrahlung(e,med1,cuts));
    brems.push_back(new Bremsstrahlung(e,med2,cuts));
    brems.push_back(new Bremsstrahlung(e,med3,cuts));

    brems.push_back(new Bremsstrahlung(mu,med1,cuts));
    brems.push_back(new Bremsstrahlung(mu,med2,cuts));
    brems.push_back(new Bremsstrahlung(mu,med3,cuts));
    brems.push_back(new Bremsstrahlung(tau,med1,cuts));
    brems.push_back(new Bremsstrahlung(tau,med2,cuts));
    brems.push_back(new Bremsstrahlung(tau,med3,cuts));
    brems.push_back(new Bremsstrahlung(e,med1,cuts));
    brems.push_back(new Bremsstrahlung(e,med2,cuts));
    brems.push_back(new Bremsstrahlung(e,med3,cuts));

    brems.push_back(new Bremsstrahlung(mu,med1,cuts));
    brems.push_back(new Bremsstrahlung(mu,med2,cuts));
    brems.push_back(new Bremsstrahlung(mu,med3,cuts));
    brems.push_back(new Bremsstrahlung(tau,med1,cuts));
    brems.push_back(new Bremsstrahlung(tau,med2,cuts));
    brems.push_back(new Bremsstrahlung(tau,med3,cuts));
    brems.push_back(new Bremsstrahlung(e,med1,cuts));
    brems.push_back(new Bremsstrahlung(e,med2,cuts));
    brems.push_back(new Bremsstrahlung(e,med3,cuts));

    brems.push_back(new Bremsstrahlung(mu,med1,cuts));
    brems.push_back(new Bremsstrahlung(mu,med2,cuts));
    brems.push_back(new Bremsstrahlung(mu,med3,cuts));
    brems.push_back(new Bremsstrahlung(tau,med1,cuts));
    brems.push_back(new Bremsstrahlung(tau,med2,cuts));
    brems.push_back(new Bremsstrahlung(tau,med3,cuts));
    brems.push_back(new Bremsstrahlung(e,med1,cuts));
    brems.push_back(new Bremsstrahlung(e,med2,cuts));
    brems.push_back(new Bremsstrahlung(e,med3,cuts));

    for(unsigned int i = 0; i<brems.size() ; i++)
    {
        if(i<9) brems.at(i)->SetParametrization(1);
        else if(i<18) brems.at(i)->SetParametrization(2);
        else if(i<27) brems.at(i)->SetParametrization(3);
        else if(i<36) brems.at(i)->SetParametrization(4);

    }

    vector<TGraph*> graphs;
    graphs.resize(brems.size());


    vector<double> dEdx;
    dEdx.resize(brems.size());

    double energy;

    vector<int> counter(graphs.size(),0);

    stringstream graph_name;
    stringstream graph_title;

    vector<TCanvas*> cans;
    cans.resize(graphs.size());

    for(unsigned int i = 0 ; i < brems.size() ; i++)
    {
        brems.at(i)->EnableDEdxInterpolation();
        graphs.at(i)    =   new TGraph();

        graph_name<<brems.at(i)->GetParticle()->GetName()<<"_"<<brems.at(i)->GetMedium()->GetName();
        graph_name<<"_"<<brems.at(i)->GetParametrization();
        graph_title<<brems.at(i)->GetParticle()->GetName()<<" in "<<brems.at(i)->GetMedium()->GetName();

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

            dEdx.at(i)  =   brems.at(i)->CalculatedEdx();

            if(dEdx.at(i) !=0 )
            {
                graphs.at(i)->SetPoint(counter.at(i),energy,dEdx.at(i)/energy);
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
        graphs.at(i)->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
        graphs.at(i)->Draw("AP");
        cans.at(i)->Write();

    }

    TCanvas* muons_water      =   new TCanvas("muons_water","muons in water",1024,768);
    TCanvas* taus_water       =   new TCanvas("taus_water","taus in water",1024,768);
    TCanvas* electrons_water  =   new TCanvas("electrons_water","electrons in water",1024,768);

    TCanvas* muons_hydrogen      =   new TCanvas("muons_hydrogen","muons in hydrogen",1024,768);
    TCanvas* taus_hydrogen       =   new TCanvas("taus_hydrogen","taus in hydrogen",1024,768);
    TCanvas* electrons_hydrogen  =   new TCanvas("electrons_hydrogen","electrons in hydrogen",1024,768);

    TCanvas* muons_uranium      =   new TCanvas("muons_uranium","muons in uranium",1024,768);
    TCanvas* taus_uranium       =   new TCanvas("taus_uranium","taus in uranium",1024,768);
    TCanvas* electrons_uranium  =   new TCanvas("electrons_uranium","electrons in uranium",1024,768);

    TMultiGraph* muons_water_gr      =   new TMultiGraph("muons_water","muons in water");
    TMultiGraph* taus_water_gr       =   new TMultiGraph("taus_water","taus in water");
    TMultiGraph* electrons_water_gr  =   new TMultiGraph("electrons_water","electrons in water");

    TMultiGraph* muons_hydrogen_gr      =   new TMultiGraph("muons_hydrogen","muons in hydrogen");
    TMultiGraph* taus_hydrogen_gr       =   new TMultiGraph("taus_hydrogen","taus in hydrogen");
    TMultiGraph* electrons_hydrogen_gr  =   new TMultiGraph("electrons_hydrogen","electrons in hydrogen");

    TMultiGraph* muons_uranium_gr      =   new TMultiGraph("muons_uranium","muons in uranium");
    TMultiGraph* taus_uranium_gr       =   new TMultiGraph("taus_uranium","taus in uranium");
    TMultiGraph* electrons_uranium_gr  =   new TMultiGraph("electrons_uranium","electrons in uranium");


    TCanvas* hydrogen_para1   =   new TCanvas("hydrogen_para1","hydrogen parametrization KKP",1024,768);
    TCanvas* water_para1      =   new TCanvas("water_para1","water parametrization KKP",1024,768);
    TCanvas* uranium_para1    =   new TCanvas("uranium_para1","uranium parametrization KKP",1024,768);

    TCanvas* hydrogen_para2   =   new TCanvas("hydrogen_para2","hydrogen parametrization ABB",1024,768);
    TCanvas* water_para2      =   new TCanvas("water_para2","water parametrization ABB",1024,768);
    TCanvas* uranium_para2    =   new TCanvas("uranium_para2","uranium parametrization ABB",1024,768);

    TCanvas* hydrogen_para3   =   new TCanvas("hydrogen_para3","hydrogen parametrization PS",1024,768);
    TCanvas* water_para3      =   new TCanvas("water_para3","water parametrization PS",1024,768);
    TCanvas* uranium_para3    =   new TCanvas("uranium_para3","uranium parametrization PS",1024,768);

    TCanvas* hydrogen_para4   =   new TCanvas("hydrogen_para4","hydrogen parametrization CSC",1024,768);
    TCanvas* water_para4      =   new TCanvas("water_para4","water parametrization CSC",1024,768);
    TCanvas* uranium_para4    =   new TCanvas("uranium_para4","uranium parametrization CSC",1024,768);

    TMultiGraph* hydrogen_para1_gr   =   new TMultiGraph("hydrogen_para1","hydrogen, parametrization = KKP");
    TMultiGraph* water_para1_gr      =   new TMultiGraph("water_para1","water, parametrization = KKP");
    TMultiGraph* uranium_para1_gr    =   new TMultiGraph("uranium_para1","uranium, parametrization = KKP");

    TMultiGraph* hydrogen_para2_gr   =   new TMultiGraph("hydrogen_para2","hydrogen, parametrization = ABB");
    TMultiGraph* water_para2_gr      =   new TMultiGraph("water_para2","water, parametrization = ABB");
    TMultiGraph* uranium_para2_gr    =   new TMultiGraph("uranium_para2","uranium, parametrization = ABB");

    TMultiGraph* hydrogen_para3_gr   =   new TMultiGraph("hydrogen_para3","hydrogen, parametrization = PS");
    TMultiGraph* water_para3_gr      =   new TMultiGraph("water_para3","water, parametrization = PS");
    TMultiGraph* uranium_para3_gr    =   new TMultiGraph("uranium_para3","uranium, parametrization = PS");

    TMultiGraph* hydrogen_para4_gr   =   new TMultiGraph("hydrogen_para4","hydrogen, parametrization = CSC");
    TMultiGraph* water_para4_gr      =   new TMultiGraph("water_para4","water, parametrization = CSC");
    TMultiGraph* uranium_para4_gr    =   new TMultiGraph("uranium_para4","uranium, parametrization = CSC");


    TLegend* muons_leg_water = new TLegend(0.7,0.2,0.9,0.4);
    muons_leg_water->SetFillColor(0);

    TLegend* muons_leg_uranium = new TLegend(0.7,0.2,0.9,0.4);
    muons_leg_uranium->SetFillColor(0);

    TLegend* muons_leg_hydrogen = new TLegend(0.7,0.2,0.9,0.4);
    muons_leg_hydrogen->SetFillColor(0);


    TLegend* taus_leg_water = new TLegend(0.7,0.2,0.9,0.4);
    taus_leg_water->SetFillColor(0);

    TLegend* taus_leg_uranium = new TLegend(0.7,0.2,0.9,0.4);
    taus_leg_uranium->SetFillColor(0);

    TLegend* taus_leg_hydrogen = new TLegend(0.7,0.2,0.9,0.4);
    taus_leg_hydrogen->SetFillColor(0);


    TLegend* electrons_leg_water = new TLegend(0.7,0.2,0.9,0.4);
    electrons_leg_water->SetFillColor(0);

    TLegend* electrons_leg_uranium = new TLegend(0.7,0.2,0.9,0.4);
    electrons_leg_uranium->SetFillColor(0);

    TLegend* electrons_leg_hydrogen = new TLegend(0.7,0.2,0.9,0.4);
    electrons_leg_hydrogen->SetFillColor(0);


    TLegend* hydrogen_para1_leg = new TLegend(0.7,0.2,0.9,0.4);
    hydrogen_para1_leg->SetFillColor(0);

    TLegend* hydrogen_para2_leg = new TLegend(0.7,0.2,0.9,0.4);
    hydrogen_para2_leg->SetFillColor(0);

    TLegend* hydrogen_para3_leg = new TLegend(0.7,0.2,0.9,0.4);
    hydrogen_para3_leg->SetFillColor(0);

    TLegend* hydrogen_para4_leg = new TLegend(0.7,0.2,0.9,0.4);
    hydrogen_para4_leg->SetFillColor(0);


    TLegend* water_para1_leg = new TLegend(0.7,0.2,0.9,0.4);
    water_para1_leg->SetFillColor(0);

    TLegend* water_para2_leg = new TLegend(0.7,0.2,0.9,0.4);
    water_para2_leg->SetFillColor(0);

    TLegend* water_para3_leg = new TLegend(0.7,0.2,0.9,0.4);
    water_para3_leg->SetFillColor(0);

    TLegend* water_para4_leg = new TLegend(0.7,0.2,0.9,0.4);
    water_para4_leg->SetFillColor(0);


    TLegend* uranium_para1_leg = new TLegend(0.7,0.2,0.9,0.4);
    uranium_para1_leg->SetFillColor(0);

    TLegend* uranium_para2_leg = new TLegend(0.7,0.2,0.9,0.4);
    uranium_para2_leg->SetFillColor(0);

    TLegend* uranium_para3_leg = new TLegend(0.7,0.2,0.9,0.4);
    uranium_para3_leg->SetFillColor(0);

    TLegend* uranium_para4_leg = new TLegend(0.7,0.2,0.9,0.4);
    uranium_para4_leg->SetFillColor(0);

    stringstream leg_entry;

    for(unsigned int i = 0 ; i < graphs.size() ; i++)
    {
        if(brems.at(i)->GetParticle()->GetName().compare("mu")==0 && brems.at(i)->GetMedium()->GetName().compare("water")==0)
        {
            muons_water_gr->Add(graphs.at(i),"P");

            if(brems.at(i)->GetParametrization()==1)
            {
                leg_entry<<"KKP";
            }
            if(brems.at(i)->GetParametrization()==2)
            {
                leg_entry<<"ABB";
            }
            if(brems.at(i)->GetParametrization()==3)
            {
                leg_entry<<"PS";
            }
            if(brems.at(i)->GetParametrization()==4)
            {
                leg_entry<<"CSC";
            }

            muons_leg_water->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(muons_water_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(muons_water_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(muons_water_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(muons_water_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }

            leg_entry.str("");
            leg_entry.clear();

        }
        if(brems.at(i)->GetParticle()->GetName().compare("mu")==0 && brems.at(i)->GetMedium()->GetName().compare("hydrogen")==0)
        {
            muons_hydrogen_gr->Add(graphs.at(i),"P");

            if(brems.at(i)->GetParametrization()==1)
            {
                leg_entry<<"KKP";
            }
            if(brems.at(i)->GetParametrization()==2)
            {
                leg_entry<<"ABB";
            }
            if(brems.at(i)->GetParametrization()==3)
            {
                leg_entry<<"PS";
            }
            if(brems.at(i)->GetParametrization()==4)
            {
                leg_entry<<"CSC";
            }

            muons_leg_hydrogen->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(muons_hydrogen_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(muons_hydrogen_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(muons_hydrogen_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(muons_hydrogen_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }

            leg_entry.str("");
            leg_entry.clear();

        }
        if(brems.at(i)->GetParticle()->GetName().compare("mu")==0 && brems.at(i)->GetMedium()->GetName().compare("uranium")==0)
        {
            muons_uranium_gr->Add(graphs.at(i),"P");

            if(brems.at(i)->GetParametrization()==1)
            {
                leg_entry<<"KKP";
            }
            if(brems.at(i)->GetParametrization()==2)
            {
                leg_entry<<"ABB";
            }
            if(brems.at(i)->GetParametrization()==3)
            {
                leg_entry<<"PS";
            }
            if(brems.at(i)->GetParametrization()==4)
            {
                leg_entry<<"CSC";
            }

            muons_leg_uranium->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(muons_uranium_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(muons_uranium_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(muons_uranium_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(muons_uranium_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }

            leg_entry.str("");
            leg_entry.clear();

        }


        if(brems.at(i)->GetParticle()->GetName().compare("tau")==0 && brems.at(i)->GetMedium()->GetName().compare("water")==0)
        {
            taus_water_gr->Add(graphs.at(i),"P");

            if(brems.at(i)->GetParametrization()==1)
            {
                leg_entry<<"KKP";
            }
            if(brems.at(i)->GetParametrization()==2)
            {
                leg_entry<<"ABB";
            }
            if(brems.at(i)->GetParametrization()==3)
            {
                leg_entry<<"PS";
            }
            if(brems.at(i)->GetParametrization()==4)
            {
                leg_entry<<"CSC";
            }

            taus_leg_water->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(taus_water_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(taus_water_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(taus_water_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(taus_water_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }

            leg_entry.str("");
            leg_entry.clear();

        }
        if(brems.at(i)->GetParticle()->GetName().compare("tau")==0 && brems.at(i)->GetMedium()->GetName().compare("hydrogen")==0)
        {
            taus_hydrogen_gr->Add(graphs.at(i),"P");

            if(brems.at(i)->GetParametrization()==1)
            {
                leg_entry<<"KKP";
            }
            if(brems.at(i)->GetParametrization()==2)
            {
                leg_entry<<"ABB";
            }
            if(brems.at(i)->GetParametrization()==3)
            {
                leg_entry<<"PS";
            }
            if(brems.at(i)->GetParametrization()==4)
            {
                leg_entry<<"CSC";
            }

            taus_leg_hydrogen->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(taus_hydrogen_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(taus_hydrogen_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(taus_hydrogen_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(taus_hydrogen_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }

            leg_entry.str("");
            leg_entry.clear();

        }
        if(brems.at(i)->GetParticle()->GetName().compare("tau")==0 && brems.at(i)->GetMedium()->GetName().compare("uranium")==0)
        {
            taus_uranium_gr->Add(graphs.at(i),"P");

            if(brems.at(i)->GetParametrization()==1)
            {
                leg_entry<<"KKP";
            }
            if(brems.at(i)->GetParametrization()==2)
            {
                leg_entry<<"ABB";
            }
            if(brems.at(i)->GetParametrization()==3)
            {
                leg_entry<<"PS";
            }
            if(brems.at(i)->GetParametrization()==4)
            {
                leg_entry<<"CSC";
            }

            taus_leg_uranium->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(taus_uranium_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(taus_uranium_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(taus_uranium_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(taus_uranium_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }

            leg_entry.str("");
            leg_entry.clear();

        }


        if(brems.at(i)->GetParticle()->GetName().compare("e")==0 && brems.at(i)->GetMedium()->GetName().compare("water")==0)
        {
            electrons_water_gr->Add(graphs.at(i),"P");

            if(brems.at(i)->GetParametrization()==1)
            {
                leg_entry<<"KKP";
            }
            if(brems.at(i)->GetParametrization()==2)
            {
                leg_entry<<"ABB";
            }
            if(brems.at(i)->GetParametrization()==3)
            {
                leg_entry<<"PS";
            }
            if(brems.at(i)->GetParametrization()==4)
            {
                leg_entry<<"CSC";
            }

            electrons_leg_water->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(electrons_water_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(electrons_water_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(electrons_water_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(electrons_water_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }

            leg_entry.str("");
            leg_entry.clear();

        }
        if(brems.at(i)->GetParticle()->GetName().compare("e")==0 && brems.at(i)->GetMedium()->GetName().compare("hydrogen")==0)
        {
            electrons_hydrogen_gr->Add(graphs.at(i),"P");

            if(brems.at(i)->GetParametrization()==1)
            {
                leg_entry<<"KKP";
            }
            if(brems.at(i)->GetParametrization()==2)
            {
                leg_entry<<"ABB";
            }
            if(brems.at(i)->GetParametrization()==3)
            {
                leg_entry<<"PS";
            }
            if(brems.at(i)->GetParametrization()==4)
            {
                leg_entry<<"CSC";
            }

            electrons_leg_hydrogen->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(electrons_hydrogen_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(electrons_hydrogen_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(electrons_hydrogen_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(electrons_hydrogen_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }

            leg_entry.str("");
            leg_entry.clear();

        }
        if(brems.at(i)->GetParticle()->GetName().compare("e")==0 && brems.at(i)->GetMedium()->GetName().compare("uranium")==0)
        {
            electrons_uranium_gr->Add(graphs.at(i),"P");

            if(brems.at(i)->GetParametrization()==1)
            {
                leg_entry<<"KKP";
            }
            if(brems.at(i)->GetParametrization()==2)
            {
                leg_entry<<"ABB";
            }
            if(brems.at(i)->GetParametrization()==3)
            {
                leg_entry<<"PS";
            }
            if(brems.at(i)->GetParametrization()==4)
            {
                leg_entry<<"CSC";
            }

            electrons_leg_uranium->AddEntry(graphs.at(i),leg_entry.str().c_str(),"p");

            if(electrons_uranium_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(electrons_uranium_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(electrons_uranium_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
            if(electrons_uranium_gr->GetListOfGraphs()->Capacity()==4)
            {
                graphs.at(i)->SetMarkerColor(kBlack);
            }

            leg_entry.str("");
            leg_entry.clear();

        }
    }

    muons_water->cd();
    muons_water->SetLogx();
    muons_water->SetLogy();
    muons_water->SetGridx();
    muons_water->SetGridy();
    muons_water_gr->Draw("A");
    muons_water_gr->GetXaxis()->SetTitle("energy [MeV]");
    muons_water_gr->GetXaxis()->SetTitleOffset(1.1);
    muons_water_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    muons_leg_water->Draw("Same");
    muons_water->Write();

    muons_hydrogen->cd();
    muons_hydrogen->SetLogx();
    muons_hydrogen->SetLogy();
    muons_hydrogen->SetGridx();
    muons_hydrogen->SetGridy();
    muons_hydrogen_gr->Draw("A");
    muons_hydrogen_gr->GetXaxis()->SetTitle("energy [MeV]");
    muons_hydrogen_gr->GetXaxis()->SetTitleOffset(1.1);
    muons_hydrogen_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    muons_leg_hydrogen->Draw("Same");
    muons_hydrogen->Write();

    muons_uranium->cd();
    muons_uranium->SetLogx();
    muons_uranium->SetLogy();
    muons_uranium->SetGridx();
    muons_uranium->SetGridy();
    muons_uranium_gr->Draw("A");
    muons_uranium_gr->GetXaxis()->SetTitle("energy [MeV]");
    muons_uranium_gr->GetXaxis()->SetTitleOffset(1.1);
    muons_uranium_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    muons_leg_uranium->Draw("Same");
    muons_uranium->Write();


    taus_water->cd();
    taus_water->SetLogx();
    taus_water->SetLogy();
    taus_water->SetGridx();
    taus_water->SetGridy();
    taus_water_gr->Draw("A");
    taus_water_gr->GetXaxis()->SetTitle("energy [MeV]");
    taus_water_gr->GetXaxis()->SetTitleOffset(1.1);
    taus_water_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    taus_leg_water->Draw("Same");
    taus_water->Write();

    taus_hydrogen->cd();
    taus_hydrogen->SetLogx();
    taus_hydrogen->SetLogy();
    taus_hydrogen->SetGridx();
    taus_hydrogen->SetGridy();
    taus_hydrogen_gr->Draw("A");
    taus_hydrogen_gr->GetXaxis()->SetTitle("energy [MeV]");
    taus_hydrogen_gr->GetXaxis()->SetTitleOffset(1.1);
    taus_hydrogen_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    taus_leg_hydrogen->Draw("Same");
    taus_hydrogen->Write();

    taus_uranium->cd();
    taus_uranium->SetLogx();
    taus_uranium->SetLogy();
    taus_uranium->SetGridx();
    taus_uranium->SetGridy();
    taus_uranium_gr->Draw("A");
    taus_uranium_gr->GetXaxis()->SetTitle("energy [MeV]");
    taus_uranium_gr->GetXaxis()->SetTitleOffset(1.1);
    taus_uranium_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    taus_leg_uranium->Draw("Same");
    taus_uranium->Write();


    electrons_water->cd();
    electrons_water->SetLogx();
    electrons_water->SetLogy();
    electrons_water->SetGridx();
    electrons_water->SetGridy();
    electrons_water_gr->Draw("A");
    electrons_water_gr->GetXaxis()->SetTitle("energy [MeV]");
    electrons_water_gr->GetXaxis()->SetTitleOffset(1.1);
    electrons_water_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    electrons_leg_water->Draw("Same");
    electrons_water->Write();

    electrons_hydrogen->cd();
    electrons_hydrogen->SetLogx();
    electrons_hydrogen->SetLogy();
    electrons_hydrogen->SetGridx();
    electrons_hydrogen->SetGridy();
    electrons_hydrogen_gr->Draw("A");
    electrons_hydrogen_gr->GetXaxis()->SetTitle("energy [MeV]");
    electrons_hydrogen_gr->GetXaxis()->SetTitleOffset(1.1);
    electrons_hydrogen_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    electrons_leg_hydrogen->Draw("Same");
    electrons_hydrogen->Write();

    electrons_uranium->cd();
    electrons_uranium->SetLogx();
    electrons_uranium->SetLogy();
    electrons_uranium->SetGridx();
    electrons_uranium->SetGridy();
    electrons_uranium_gr->Draw("A");
    electrons_uranium_gr->GetXaxis()->SetTitle("energy [MeV]");
    electrons_uranium_gr->GetXaxis()->SetTitleOffset(1.1);
    electrons_uranium_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    electrons_leg_uranium->Draw("Same");
    electrons_uranium->Write();

    for(unsigned int i = 0 ; i < graphs.size() ; i++)
    {
        if(brems.at(i)->GetMedium()->GetName().compare("hydrogen")==0 && brems.at(i)->GetParametrization()==1)
        {
            hydrogen_para1_gr->Add(graphs.at(i),"P");
            hydrogen_para1_leg->AddEntry(graphs.at(i),brems.at(i)->GetParticle()->GetName().c_str(),"p");

            if(hydrogen_para1_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(hydrogen_para1_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(hydrogen_para1_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
        }
        if(brems.at(i)->GetMedium()->GetName().compare("hydrogen")==0 && brems.at(i)->GetParametrization()==2)
        {
            hydrogen_para2_gr->Add(graphs.at(i),"P");
            hydrogen_para2_leg->AddEntry(graphs.at(i),brems.at(i)->GetParticle()->GetName().c_str(),"p");

            if(hydrogen_para2_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(hydrogen_para2_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(hydrogen_para2_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
        }
        if(brems.at(i)->GetMedium()->GetName().compare("hydrogen")==0 && brems.at(i)->GetParametrization()==3)
        {
            hydrogen_para3_gr->Add(graphs.at(i),"P");
            hydrogen_para3_leg->AddEntry(graphs.at(i),brems.at(i)->GetParticle()->GetName().c_str(),"p");

            if(hydrogen_para3_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(hydrogen_para3_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(hydrogen_para3_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
        }
        if(brems.at(i)->GetMedium()->GetName().compare("hydrogen")==0 && brems.at(i)->GetParametrization()==4)
        {
            hydrogen_para4_gr->Add(graphs.at(i),"P");
            hydrogen_para4_leg->AddEntry(graphs.at(i),brems.at(i)->GetParticle()->GetName().c_str(),"p");

            if(hydrogen_para4_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(hydrogen_para4_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(hydrogen_para4_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
        }

        if(brems.at(i)->GetMedium()->GetName().compare("water")==0 && brems.at(i)->GetParametrization()==1)
        {
            water_para1_gr->Add(graphs.at(i),"P");
            water_para1_leg->AddEntry(graphs.at(i),brems.at(i)->GetParticle()->GetName().c_str(),"p");

            if(water_para1_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(water_para1_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(water_para1_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
        }
        if(brems.at(i)->GetMedium()->GetName().compare("water")==0 && brems.at(i)->GetParametrization()==2)
        {
            water_para2_gr->Add(graphs.at(i),"P");
            water_para2_leg->AddEntry(graphs.at(i),brems.at(i)->GetParticle()->GetName().c_str(),"p");

            if(water_para2_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(water_para2_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(water_para2_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
        }
        if(brems.at(i)->GetMedium()->GetName().compare("water")==0 && brems.at(i)->GetParametrization()==3)
        {
            water_para3_gr->Add(graphs.at(i),"P");
            water_para3_leg->AddEntry(graphs.at(i),brems.at(i)->GetParticle()->GetName().c_str(),"p");

            if(water_para3_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(water_para3_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(water_para3_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
        }
        if(brems.at(i)->GetMedium()->GetName().compare("water")==0 && brems.at(i)->GetParametrization()==4)
        {
            water_para4_gr->Add(graphs.at(i),"P");
            water_para4_leg->AddEntry(graphs.at(i),brems.at(i)->GetParticle()->GetName().c_str(),"p");

            if(water_para4_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(water_para4_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(water_para4_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
        }


        if(brems.at(i)->GetMedium()->GetName().compare("uranium")==0 && brems.at(i)->GetParametrization()==1)
        {
            uranium_para1_gr->Add(graphs.at(i),"P");
            uranium_para1_leg->AddEntry(graphs.at(i),brems.at(i)->GetParticle()->GetName().c_str(),"p");

            if(uranium_para1_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(uranium_para1_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(uranium_para1_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
        }
        if(brems.at(i)->GetMedium()->GetName().compare("uranium")==0 && brems.at(i)->GetParametrization()==2)
        {
            uranium_para2_gr->Add(graphs.at(i),"P");
            uranium_para2_leg->AddEntry(graphs.at(i),brems.at(i)->GetParticle()->GetName().c_str(),"p");

            if(uranium_para2_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(uranium_para2_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(uranium_para2_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
        }
        if(brems.at(i)->GetMedium()->GetName().compare("uranium")==0 && brems.at(i)->GetParametrization()==3)
        {
            uranium_para3_gr->Add(graphs.at(i),"P");
            uranium_para3_leg->AddEntry(graphs.at(i),brems.at(i)->GetParticle()->GetName().c_str(),"p");

            if(uranium_para3_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(uranium_para3_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(uranium_para3_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
        }
        if(brems.at(i)->GetMedium()->GetName().compare("uranium")==0 && brems.at(i)->GetParametrization()==4)
        {
            uranium_para4_gr->Add(graphs.at(i),"P");
            uranium_para4_leg->AddEntry(graphs.at(i),brems.at(i)->GetParticle()->GetName().c_str(),"p");

            if(uranium_para4_gr->GetListOfGraphs()->Capacity()==1)
            {
                graphs.at(i)->SetMarkerColor(kRed);
            }
            if(uranium_para4_gr->GetListOfGraphs()->Capacity()==2)
            {
                graphs.at(i)->SetMarkerColor(kBlue);
            }
            if(uranium_para4_gr->GetListOfGraphs()->Capacity()==3)
            {
                graphs.at(i)->SetMarkerColor(kGreen);
            }
        }
    }


    hydrogen_para1->cd();
    hydrogen_para1->SetLogx();
    hydrogen_para1->SetLogy();
    hydrogen_para1->SetGridx();
    hydrogen_para1->SetGridy();
    hydrogen_para1_gr->Draw("A");
    hydrogen_para1_gr->GetXaxis()->SetTitle("energy [MeV]");
    hydrogen_para1_gr->GetXaxis()->SetTitleOffset(1.1);
    hydrogen_para1_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    hydrogen_para1_leg->Draw("Same");
    hydrogen_para1->Write();

    hydrogen_para2->cd();
    hydrogen_para2->SetLogx();
    hydrogen_para2->SetLogy();
    hydrogen_para2->SetGridx();
    hydrogen_para2->SetGridy();
    hydrogen_para2_gr->Draw("A");
    hydrogen_para2_gr->GetXaxis()->SetTitle("energy [MeV]");
    hydrogen_para2_gr->GetXaxis()->SetTitleOffset(1.1);
    hydrogen_para2_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    hydrogen_para2_leg->Draw("Same");
    hydrogen_para2->Write();

    hydrogen_para3->cd();
    hydrogen_para3->SetLogx();
    hydrogen_para3->SetLogy();
    hydrogen_para3->SetGridx();
    hydrogen_para3->SetGridy();
    hydrogen_para3_gr->Draw("A");
    hydrogen_para3_gr->GetXaxis()->SetTitle("energy [MeV]");
    hydrogen_para3_gr->GetXaxis()->SetTitleOffset(1.1);
    hydrogen_para3_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    hydrogen_para3_leg->Draw("Same");
    hydrogen_para3->Write();

    hydrogen_para4->cd();
    hydrogen_para4->SetLogx();
    hydrogen_para4->SetLogy();
    hydrogen_para4->SetGridx();
    hydrogen_para4->SetGridy();
    hydrogen_para4_gr->Draw("A");
    hydrogen_para4_gr->GetXaxis()->SetTitle("energy [MeV]");
    hydrogen_para4_gr->GetXaxis()->SetTitleOffset(1.1);
    hydrogen_para4_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    hydrogen_para4_leg->Draw("Same");
    hydrogen_para4->Write();


    water_para1->cd();
    water_para1->SetLogx();
    water_para1->SetLogy();
    water_para1->SetGridx();
    water_para1->SetGridy();
    water_para1_gr->Draw("A");
    water_para1_gr->GetXaxis()->SetTitle("energy [MeV]");
    water_para1_gr->GetXaxis()->SetTitleOffset(1.1);
    water_para1_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    water_para1_leg->Draw("Same");
    water_para1->Write();

    water_para2->cd();
    water_para2->SetLogx();
    water_para2->SetLogy();
    water_para2->SetGridx();
    water_para2->SetGridy();
    water_para2_gr->Draw("A");
    water_para2_gr->GetXaxis()->SetTitle("energy [MeV]");
    water_para2_gr->GetXaxis()->SetTitleOffset(1.1);
    water_para2_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    water_para2_leg->Draw("Same");
    water_para2->Write();

    water_para3->cd();
    water_para3->SetLogx();
    water_para3->SetLogy();
    water_para3->SetGridx();
    water_para3->SetGridy();
    water_para3_gr->Draw("A");
    water_para3_gr->GetXaxis()->SetTitle("energy [MeV]");
    water_para3_gr->GetXaxis()->SetTitleOffset(1.1);
    water_para3_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    water_para3_leg->Draw("Same");
    water_para3->Write();

    water_para4->cd();
    water_para4->SetLogx();
    water_para4->SetLogy();
    water_para4->SetGridx();
    water_para4->SetGridy();
    water_para4_gr->Draw("A");
    water_para4_gr->GetXaxis()->SetTitle("energy [MeV]");
    water_para4_gr->GetXaxis()->SetTitleOffset(1.1);
    water_para4_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    water_para4_leg->Draw("Same");
    water_para4->Write();


    uranium_para1->cd();
    uranium_para1->SetLogx();
    uranium_para1->SetLogy();
    uranium_para1->SetGridx();
    uranium_para1->SetGridy();
    uranium_para1_gr->Draw("A");
    uranium_para1_gr->GetXaxis()->SetTitle("energy [MeV]");
    uranium_para1_gr->GetXaxis()->SetTitleOffset(1.1);
    uranium_para1_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    uranium_para1_leg->Draw("Same");
    uranium_para1->Write();

    uranium_para2->cd();
    uranium_para2->SetLogx();
    uranium_para2->SetLogy();
    uranium_para2->SetGridx();
    uranium_para2->SetGridy();
    uranium_para2_gr->Draw("A");
    uranium_para2_gr->GetXaxis()->SetTitle("energy [MeV]");
    uranium_para2_gr->GetXaxis()->SetTitleOffset(1.1);
    uranium_para2_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    uranium_para2_leg->Draw("Same");
    uranium_para2->Write();

    uranium_para3->cd();
    uranium_para3->SetLogx();
    uranium_para3->SetLogy();
    uranium_para3->SetGridx();
    uranium_para3->SetGridy();
    uranium_para3_gr->Draw("A");
    uranium_para3_gr->GetXaxis()->SetTitle("energy [MeV]");
    uranium_para3_gr->GetXaxis()->SetTitleOffset(1.1);
    uranium_para3_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    uranium_para3_leg->Draw("Same");
    uranium_para3->Write();

    uranium_para4->cd();
    uranium_para4->SetLogx();
    uranium_para4->SetLogy();
    uranium_para4->SetGridx();
    uranium_para4->SetGridy();
    uranium_para4_gr->Draw("A");
    uranium_para4_gr->GetXaxis()->SetTitle("energy [MeV]");
    uranium_para4_gr->GetXaxis()->SetTitleOffset(1.1);
    uranium_para4_gr->GetYaxis()->SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
    uranium_para4_leg->Draw("Same");
    uranium_para4->Write();

    file->Close();

    return 0;
}
















