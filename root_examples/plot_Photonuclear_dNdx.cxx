/** \file
 * This is an example to plot dN/dx for Photonuclear.\n
 * If you have ROOT installed this file will be compilied when running make\n
 * and you can find an executable in your build directory in root_examples:\n
 * -->   Plot_Photonuclear_dNdx   <--\n
 * A ROOT-file named Photonuclear_dNdx.root will be generated,
 * which contains TCanvas's with TGraphs dN/dx vs. energy.\n
 * These TGraphs are generated for muons, taus and electrons in
 * water, hydrogen and uranium. All 14 parametrizations are shown:\n
 * Kokoulin\n
 * Kokoulin + hard component\n
 * Rhode\n
 * Rhode + hard component\n
 * Bezrukov/Bugaev\n
 * Bezrukov/Bugaev + hard component\n
 * Zeus\n
 * Zeus + hard component\n
 * ALLM 91 shadow=Dutta\n
 * ALLM 91 shadow=Butkevich\n
 * ALLM 97 shadow=Dutta\n
 * ALLM 97 shadow=Butkevich\n
 * Butkevich/Mikhailov shadow=Dutta\n
 * Butkevich/Mikhailov shadow=Butkevich\n
 * You will also find summary plots for every particle and medium.
 * @brief Example to plot dN/dx for Photonuclear
 * @author Jan-Hendrik Köhne
 */

#include "TGraph.h"
#include "TCanvas.h"
#include "PROPOSAL/Photonuclear.h"
#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/Output.h"
#include "TFile.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include <cmath>
#include <sstream>

using namespace std;
using namespace PROPOSAL;

int main()
{

    cout<<"\n-----------------------------------------------------------------\n"
        <<"This is an example to plot dN/dx for Photonuclear.\n"
        <<"A ROOT-file named Photonuclear_dNdx.root will be generated,\n"
        <<"which contains TCanvas's with TGraphs dN/dx vs. energy.\n"
        <<"These TGraphs are generated for muons, taus and electrons in\n"
        <<"water, hydrogen and uranium. All 14 parametrizations are shown:\n"
        <<"Kokoulin\n"
        <<"Kokoulin + hard component\n"
        <<"Rhode\n"
        <<"Rhode + hard component\n"
        <<"Bezrukov/Bugaev\n"
        <<"Bezrukov/Bugaev + hard component\n"
        <<"Zeus\n"
        <<"Zeus + hard component\n"
        <<"ALLM 91 shadow=Dutta\n"
        <<"ALLM 91 shadow=Butkevich\n"
        <<"ALLM 97 shadow=Dutta\n"
        <<"ALLM 97 shadow=Butkevich\n"
        <<"Butkevich/Mikhailov shadow=Dutta\n"
        <<"Butkevich/Mikhailov shadow=Butkevich\n"
        <<"You will also find summary plots for every particle and medium.\n"
        <<"-----------------------------------------------------------------\n"
        <<endl;

    TFile *file     =   new TFile("Photonuclear_dNdx.root","RECREATE");

    PROPOSALParticle *mu    =   new PROPOSALParticle(ParticleType::MuMinus);
    PROPOSALParticle *tau   =   new PROPOSALParticle(ParticleType::TauMinus);
    PROPOSALParticle *e     =   new PROPOSALParticle(ParticleType::EMinus);

    Medium  *med1   =   new Medium("hydrogen",1.);
    Medium  *med2   =   new Medium("water",1.);
    Medium  *med3   =   new Medium("uranium",1.);

    EnergyCutSettings*  cuts    =   new EnergyCutSettings(0.1,-1);

    vector<Photonuclear*> photo;
    for(int i =0 ; i<14; i++)
    {
        photo.push_back(new Photonuclear(mu,med1,cuts));
        photo.push_back(new Photonuclear(mu,med2,cuts));
        photo.push_back(new Photonuclear(mu,med3,cuts));
        photo.push_back(new Photonuclear(tau,med1,cuts));
        photo.push_back(new Photonuclear(tau,med2,cuts));
        photo.push_back(new Photonuclear(tau,med3,cuts));
        photo.push_back(new Photonuclear(e,med1,cuts));
        photo.push_back(new Photonuclear(e,med2,cuts));
        photo.push_back(new Photonuclear(e,med3,cuts));
    }

    for(unsigned int i = 0; i<photo.size() ; i++)
    {
        if(i<9) photo.at(i)->SetParametrization(ParametrizationType::PhotoKokoulinShadowBezrukovSoft);
        else if(i<18) photo.at(i)->SetParametrization(ParametrizationType::PhotoKokoulinShadowBezrukovHard);
        else if(i<27) photo.at(i)->SetParametrization(ParametrizationType::PhotoRhodeShadowBezrukovSoft);
        else if(i<36) photo.at(i)->SetParametrization(ParametrizationType::PhotoRhodeShadowBezrukovHard);
        else if(i<45) photo.at(i)->SetParametrization(ParametrizationType::PhotoBezrukovBugaevShadowBezrukovSoft);
        else if(i<54) photo.at(i)->SetParametrization(ParametrizationType::PhotoBezrukovBugaevShadowBezrukovHard);
        else if(i<63) photo.at(i)->SetParametrization(ParametrizationType::PhotoZeusShadowBezrukovSoft);
        else if(i<72) photo.at(i)->SetParametrization(ParametrizationType::PhotoZeusShadowBezrukovHard);
        else if(i<81) photo.at(i)->SetParametrization(ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowDutta);
        else if(i<90) photo.at(i)->SetParametrization(ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowButkevich);
        else if(i<99) photo.at(i)->SetParametrization(ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowDutta);
        else if(i<108) photo.at(i)->SetParametrization(ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich);
        else if(i<117) photo.at(i)->SetParametrization(ParametrizationType::PhotoButkevichMikhailovShadowDutta);
        else if(i<126) photo.at(i)->SetParametrization(ParametrizationType::PhotoButkevichMikhailovShadowButkevich);
    }

    vector<TGraph*> graphs;
    graphs.resize(photo.size());


    vector<double> dNdx;
    dNdx.resize(photo.size());

    double energy;

    vector<int> counter(graphs.size(),0);

    stringstream graph_name;
    stringstream graph_title;

    vector<TCanvas*> cans;
    cans.resize(graphs.size());

    for(unsigned int i = 0 ; i < photo.size() ; i++)
    {
        photo.at(i)->EnableDNdxInterpolation();
        graphs.at(i)    =   new TGraph();

        graph_name<<photo.at(i)->GetParticle()->GetName()<<"_"<<photo.at(i)->GetMedium()->GetName();
        graph_name<<"_"<<photo.at(i)->GetParametrization();
        graph_title<<photo.at(i)->GetParticle()->GetName()<<" in "<<photo.at(i)->GetMedium()->GetName();


        switch (photo.at(i)->GetParametrization())
        {
            case ParametrizationType::PhotoKokoulinShadowBezrukovSoft:
                graph_title<<", parametrization: Kokoulin";
                graphs.at(i)->SetMarkerStyle(2);
                break;
            case ParametrizationType::PhotoKokoulinShadowBezrukovHard:
                graph_title<<", parametrization: Kokoulin + hard component";
                graphs.at(i)->SetMarkerStyle(30);
                break;
            case ParametrizationType::PhotoRhodeShadowBezrukovSoft:
                graph_title<<", parametrization: Rhode";
                graphs.at(i)->SetMarkerStyle(4);
                break;
            case ParametrizationType::PhotoRhodeShadowBezrukovHard:
                graph_title<<", parametrization: Rhode + hard component";
                graphs.at(i)->SetMarkerStyle(5);
                break;
            case ParametrizationType::PhotoBezrukovBugaevShadowBezrukovSoft:
                graph_title<<", parametrization: Bezrukov/Bugaev";
                graphs.at(i)->SetMarkerStyle(20);
                break;
            case ParametrizationType::PhotoBezrukovBugaevShadowBezrukovHard:
                graph_title<<", parametrization: Bezrukov/Bugaev + hard component";
                graphs.at(i)->SetMarkerStyle(21);
                break;
            case ParametrizationType::PhotoZeusShadowBezrukovSoft:
                graph_title<<", parametrization: Zeus";
                graphs.at(i)->SetMarkerStyle(25);
                break;
            case ParametrizationType::PhotoZeusShadowBezrukovHard:
                graph_title<<", parametrization: Zeus + hard component";
                graphs.at(i)->SetMarkerStyle(26);
                break;
            case ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowDutta:
                graph_title<<", parametrization: ALLM 91 shadow=Dutta";
                graphs.at(i)->SetMarkerStyle(22);
                break;
            case ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowButkevich:
                graph_title<<", parametrization: ALLM 91 shadow=Butkevich";
                graphs.at(i)->SetMarkerStyle(27);
                break;
            case ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowDutta:
                graph_title<<", parametrization: ALLM 97 shadow=Dutta";
                graphs.at(i)->SetMarkerStyle(29);
                break;
            case ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich:
                graph_title<<", parametrization: ALLM 97 shadow=Butkevich";
                graphs.at(i)->SetMarkerStyle(3);
                break;
            case ParametrizationType::PhotoButkevichMikhailovShadowDutta:
                graph_title<<", parametrization: Butkevich/Mikhailov shadow=Dutta";
                graphs.at(i)->SetMarkerStyle(23);
                break;
            case ParametrizationType::PhotoButkevichMikhailovShadowButkevich:
                graph_title<<", parametrization: Butkevich/Mikhailov shadow=Butkevich";
                graphs.at(i)->SetMarkerStyle(28);
                break;
            default:
                log_fatal("Wrong Nuclear Interaction Parametrization Type '%i'. The default type is '%i'"
                    , photo.at(i)->GetParametrization(), ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich);
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
        for(unsigned int i = 0 ; i < photo.size() ; i++)
        {
            energy  =   photo.at(i)->GetParticle()->GetMass()   +   pow(10,log_energy);

            photo.at(i)->GetParticle()->SetEnergy(energy);

            dNdx.at(i)  =   photo.at(i)->CalculatedNdx();

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


    TCanvas* hydrogen_para1   =   new TCanvas("hydrogen_para1","hydrogen parametrization Kokoulin",1024,768);
    TCanvas* water_para1      =   new TCanvas("water_para1","water parametrization Kokoulin",1024,768);
    TCanvas* uranium_para1    =   new TCanvas("uranium_para1","uranium parametrization Kokoulin",1024,768);

    TCanvas* hydrogen_para2   =   new TCanvas("hydrogen_para2","hydrogen parametrization Kokoulin + hard",1024,768);
    TCanvas* water_para2      =   new TCanvas("water_para2","water parametrization Kokoulin + hard",1024,768);
    TCanvas* uranium_para2    =   new TCanvas("uranium_para2","uranium parametrization Kokoulin + hard",1024,768);

    TCanvas* hydrogen_para3   =   new TCanvas("hydrogen_para3","hydrogen parametrization Rhode",1024,768);
    TCanvas* water_para3      =   new TCanvas("water_para3","water parametrization Rhode",1024,768);
    TCanvas* uranium_para3    =   new TCanvas("uranium_para3","uranium parametrization Rhode",1024,768);

    TCanvas* hydrogen_para4   =   new TCanvas("hydrogen_para4","hydrogen parametrization Rhode + hard",1024,768);
    TCanvas* water_para4      =   new TCanvas("water_para4","water parametrization Rhode + hard",1024,768);
    TCanvas* uranium_para4    =   new TCanvas("uranium_para4","uranium parametrization Rhode + hard",1024,768);

    TCanvas* hydrogen_para5   =   new TCanvas("hydrogen_para5","hydrogen parametrization Bezrukov/Bugaev",1024,768);
    TCanvas* water_para5      =   new TCanvas("water_para5","water parametrization Bezrukov/Bugaev",1024,768);
    TCanvas* uranium_para5    =   new TCanvas("uranium_para5","uranium parametrization Bezrukov/Bugaev",1024,768);

    TCanvas* hydrogen_para6   =   new TCanvas("hydrogen_para6","hydrogen parametrization Bezrukov/Bugaev + hard",1024,768);
    TCanvas* water_para6      =   new TCanvas("water_para6","water parametrization Bezrukov/Bugaev + hard",1024,768);
    TCanvas* uranium_para6    =   new TCanvas("uranium_para6","uranium parametrization Bezrukov/Bugaev + hard",1024,768);

    TCanvas* hydrogen_para7   =   new TCanvas("hydrogen_para7","hydrogen parametrization Zeus",1024,768);
    TCanvas* water_para7      =   new TCanvas("water_para7","water parametrization Zeus",1024,768);
    TCanvas* uranium_para7    =   new TCanvas("uranium_para7","uranium parametrization Zeus",1024,768);

    TCanvas* hydrogen_para8   =   new TCanvas("hydrogen_para8","hydrogen parametrization Zeus + hard",1024,768);
    TCanvas* water_para8      =   new TCanvas("water_para8","water parametrization Zeus + hard",1024,768);
    TCanvas* uranium_para8    =   new TCanvas("uranium_para8","uranium parametrization Zeus + hard",1024,768);

    TCanvas* hydrogen_para9   =   new TCanvas("hydrogen_para9","hydrogen parametrization ALLM 91 shadow=Dutta",1024,768);
    TCanvas* water_para9      =   new TCanvas("water_para9","water parametrization ALLM 91 shadow=Dutta",1024,768);
    TCanvas* uranium_para9    =   new TCanvas("uranium_para9","uranium parametrization ALLM 91 shadow=Dutta",1024,768);

    TCanvas* hydrogen_para10   =   new TCanvas("hydrogen_para10","hydrogen parametrization ALLM 91 shadow=Butkevich",1024,768);
    TCanvas* water_para10      =   new TCanvas("water_para10","water parametrization ALLM 91 shadow=Butkevich",1024,768);
    TCanvas* uranium_para10    =   new TCanvas("uranium_para10","uranium parametrization ALLM 91 shadow=Butkevich",1024,768);

    TCanvas* hydrogen_para11   =   new TCanvas("hydrogen_para11","hydrogen parametrization ALLM 97 shadow=Dutta",1024,768);
    TCanvas* water_para11      =   new TCanvas("water_para11","water parametrization ALLM 97 shadow=Dutta",1024,768);
    TCanvas* uranium_para11    =   new TCanvas("uranium_para11","uranium parametrization ALLM 97 shadow=Dutta",1024,768);

    TCanvas* hydrogen_para12   =   new TCanvas("hydrogen_para12","hydrogen parametrization ALLM 97 shadow=Butkevich",1024,768);
    TCanvas* water_para12      =   new TCanvas("water_para12","water parametrization ALLM 97 shadow=Butkevich",1024,768);
    TCanvas* uranium_para12    =   new TCanvas("uranium_para12","uranium parametrization ALLM 97 shadow=Butkevich",1024,768);

    TCanvas* hydrogen_para13   =   new TCanvas("hydrogen_para13","hydrogen parametrization Butkevich/Mikhailov shadow=Dutta",1024,768);
    TCanvas* water_para13      =   new TCanvas("water_para13","water parametrization Butkevich/Mikhailov shadow=Dutta",1024,768);
    TCanvas* uranium_para13    =   new TCanvas("uranium_para13","uranium parametrization Butkevich/Mikhailov shadow=Dutta",1024,768);

    TCanvas* hydrogen_para14   =   new TCanvas("hydrogen_para14","hydrogen parametrization Butkevich/Mikhailov shadow=Butkevich",1024,768);
    TCanvas* water_para14      =   new TCanvas("water_para14","water parametrization Butkevich/Mikhailov shadow=Butkevich",1024,768);
    TCanvas* uranium_para14    =   new TCanvas("uranium_para14","uranium parametrization Butkevich/Mikhailov shadow=Butkevich",1024,768);


    TMultiGraph* hydrogen_para1_gr   =   new TMultiGraph("hydrogen_para1","hydrogen, parametrization = Kokoulin");
    TMultiGraph* water_para1_gr      =   new TMultiGraph("water_para1","water, parametrization = Kokoulin");
    TMultiGraph* uranium_para1_gr    =   new TMultiGraph("uranium_para1","uranium, parametrization = Kokoulin");

    TMultiGraph* hydrogen_para2_gr   =   new TMultiGraph("hydrogen_para2","hydrogen, parametrization = Kokoulin + hard");
    TMultiGraph* water_para2_gr      =   new TMultiGraph("water_para2","water, parametrization = Kokoulin + hard");
    TMultiGraph* uranium_para2_gr    =   new TMultiGraph("uranium_para2","uranium, parametrization = Kokoulin + hard");

    TMultiGraph* hydrogen_para3_gr   =   new TMultiGraph("hydrogen_para3","hydrogen, parametrization = Rhode");
    TMultiGraph* water_para3_gr      =   new TMultiGraph("water_para3","water, parametrization = Rhode");
    TMultiGraph* uranium_para3_gr    =   new TMultiGraph("uranium_para3","uranium, parametrization = Rhode");

    TMultiGraph* hydrogen_para4_gr   =   new TMultiGraph("hydrogen_para4","hydrogen, parametrization = Rhode + hard");
    TMultiGraph* water_para4_gr      =   new TMultiGraph("water_para4","water, parametrization = Rhode + hard");
    TMultiGraph* uranium_para4_gr    =   new TMultiGraph("uranium_para4","uranium, parametrization = Rhode + hard");

    TMultiGraph* hydrogen_para5_gr   =   new TMultiGraph("hydrogen_para5","hydrogen, parametrization = Bezrukov/Bugaev");
    TMultiGraph* water_para5_gr      =   new TMultiGraph("water_para5","water, parametrization = Bezrukov/Bugaev");
    TMultiGraph* uranium_para5_gr    =   new TMultiGraph("uranium_para5","uranium, parametrization = Bezrukov/Bugaev");

    TMultiGraph* hydrogen_para6_gr   =   new TMultiGraph("hydrogen_para6","hydrogen, parametrization = Bezrukov/Bugaev + hard");
    TMultiGraph* water_para6_gr      =   new TMultiGraph("water_para6","water, parametrization = Bezrukov/Bugaev + hard");
    TMultiGraph* uranium_para6_gr    =   new TMultiGraph("uranium_para6","uranium, parametrization = Bezrukov/Bugaev + hard");

    TMultiGraph* hydrogen_para7_gr   =   new TMultiGraph("hydrogen_para7","hydrogen, parametrization = Zeus");
    TMultiGraph* water_para7_gr      =   new TMultiGraph("water_para7","water, parametrization = Zeus");
    TMultiGraph* uranium_para7_gr    =   new TMultiGraph("uranium_para7","uranium, parametrization = Zeus");

    TMultiGraph* hydrogen_para8_gr   =   new TMultiGraph("hydrogen_para8","hydrogen, parametrization = Zeus + hard");
    TMultiGraph* water_para8_gr      =   new TMultiGraph("water_para8","water, parametrization = Zeus + hard");
    TMultiGraph* uranium_para8_gr    =   new TMultiGraph("uranium_para8","uranium, parametrization = Zeus + hard");

    TMultiGraph* hydrogen_para9_gr   =   new TMultiGraph("hydrogen_para9","hydrogen, parametrization = ALLM 91 shadow=Dutta");
    TMultiGraph* water_para9_gr      =   new TMultiGraph("water_para9","water, parametrization = ALLM 91 shadow=Dutta");
    TMultiGraph* uranium_para9_gr    =   new TMultiGraph("uranium_para9","uranium, parametrization = ALLM 91 shadow=Dutta");

    TMultiGraph* hydrogen_para10_gr   =   new TMultiGraph("hydrogen_para10","hydrogen, parametrization = ALLM 91 shadow=Butkevich");
    TMultiGraph* water_para10_gr      =   new TMultiGraph("water_para10","water, parametrization = ALLM 91 shadow=Butkevich");
    TMultiGraph* uranium_para10_gr    =   new TMultiGraph("uranium_para10","uranium, parametrization = ALLM 91 shadow=Butkevich");

    TMultiGraph* hydrogen_para11_gr   =   new TMultiGraph("hydrogen_para11","hydrogen, parametrization = ALLM 97 shadow=Dutta");
    TMultiGraph* water_para11_gr      =   new TMultiGraph("water_para11","water, parametrization = ALLM 97 shadow=Dutta");
    TMultiGraph* uranium_para11_gr    =   new TMultiGraph("uranium_para11","uranium, parametrization = ALLM 97 shadow=Dutta");

    TMultiGraph* hydrogen_para12_gr   =   new TMultiGraph("hydrogen_para12","hydrogen, parametrization = ALLM 97 shadow=Butkevich");
    TMultiGraph* water_para12_gr      =   new TMultiGraph("water_para12","water, parametrization = ALLM 97 shadow=Butkevich");
    TMultiGraph* uranium_para12_gr    =   new TMultiGraph("uranium_para12","uranium, parametrization = ALLM 97 shadow=Butkevich");

    TMultiGraph* hydrogen_para13_gr   =   new TMultiGraph("hydrogen_para13","hydrogen, parametrization = Butkevich/Mikhailov shadow=Dutta");
    TMultiGraph* water_para13_gr      =   new TMultiGraph("water_para13","water, parametrization = Butkevich/Mikhailov shadow=Dutta");
    TMultiGraph* uranium_para13_gr    =   new TMultiGraph("uranium_para13","uranium, parametrization = Butkevich/Mikhailov shadow=Dutta");

    TMultiGraph* hydrogen_para14_gr   =   new TMultiGraph("hydrogen_para14","hydrogen, parametrization = Butkevich/Mikhailov shadow=Butkevich");
    TMultiGraph* water_para14_gr      =   new TMultiGraph("water_para14","water, parametrization = Butkevich/Mikhailov shadow=Butkevich");
    TMultiGraph* uranium_para14_gr    =   new TMultiGraph("uranium_para14","uranium, parametrization = Butkevich/Mikhailov shadow=Butkevich");



    TLegend* muons_leg_water = new TLegend(0.6,0.1,0.9,0.4);
    muons_leg_water->SetFillColor(0);

    TLegend* muons_leg_uranium = new TLegend(0.6,0.1,0.9,0.4);
    muons_leg_uranium->SetFillColor(0);

    TLegend* muons_leg_hydrogen = new TLegend(0.6,0.1,0.9,0.4);
    muons_leg_hydrogen->SetFillColor(0);


    TLegend* taus_leg_water = new TLegend(0.6,0.1,0.9,0.4);
    taus_leg_water->SetFillColor(0);

    TLegend* taus_leg_uranium = new TLegend(0.6,0.1,0.9,0.4);
    taus_leg_uranium->SetFillColor(0);

    TLegend* taus_leg_hydrogen = new TLegend(0.6,0.1,0.9,0.4);
    taus_leg_hydrogen->SetFillColor(0);


    TLegend* electrons_leg_water = new TLegend(0.6,0.1,0.9,0.4);
    electrons_leg_water->SetFillColor(0);

    TLegend* electrons_leg_uranium = new TLegend(0.6,0.1,0.9,0.4);
    electrons_leg_uranium->SetFillColor(0);

    TLegend* electrons_leg_hydrogen = new TLegend(0.6,0.1,0.9,0.4);
    electrons_leg_hydrogen->SetFillColor(0);


    TLegend* hydrogen_para1_leg = new TLegend(0.7,0.2,0.9,0.4);
    hydrogen_para1_leg->SetFillColor(0);

    TLegend* hydrogen_para2_leg = new TLegend(0.7,0.2,0.9,0.4);
    hydrogen_para2_leg->SetFillColor(0);

    TLegend* hydrogen_para3_leg = new TLegend(0.7,0.2,0.9,0.4);
    hydrogen_para3_leg->SetFillColor(0);

    TLegend* hydrogen_para4_leg = new TLegend(0.7,0.2,0.9,0.4);
    hydrogen_para4_leg->SetFillColor(0);

    TLegend* hydrogen_para5_leg = new TLegend(0.7,0.2,0.9,0.4);
    hydrogen_para5_leg->SetFillColor(0);

    TLegend* hydrogen_para6_leg = new TLegend(0.7,0.2,0.9,0.4);
    hydrogen_para6_leg->SetFillColor(0);

    TLegend* hydrogen_para7_leg = new TLegend(0.7,0.2,0.9,0.4);
    hydrogen_para7_leg->SetFillColor(0);

    TLegend* hydrogen_para8_leg = new TLegend(0.7,0.2,0.9,0.4);
    hydrogen_para8_leg->SetFillColor(0);

    TLegend* hydrogen_para9_leg = new TLegend(0.7,0.2,0.9,0.4);
    hydrogen_para9_leg->SetFillColor(0);

    TLegend* hydrogen_para10_leg = new TLegend(0.7,0.2,0.9,0.4);
    hydrogen_para10_leg->SetFillColor(0);

    TLegend* hydrogen_para11_leg = new TLegend(0.7,0.2,0.9,0.4);
    hydrogen_para11_leg->SetFillColor(0);

    TLegend* hydrogen_para12_leg = new TLegend(0.7,0.2,0.9,0.4);
    hydrogen_para12_leg->SetFillColor(0);

    TLegend* hydrogen_para13_leg = new TLegend(0.7,0.2,0.9,0.4);
    hydrogen_para13_leg->SetFillColor(0);

    TLegend* hydrogen_para14_leg = new TLegend(0.7,0.2,0.9,0.4);
    hydrogen_para14_leg->SetFillColor(0);



    TLegend* water_para1_leg = new TLegend(0.7,0.2,0.9,0.4);
    water_para1_leg->SetFillColor(0);

    TLegend* water_para2_leg = new TLegend(0.7,0.2,0.9,0.4);
    water_para2_leg->SetFillColor(0);

    TLegend* water_para3_leg = new TLegend(0.7,0.2,0.9,0.4);
    water_para3_leg->SetFillColor(0);

    TLegend* water_para4_leg = new TLegend(0.7,0.2,0.9,0.4);
    water_para4_leg->SetFillColor(0);

    TLegend* water_para5_leg = new TLegend(0.7,0.2,0.9,0.4);
    water_para5_leg->SetFillColor(0);

    TLegend* water_para6_leg = new TLegend(0.7,0.2,0.9,0.4);
    water_para6_leg->SetFillColor(0);

    TLegend* water_para7_leg = new TLegend(0.7,0.2,0.9,0.4);
    water_para7_leg->SetFillColor(0);

    TLegend* water_para8_leg = new TLegend(0.7,0.2,0.9,0.4);
    water_para8_leg->SetFillColor(0);

    TLegend* water_para9_leg = new TLegend(0.7,0.2,0.9,0.4);
    water_para9_leg->SetFillColor(0);

    TLegend* water_para10_leg = new TLegend(0.7,0.2,0.9,0.4);
    water_para10_leg->SetFillColor(0);

    TLegend* water_para11_leg = new TLegend(0.7,0.2,0.9,0.4);
    water_para11_leg->SetFillColor(0);

    TLegend* water_para12_leg = new TLegend(0.7,0.2,0.9,0.4);
    water_para12_leg->SetFillColor(0);

    TLegend* water_para13_leg = new TLegend(0.7,0.2,0.9,0.4);
    water_para13_leg->SetFillColor(0);

    TLegend* water_para14_leg = new TLegend(0.7,0.2,0.9,0.4);
    water_para14_leg->SetFillColor(0);


    TLegend* uranium_para1_leg = new TLegend(0.7,0.2,0.9,0.4);
    uranium_para1_leg->SetFillColor(0);

    TLegend* uranium_para2_leg = new TLegend(0.7,0.2,0.9,0.4);
    uranium_para2_leg->SetFillColor(0);

    TLegend* uranium_para3_leg = new TLegend(0.7,0.2,0.9,0.4);
    uranium_para3_leg->SetFillColor(0);

    TLegend* uranium_para4_leg = new TLegend(0.7,0.2,0.9,0.4);
    uranium_para4_leg->SetFillColor(0);

    TLegend* uranium_para5_leg = new TLegend(0.7,0.2,0.9,0.4);
    uranium_para5_leg->SetFillColor(0);

    TLegend* uranium_para6_leg = new TLegend(0.7,0.2,0.9,0.4);
    uranium_para6_leg->SetFillColor(0);

    TLegend* uranium_para7_leg = new TLegend(0.7,0.2,0.9,0.4);
    uranium_para7_leg->SetFillColor(0);

    TLegend* uranium_para8_leg = new TLegend(0.7,0.2,0.9,0.4);
    uranium_para8_leg->SetFillColor(0);

    TLegend* uranium_para9_leg = new TLegend(0.7,0.2,0.9,0.4);
    uranium_para9_leg->SetFillColor(0);

    TLegend* uranium_para10_leg = new TLegend(0.7,0.2,0.9,0.4);
    uranium_para10_leg->SetFillColor(0);

    TLegend* uranium_para11_leg = new TLegend(0.7,0.2,0.9,0.4);
    uranium_para11_leg->SetFillColor(0);

    TLegend* uranium_para12_leg = new TLegend(0.7,0.2,0.9,0.4);
    uranium_para12_leg->SetFillColor(0);

    TLegend* uranium_para13_leg = new TLegend(0.7,0.2,0.9,0.4);
    uranium_para13_leg->SetFillColor(0);

    TLegend* uranium_para14_leg = new TLegend(0.7,0.2,0.9,0.4);
    uranium_para14_leg->SetFillColor(0);

    stringstream leg_entry;

    for(unsigned int i = 0 ; i < graphs.size() ; i++)
    {
        switch (photo.at(i)->GetParametrization())
        {
            case ParametrizationType::PhotoKokoulinShadowBezrukovSoft:
                leg_entry<<"Kokoulin";
                break;
            case ParametrizationType::PhotoKokoulinShadowBezrukovHard:
                leg_entry<<"Kokoulin + hard";
                break;
            case ParametrizationType::PhotoRhodeShadowBezrukovSoft:
                leg_entry<<"Rhode";
                break;
            case ParametrizationType::PhotoRhodeShadowBezrukovHard:
                leg_entry<<"Rhode + hard";
                break;
            case ParametrizationType::PhotoBezrukovBugaevShadowBezrukovSoft:
                leg_entry<<"Bezrukov/Bugaev";
                break;
            case ParametrizationType::PhotoBezrukovBugaevShadowBezrukovHard:
                leg_entry<<"Bezrukov/Bugaev + hard";
                break;
            case ParametrizationType::PhotoZeusShadowBezrukovSoft:
                leg_entry<<"Zeus";
                break;
            case ParametrizationType::PhotoZeusShadowBezrukovHard:
                leg_entry<<"Zeus + hard";
                break;
            case ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowDutta:
                leg_entry<<"ALLM 91 shadow=Dutta";
                break;
            case ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowButkevich:
                leg_entry<<"ALLM 91 shadow=Butkevich";
                break;
            case ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowDutta:
                leg_entry<<"ALLM 97 shadow=Dutta";
                break;
            case ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich:
                leg_entry<<"ALLM 97 shadow=Butkevich";
                break;
            case ParametrizationType::PhotoButkevichMikhailovShadowDutta:
                leg_entry<<"Butkevich/Mikhailov shadow=Dutta";
                break;
            case ParametrizationType::PhotoButkevichMikhailovShadowButkevich:
                leg_entry<<"Butkevich/Mikhailovd shadow=Butkevich";
                break;
            default:
                log_fatal("Wrong Nuclear Interaction Parametrization Type '%i'."
                    , photo.at(i)->GetParametrization());
        }
        if(photo.at(i)->GetParticle()->GetType() == ParticleType::MuMinus)
        {
            if (photo.at(i)->GetMedium()->GetName().compare("water")==0)
            {
                muons_water_gr->Add(graphs.at(i),"l");
                graphs.at(i)->SetLineWidth(2);
                muons_leg_water->AddEntry(graphs.at(i),leg_entry.str().c_str(),"l");

                switch (muons_water_gr->GetListOfGraphs()->Capacity())
                {
                    case 1:
                        graphs.at(i)->SetLineColor(kRed);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 2:
                        graphs.at(i)->SetLineColor(kRed);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 3:
                        graphs.at(i)->SetLineColor(kGreen);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 4:
                        graphs.at(i)->SetLineColor(kGreen);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 5:
                        graphs.at(i)->SetLineColor(kBlue);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 6:
                        graphs.at(i)->SetLineColor(kBlue);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 7:
                        graphs.at(i)->SetLineColor(kBlack);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 8:
                        graphs.at(i)->SetLineColor(kBlack);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 9:
                        graphs.at(i)->SetLineColor(kOrange);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 10:
                        graphs.at(i)->SetLineColor(kOrange);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 11:
                        graphs.at(i)->SetLineColor(kMagenta);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 12:
                        graphs.at(i)->SetLineColor(kMagenta);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 13:
                        graphs.at(i)->SetLineColor(kCyan +2);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 14:
                        graphs.at(i)->SetLineColor(kCyan +2);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                }
            }
            else if(photo.at(i)->GetMedium()->GetName().compare("hydrogen")==0)
            {
                muons_hydrogen_gr->Add(graphs.at(i),"l");
                graphs.at(i)->SetLineWidth(2);
                muons_leg_hydrogen->AddEntry(graphs.at(i),leg_entry.str().c_str(),"l");

                switch (muons_hydrogen_gr->GetListOfGraphs()->Capacity())
                {
                    case 1:
                        graphs.at(i)->SetLineColor(kRed);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 2:
                        graphs.at(i)->SetLineColor(kRed);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 3:
                        graphs.at(i)->SetLineColor(kGreen);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 4:
                        graphs.at(i)->SetLineColor(kGreen);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 5:
                        graphs.at(i)->SetLineColor(kBlue);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 6:
                        graphs.at(i)->SetLineColor(kBlue);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 7:
                        graphs.at(i)->SetLineColor(kBlack);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 8:
                        graphs.at(i)->SetLineColor(kBlack);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 9:
                        graphs.at(i)->SetLineColor(kOrange);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 10:
                        graphs.at(i)->SetLineColor(kOrange);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 11:
                        graphs.at(i)->SetLineColor(kMagenta);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 12:
                        graphs.at(i)->SetLineColor(kMagenta);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 13:
                        graphs.at(i)->SetLineColor(kCyan +2);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 14:
                        graphs.at(i)->SetLineColor(kCyan +2);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                }
            }
            else if(photo.at(i)->GetMedium()->GetName().compare("uranium")==0)
            {
                muons_uranium_gr->Add(graphs.at(i),"l");
                graphs.at(i)->SetLineWidth(2);
                muons_leg_uranium->AddEntry(graphs.at(i),leg_entry.str().c_str(),"l");

                switch (muons_uranium_gr->GetListOfGraphs()->Capacity())
                {
                    case 1:
                        graphs.at(i)->SetLineColor(kRed);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 2:
                        graphs.at(i)->SetLineColor(kRed);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 3:
                        graphs.at(i)->SetLineColor(kGreen);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 4:
                        graphs.at(i)->SetLineColor(kGreen);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 5:
                        graphs.at(i)->SetLineColor(kBlue);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 6:
                        graphs.at(i)->SetLineColor(kBlue);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 7:
                        graphs.at(i)->SetLineColor(kBlack);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 8:
                        graphs.at(i)->SetLineColor(kBlack);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 9:
                        graphs.at(i)->SetLineColor(kOrange);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 10:
                        graphs.at(i)->SetLineColor(kOrange);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 11:
                        graphs.at(i)->SetLineColor(kMagenta);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 12:
                        graphs.at(i)->SetLineColor(kMagenta);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 13:
                        graphs.at(i)->SetLineColor(kCyan +2);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 14:
                        graphs.at(i)->SetLineColor(kCyan +2);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                }
            }
        }
        else if(photo.at(i)->GetParticle()->GetType() == ParticleType::TauMinus)
        {
            if (photo.at(i)->GetMedium()->GetName().compare("water")==0)
            {
                taus_water_gr->Add(graphs.at(i),"l");
                graphs.at(i)->SetLineWidth(2);
                taus_leg_water->AddEntry(graphs.at(i),leg_entry.str().c_str(),"l");

                switch (taus_water_gr->GetListOfGraphs()->Capacity())
                {
                    case 1:
                        graphs.at(i)->SetLineColor(kRed);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 2:
                        graphs.at(i)->SetLineColor(kRed);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 3:
                        graphs.at(i)->SetLineColor(kGreen);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 4:
                        graphs.at(i)->SetLineColor(kGreen);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 5:
                        graphs.at(i)->SetLineColor(kBlue);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 6:
                        graphs.at(i)->SetLineColor(kBlue);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 7:
                        graphs.at(i)->SetLineColor(kBlack);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 8:
                        graphs.at(i)->SetLineColor(kBlack);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 9:
                        graphs.at(i)->SetLineColor(kOrange);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 10:
                        graphs.at(i)->SetLineColor(kOrange);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 11:
                        graphs.at(i)->SetLineColor(kMagenta);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 12:
                        graphs.at(i)->SetLineColor(kMagenta);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 13:
                        graphs.at(i)->SetLineColor(kCyan +2);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 14:
                        graphs.at(i)->SetLineColor(kCyan +2);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                }
            }
            else if(photo.at(i)->GetMedium()->GetName().compare("hydrogen")==0)
            {
                taus_hydrogen_gr->Add(graphs.at(i),"l");
                graphs.at(i)->SetLineWidth(2);
                taus_leg_hydrogen->AddEntry(graphs.at(i),leg_entry.str().c_str(),"l");

                switch (taus_hydrogen_gr->GetListOfGraphs()->Capacity())
                {
                    case 1:
                        graphs.at(i)->SetLineColor(kRed);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 2:
                        graphs.at(i)->SetLineColor(kRed);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 3:
                        graphs.at(i)->SetLineColor(kGreen);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 4:
                        graphs.at(i)->SetLineColor(kGreen);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 5:
                        graphs.at(i)->SetLineColor(kBlue);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 6:
                        graphs.at(i)->SetLineColor(kBlue);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 7:
                        graphs.at(i)->SetLineColor(kBlack);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 8:
                        graphs.at(i)->SetLineColor(kBlack);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 9:
                        graphs.at(i)->SetLineColor(kOrange);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 10:
                        graphs.at(i)->SetLineColor(kOrange);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 11:
                        graphs.at(i)->SetLineColor(kMagenta);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 12:
                        graphs.at(i)->SetLineColor(kMagenta);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 13:
                        graphs.at(i)->SetLineColor(kCyan +2);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 14:
                        graphs.at(i)->SetLineColor(kCyan +2);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                }
            }
            else if(photo.at(i)->GetMedium()->GetName().compare("uranium")==0)
            {
                taus_uranium_gr->Add(graphs.at(i),"l");
                graphs.at(i)->SetLineWidth(2);
                taus_leg_uranium->AddEntry(graphs.at(i),leg_entry.str().c_str(),"l");

                switch (taus_uranium_gr->GetListOfGraphs()->Capacity())
                {
                    case 1:
                        graphs.at(i)->SetLineColor(kRed);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 2:
                        graphs.at(i)->SetLineColor(kRed);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 3:
                        graphs.at(i)->SetLineColor(kGreen);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 4:
                        graphs.at(i)->SetLineColor(kGreen);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 5:
                        graphs.at(i)->SetLineColor(kBlue);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 6:
                        graphs.at(i)->SetLineColor(kBlue);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 7:
                        graphs.at(i)->SetLineColor(kBlack);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 8:
                        graphs.at(i)->SetLineColor(kBlack);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 9:
                        graphs.at(i)->SetLineColor(kOrange);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 10:
                        graphs.at(i)->SetLineColor(kOrange);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 11:
                        graphs.at(i)->SetLineColor(kMagenta);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 12:
                        graphs.at(i)->SetLineColor(kMagenta);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 13:
                        graphs.at(i)->SetLineColor(kCyan +2);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 14:
                        graphs.at(i)->SetLineColor(kCyan +2);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                }
            }
        }
        else if(photo.at(i)->GetParticle()->GetType() == ParticleType::EMinus)
        {
            if (photo.at(i)->GetMedium()->GetName().compare("water")==0)
            {
                electrons_water_gr->Add(graphs.at(i),"l");
                graphs.at(i)->SetLineWidth(2);
                electrons_leg_water->AddEntry(graphs.at(i),leg_entry.str().c_str(),"l");

                switch (electrons_water_gr->GetListOfGraphs()->Capacity())
                {
                    case 1:
                        graphs.at(i)->SetLineColor(kRed);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 2:
                        graphs.at(i)->SetLineColor(kRed);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 3:
                        graphs.at(i)->SetLineColor(kGreen);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 4:
                        graphs.at(i)->SetLineColor(kGreen);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 5:
                        graphs.at(i)->SetLineColor(kBlue);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 6:
                        graphs.at(i)->SetLineColor(kBlue);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 7:
                        graphs.at(i)->SetLineColor(kBlack);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 8:
                        graphs.at(i)->SetLineColor(kBlack);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 9:
                        graphs.at(i)->SetLineColor(kOrange);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 10:
                        graphs.at(i)->SetLineColor(kOrange);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 11:
                        graphs.at(i)->SetLineColor(kMagenta);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 12:
                        graphs.at(i)->SetLineColor(kMagenta);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 13:
                        graphs.at(i)->SetLineColor(kCyan +2);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 14:
                        graphs.at(i)->SetLineColor(kCyan +2);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                }
            }
            else if(photo.at(i)->GetMedium()->GetName().compare("hydrogen")==0)
            {
                electrons_hydrogen_gr->Add(graphs.at(i),"l");
                graphs.at(i)->SetLineWidth(2);
                electrons_leg_hydrogen->AddEntry(graphs.at(i),leg_entry.str().c_str(),"l");

                switch (electrons_hydrogen_gr->GetListOfGraphs()->Capacity())
                {
                    case 1:
                        graphs.at(i)->SetLineColor(kRed);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 2:
                        graphs.at(i)->SetLineColor(kRed);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 3:
                        graphs.at(i)->SetLineColor(kGreen);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 4:
                        graphs.at(i)->SetLineColor(kGreen);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 5:
                        graphs.at(i)->SetLineColor(kBlue);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 6:
                        graphs.at(i)->SetLineColor(kBlue);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 7:
                        graphs.at(i)->SetLineColor(kBlack);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 8:
                        graphs.at(i)->SetLineColor(kBlack);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 9:
                        graphs.at(i)->SetLineColor(kOrange);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 10:
                        graphs.at(i)->SetLineColor(kOrange);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 11:
                        graphs.at(i)->SetLineColor(kMagenta);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 12:
                        graphs.at(i)->SetLineColor(kMagenta);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 13:
                        graphs.at(i)->SetLineColor(kCyan +2);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 14:
                        graphs.at(i)->SetLineColor(kCyan +2);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                }
            }
            else if(photo.at(i)->GetMedium()->GetName().compare("uranium")==0)
            {
                electrons_uranium_gr->Add(graphs.at(i),"l");
                graphs.at(i)->SetLineWidth(2);
                electrons_leg_uranium->AddEntry(graphs.at(i),leg_entry.str().c_str(),"l");

                switch (electrons_uranium_gr->GetListOfGraphs()->Capacity())
                {
                    case 1:
                        graphs.at(i)->SetLineColor(kRed);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 2:
                        graphs.at(i)->SetLineColor(kRed);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 3:
                        graphs.at(i)->SetLineColor(kGreen);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 4:
                        graphs.at(i)->SetLineColor(kGreen);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 5:
                        graphs.at(i)->SetLineColor(kBlue);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 6:
                        graphs.at(i)->SetLineColor(kBlue);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 7:
                        graphs.at(i)->SetLineColor(kBlack);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 8:
                        graphs.at(i)->SetLineColor(kBlack);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 9:
                        graphs.at(i)->SetLineColor(kOrange);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 10:
                        graphs.at(i)->SetLineColor(kOrange);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 11:
                        graphs.at(i)->SetLineColor(kMagenta);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 12:
                        graphs.at(i)->SetLineColor(kMagenta);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                    case 13:
                        graphs.at(i)->SetLineColor(kCyan +2);
                        graphs.at(i)->SetLineStyle(1);
                        break;
                    case 14:
                        graphs.at(i)->SetLineColor(kCyan +2);
                        graphs.at(i)->SetLineStyle(2);
                        break;
                }
            }
        }
        leg_entry.str("");
        leg_entry.clear();
    }
    muons_water->cd();
    muons_water->SetLogx();
    muons_water->SetLogy();
    muons_water->SetGridx();
    muons_water->SetGridy();
    muons_water_gr->Draw("A");
    muons_water_gr->GetXaxis()->SetTitle("energy [MeV]");
    muons_water_gr->GetXaxis()->SetTitleOffset(1.1);
    muons_water_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
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
    muons_hydrogen_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
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
    muons_uranium_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
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
    taus_water_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
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
    taus_hydrogen_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
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
    taus_uranium_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
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
    electrons_water_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
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
    electrons_hydrogen_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
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
    electrons_uranium_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    electrons_leg_uranium->Draw("Same");
    electrons_uranium->Write();

    for(unsigned int i = 0 ; i < graphs.size() ; i++)
    {
        if(photo.at(i)->GetMedium()->GetName().compare("hydrogen")==0)
        {
            switch (photo.at(i)->GetParametrization())
            {
                case ParametrizationType::PhotoKokoulinShadowBezrukovSoft:
                    hydrogen_para1_gr->Add(graphs.at(i),"P");
                    hydrogen_para1_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (hydrogen_para1_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoKokoulinShadowBezrukovHard:
                    hydrogen_para2_gr->Add(graphs.at(i),"P");
                    hydrogen_para2_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (hydrogen_para2_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoRhodeShadowBezrukovSoft:
                    hydrogen_para3_gr->Add(graphs.at(i),"P");
                    hydrogen_para3_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (hydrogen_para3_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoRhodeShadowBezrukovHard:
                    hydrogen_para4_gr->Add(graphs.at(i),"P");
                    hydrogen_para4_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (hydrogen_para4_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoBezrukovBugaevShadowBezrukovSoft:
                    hydrogen_para5_gr->Add(graphs.at(i),"P");
                    hydrogen_para5_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (hydrogen_para5_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoBezrukovBugaevShadowBezrukovHard:
                    hydrogen_para6_gr->Add(graphs.at(i),"P");
                    hydrogen_para6_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (hydrogen_para6_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoZeusShadowBezrukovSoft:
                    hydrogen_para7_gr->Add(graphs.at(i),"P");
                    hydrogen_para7_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (hydrogen_para7_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoZeusShadowBezrukovHard:
                    hydrogen_para8_gr->Add(graphs.at(i),"P");
                    hydrogen_para8_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (hydrogen_para8_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowDutta:
                    hydrogen_para9_gr->Add(graphs.at(i),"P");
                    hydrogen_para9_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (hydrogen_para9_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowButkevich:
                    hydrogen_para10_gr->Add(graphs.at(i),"P");
                    hydrogen_para10_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (hydrogen_para10_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowDutta:
                    hydrogen_para11_gr->Add(graphs.at(i),"P");
                    hydrogen_para11_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (hydrogen_para11_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich:
                    hydrogen_para12_gr->Add(graphs.at(i),"P");
                    hydrogen_para12_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (hydrogen_para12_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoButkevichMikhailovShadowDutta:
                    hydrogen_para13_gr->Add(graphs.at(i),"P");
                    hydrogen_para13_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (hydrogen_para13_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoButkevichMikhailovShadowButkevich:
                    hydrogen_para14_gr->Add(graphs.at(i),"P");
                    hydrogen_para14_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (hydrogen_para14_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                default:
                    log_fatal("Wrong Nuclear Interaction Parametrization Type '%i'."
                        , photo.at(i)->GetParametrization());
            }
        }
        else if(photo.at(i)->GetMedium()->GetName().compare("water")==0)
        {
            switch (photo.at(i)->GetParametrization())
            {
                case ParametrizationType::PhotoKokoulinShadowBezrukovSoft:
                    water_para1_gr->Add(graphs.at(i),"P");
                    water_para1_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (water_para1_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoKokoulinShadowBezrukovHard:
                    water_para2_gr->Add(graphs.at(i),"P");
                    water_para2_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (water_para2_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoRhodeShadowBezrukovSoft:
                    water_para3_gr->Add(graphs.at(i),"P");
                    water_para3_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (water_para3_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoRhodeShadowBezrukovHard:
                    water_para4_gr->Add(graphs.at(i),"P");
                    water_para4_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (water_para4_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoBezrukovBugaevShadowBezrukovSoft:
                    water_para5_gr->Add(graphs.at(i),"P");
                    water_para5_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (water_para5_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoBezrukovBugaevShadowBezrukovHard:
                    water_para6_gr->Add(graphs.at(i),"P");
                    water_para6_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (water_para6_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoZeusShadowBezrukovSoft:
                    water_para7_gr->Add(graphs.at(i),"P");
                    water_para7_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (water_para7_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoZeusShadowBezrukovHard:
                    water_para8_gr->Add(graphs.at(i),"P");
                    water_para8_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (water_para8_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowDutta:
                    water_para9_gr->Add(graphs.at(i),"P");
                    water_para9_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (water_para9_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowButkevich:
                    water_para10_gr->Add(graphs.at(i),"P");
                    water_para10_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (water_para10_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowDutta:
                    water_para11_gr->Add(graphs.at(i),"P");
                    water_para11_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (water_para11_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich:
                    water_para12_gr->Add(graphs.at(i),"P");
                    water_para12_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (water_para12_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoButkevichMikhailovShadowDutta:
                    water_para13_gr->Add(graphs.at(i),"P");
                    water_para13_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (water_para13_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoButkevichMikhailovShadowButkevich:
                    water_para14_gr->Add(graphs.at(i),"P");
                    water_para14_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (water_para14_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                default:
                    log_fatal("Wrong Nuclear Interaction Parametrization Type '%i'."
                        , photo.at(i)->GetParametrization());
            }
        }
        else if(photo.at(i)->GetMedium()->GetName().compare("uranium")==0)
        {
            switch (photo.at(i)->GetParametrization())
            {
                case ParametrizationType::PhotoKokoulinShadowBezrukovSoft:
                    uranium_para1_gr->Add(graphs.at(i),"P");
                    uranium_para1_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (uranium_para1_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoKokoulinShadowBezrukovHard:
                    uranium_para2_gr->Add(graphs.at(i),"P");
                    uranium_para2_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (uranium_para2_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoRhodeShadowBezrukovSoft:
                    uranium_para3_gr->Add(graphs.at(i),"P");
                    uranium_para3_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (uranium_para3_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoRhodeShadowBezrukovHard:
                    uranium_para4_gr->Add(graphs.at(i),"P");
                    uranium_para4_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (uranium_para4_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoBezrukovBugaevShadowBezrukovSoft:
                    uranium_para5_gr->Add(graphs.at(i),"P");
                    uranium_para5_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (uranium_para5_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoBezrukovBugaevShadowBezrukovHard:
                    uranium_para6_gr->Add(graphs.at(i),"P");
                    uranium_para6_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (uranium_para6_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoZeusShadowBezrukovSoft:
                    uranium_para7_gr->Add(graphs.at(i),"P");
                    uranium_para7_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (uranium_para7_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoZeusShadowBezrukovHard:
                    uranium_para8_gr->Add(graphs.at(i),"P");
                    uranium_para8_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (uranium_para8_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowDutta:
                    uranium_para9_gr->Add(graphs.at(i),"P");
                    uranium_para9_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (uranium_para9_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowButkevich:
                    uranium_para10_gr->Add(graphs.at(i),"P");
                    uranium_para10_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (uranium_para10_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowDutta:
                    uranium_para11_gr->Add(graphs.at(i),"P");
                    uranium_para11_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (uranium_para11_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich:
                    uranium_para12_gr->Add(graphs.at(i),"P");
                    uranium_para12_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (uranium_para12_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoButkevichMikhailovShadowDutta:
                    uranium_para13_gr->Add(graphs.at(i),"P");
                    uranium_para13_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (uranium_para13_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                case ParametrizationType::PhotoButkevichMikhailovShadowButkevich:
                    uranium_para14_gr->Add(graphs.at(i),"P");
                    uranium_para14_leg->AddEntry(graphs.at(i),photo.at(i)->GetParticle()->GetName().c_str(),"p");
                    switch (uranium_para14_gr->GetListOfGraphs()->Capacity())
                    {
                        case 1:
                            graphs.at(i)->SetMarkerColor(kRed);
                            break;
                        case 2:
                            graphs.at(i)->SetMarkerColor(kBlue);
                            break;
                        case 3:
                            graphs.at(i)->SetMarkerColor(kGreen);
                            break;
                    }
                    break;
                default:
                    log_fatal("Wrong Nuclear Interaction Parametrization Type '%i'."
                        , photo.at(i)->GetParametrization());
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
    hydrogen_para1_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
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
    hydrogen_para2_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
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
    hydrogen_para3_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
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
    hydrogen_para4_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    hydrogen_para4_leg->Draw("Same");
    hydrogen_para4->Write();

    hydrogen_para5->cd();
    hydrogen_para5->SetLogx();
    hydrogen_para5->SetLogy();
    hydrogen_para5->SetGridx();
    hydrogen_para5->SetGridy();
    hydrogen_para5_gr->Draw("A");
    hydrogen_para5_gr->GetXaxis()->SetTitle("energy [MeV]");
    hydrogen_para5_gr->GetXaxis()->SetTitleOffset(1.1);
    hydrogen_para5_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    hydrogen_para5_leg->Draw("Same");
    hydrogen_para5->Write();

    hydrogen_para6->cd();
    hydrogen_para6->SetLogx();
    hydrogen_para6->SetLogy();
    hydrogen_para6->SetGridx();
    hydrogen_para6->SetGridy();
    hydrogen_para6_gr->Draw("A");
    hydrogen_para6_gr->GetXaxis()->SetTitle("energy [MeV]");
    hydrogen_para6_gr->GetXaxis()->SetTitleOffset(1.1);
    hydrogen_para6_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    hydrogen_para6_leg->Draw("Same");
    hydrogen_para6->Write();

    hydrogen_para7->cd();
    hydrogen_para7->SetLogx();
    hydrogen_para7->SetLogy();
    hydrogen_para7->SetGridx();
    hydrogen_para7->SetGridy();
    hydrogen_para7_gr->Draw("A");
    hydrogen_para7_gr->GetXaxis()->SetTitle("energy [MeV]");
    hydrogen_para7_gr->GetXaxis()->SetTitleOffset(1.1);
    hydrogen_para7_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    hydrogen_para7_leg->Draw("Same");
    hydrogen_para7->Write();

    hydrogen_para8->cd();
    hydrogen_para8->SetLogx();
    hydrogen_para8->SetLogy();
    hydrogen_para8->SetGridx();
    hydrogen_para8->SetGridy();
    hydrogen_para8_gr->Draw("A");
    hydrogen_para8_gr->GetXaxis()->SetTitle("energy [MeV]");
    hydrogen_para8_gr->GetXaxis()->SetTitleOffset(1.1);
    hydrogen_para8_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    hydrogen_para8_leg->Draw("Same");
    hydrogen_para8->Write();

    hydrogen_para9->cd();
    hydrogen_para9->SetLogx();
    hydrogen_para9->SetLogy();
    hydrogen_para9->SetGridx();
    hydrogen_para9->SetGridy();
    hydrogen_para9_gr->Draw("A");
    hydrogen_para9_gr->GetXaxis()->SetTitle("energy [MeV]");
    hydrogen_para9_gr->GetXaxis()->SetTitleOffset(1.1);
    hydrogen_para9_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    hydrogen_para9_leg->Draw("Same");
    hydrogen_para9->Write();

    hydrogen_para10->cd();
    hydrogen_para10->SetLogx();
    hydrogen_para10->SetLogy();
    hydrogen_para10->SetGridx();
    hydrogen_para10->SetGridy();
    hydrogen_para10_gr->Draw("A");
    hydrogen_para10_gr->GetXaxis()->SetTitle("energy [MeV]");
    hydrogen_para10_gr->GetXaxis()->SetTitleOffset(1.1);
    hydrogen_para10_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    hydrogen_para10_leg->Draw("Same");
    hydrogen_para10->Write();

    hydrogen_para11->cd();
    hydrogen_para11->SetLogx();
    hydrogen_para11->SetLogy();
    hydrogen_para11->SetGridx();
    hydrogen_para11->SetGridy();
    hydrogen_para11_gr->Draw("A");
    hydrogen_para11_gr->GetXaxis()->SetTitle("energy [MeV]");
    hydrogen_para11_gr->GetXaxis()->SetTitleOffset(1.1);
    hydrogen_para11_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    hydrogen_para11_leg->Draw("Same");
    hydrogen_para11->Write();

    hydrogen_para12->cd();
    hydrogen_para12->SetLogx();
    hydrogen_para12->SetLogy();
    hydrogen_para12->SetGridx();
    hydrogen_para12->SetGridy();
    hydrogen_para12_gr->Draw("A");
    hydrogen_para12_gr->GetXaxis()->SetTitle("energy [MeV]");
    hydrogen_para12_gr->GetXaxis()->SetTitleOffset(1.1);
    hydrogen_para12_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    hydrogen_para12_leg->Draw("Same");
    hydrogen_para12->Write();

    hydrogen_para13->cd();
    hydrogen_para13->SetLogx();
    hydrogen_para13->SetLogy();
    hydrogen_para13->SetGridx();
    hydrogen_para13->SetGridy();
    hydrogen_para13_gr->Draw("A");
    hydrogen_para13_gr->GetXaxis()->SetTitle("energy [MeV]");
    hydrogen_para13_gr->GetXaxis()->SetTitleOffset(1.1);
    hydrogen_para13_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    hydrogen_para13_leg->Draw("Same");
    hydrogen_para13->Write();

    hydrogen_para14->cd();
    hydrogen_para14->SetLogx();
    hydrogen_para14->SetLogy();
    hydrogen_para14->SetGridx();
    hydrogen_para14->SetGridy();
    hydrogen_para14_gr->Draw("A");
    hydrogen_para14_gr->GetXaxis()->SetTitle("energy [MeV]");
    hydrogen_para14_gr->GetXaxis()->SetTitleOffset(1.1);
    hydrogen_para14_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    hydrogen_para14_leg->Draw("Same");
    hydrogen_para14->Write();

    water_para1->cd();
    water_para1->SetLogx();
    water_para1->SetLogy();
    water_para1->SetGridx();
    water_para1->SetGridy();
    water_para1_gr->Draw("A");
    water_para1_gr->GetXaxis()->SetTitle("energy [MeV]");
    water_para1_gr->GetXaxis()->SetTitleOffset(1.1);
    water_para1_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
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
    water_para2_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
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
    water_para3_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
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
    water_para4_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    water_para4_leg->Draw("Same");
    water_para4->Write();

    water_para5->cd();
    water_para5->SetLogx();
    water_para5->SetLogy();
    water_para5->SetGridx();
    water_para5->SetGridy();
    water_para5_gr->Draw("A");
    water_para5_gr->GetXaxis()->SetTitle("energy [MeV]");
    water_para5_gr->GetXaxis()->SetTitleOffset(1.1);
    water_para5_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    water_para5_leg->Draw("Same");
    water_para5->Write();

    water_para6->cd();
    water_para6->SetLogx();
    water_para6->SetLogy();
    water_para6->SetGridx();
    water_para6->SetGridy();
    water_para6_gr->Draw("A");
    water_para6_gr->GetXaxis()->SetTitle("energy [MeV]");
    water_para6_gr->GetXaxis()->SetTitleOffset(1.1);
    water_para6_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    water_para6_leg->Draw("Same");
    water_para6->Write();

    water_para7->cd();
    water_para7->SetLogx();
    water_para7->SetLogy();
    water_para7->SetGridx();
    water_para7->SetGridy();
    water_para7_gr->Draw("A");
    water_para7_gr->GetXaxis()->SetTitle("energy [MeV]");
    water_para7_gr->GetXaxis()->SetTitleOffset(1.1);
    water_para7_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    water_para7_leg->Draw("Same");
    water_para7->Write();

    water_para8->cd();
    water_para8->SetLogx();
    water_para8->SetLogy();
    water_para8->SetGridx();
    water_para8->SetGridy();
    water_para8_gr->Draw("A");
    water_para8_gr->GetXaxis()->SetTitle("energy [MeV]");
    water_para8_gr->GetXaxis()->SetTitleOffset(1.1);
    water_para8_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    water_para8_leg->Draw("Same");
    water_para8->Write();

    water_para9->cd();
    water_para9->SetLogx();
    water_para9->SetLogy();
    water_para9->SetGridx();
    water_para9->SetGridy();
    water_para9_gr->Draw("A");
    water_para9_gr->GetXaxis()->SetTitle("energy [MeV]");
    water_para9_gr->GetXaxis()->SetTitleOffset(1.1);
    water_para9_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    water_para9_leg->Draw("Same");
    water_para9->Write();

    water_para10->cd();
    water_para10->SetLogx();
    water_para10->SetLogy();
    water_para10->SetGridx();
    water_para10->SetGridy();
    water_para10_gr->Draw("A");
    water_para10_gr->GetXaxis()->SetTitle("energy [MeV]");
    water_para10_gr->GetXaxis()->SetTitleOffset(1.1);
    water_para10_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    water_para10_leg->Draw("Same");
    water_para10->Write();

    water_para11->cd();
    water_para11->SetLogx();
    water_para11->SetLogy();
    water_para11->SetGridx();
    water_para11->SetGridy();
    water_para11_gr->Draw("A");
    water_para11_gr->GetXaxis()->SetTitle("energy [MeV]");
    water_para11_gr->GetXaxis()->SetTitleOffset(1.1);
    water_para11_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    water_para11_leg->Draw("Same");
    water_para11->Write();

    water_para12->cd();
    water_para12->SetLogx();
    water_para12->SetLogy();
    water_para12->SetGridx();
    water_para12->SetGridy();
    water_para12_gr->Draw("A");
    water_para12_gr->GetXaxis()->SetTitle("energy [MeV]");
    water_para12_gr->GetXaxis()->SetTitleOffset(1.1);
    water_para12_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    water_para12_leg->Draw("Same");
    water_para12->Write();

    water_para13->cd();
    water_para13->SetLogx();
    water_para13->SetLogy();
    water_para13->SetGridx();
    water_para13->SetGridy();
    water_para13_gr->Draw("A");
    water_para13_gr->GetXaxis()->SetTitle("energy [MeV]");
    water_para13_gr->GetXaxis()->SetTitleOffset(1.1);
    water_para13_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    water_para13_leg->Draw("Same");
    water_para13->Write();

    water_para14->cd();
    water_para14->SetLogx();
    water_para14->SetLogy();
    water_para14->SetGridx();
    water_para14->SetGridy();
    water_para14_gr->Draw("A");
    water_para14_gr->GetXaxis()->SetTitle("energy [MeV]");
    water_para14_gr->GetXaxis()->SetTitleOffset(1.1);
    water_para14_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    water_para14_leg->Draw("Same");
    water_para14->Write();

    uranium_para1->cd();
    uranium_para1->SetLogx();
    uranium_para1->SetLogy();
    uranium_para1->SetGridx();
    uranium_para1->SetGridy();
    uranium_para1_gr->Draw("A");
    uranium_para1_gr->GetXaxis()->SetTitle("energy [MeV]");
    uranium_para1_gr->GetXaxis()->SetTitleOffset(1.1);
    uranium_para1_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
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
    uranium_para2_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
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
    uranium_para3_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
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
    uranium_para4_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    uranium_para4_leg->Draw("Same");
    uranium_para4->Write();

    uranium_para5->cd();
    uranium_para5->SetLogx();
    uranium_para5->SetLogy();
    uranium_para5->SetGridx();
    uranium_para5->SetGridy();
    uranium_para5_gr->Draw("A");
    uranium_para5_gr->GetXaxis()->SetTitle("energy [MeV]");
    uranium_para5_gr->GetXaxis()->SetTitleOffset(1.1);
    uranium_para5_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    uranium_para5_leg->Draw("Same");
    uranium_para5->Write();

    uranium_para6->cd();
    uranium_para6->SetLogx();
    uranium_para6->SetLogy();
    uranium_para6->SetGridx();
    uranium_para6->SetGridy();
    uranium_para6_gr->Draw("A");
    uranium_para6_gr->GetXaxis()->SetTitle("energy [MeV]");
    uranium_para6_gr->GetXaxis()->SetTitleOffset(1.1);
    uranium_para6_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    uranium_para6_leg->Draw("Same");
    uranium_para6->Write();

    uranium_para7->cd();
    uranium_para7->SetLogx();
    uranium_para7->SetLogy();
    uranium_para7->SetGridx();
    uranium_para7->SetGridy();
    uranium_para7_gr->Draw("A");
    uranium_para7_gr->GetXaxis()->SetTitle("energy [MeV]");
    uranium_para7_gr->GetXaxis()->SetTitleOffset(1.1);
    uranium_para7_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    uranium_para7_leg->Draw("Same");
    uranium_para7->Write();

    uranium_para8->cd();
    uranium_para8->SetLogx();
    uranium_para8->SetLogy();
    uranium_para8->SetGridx();
    uranium_para8->SetGridy();
    uranium_para8_gr->Draw("A");
    uranium_para8_gr->GetXaxis()->SetTitle("energy [MeV]");
    uranium_para8_gr->GetXaxis()->SetTitleOffset(1.1);
    uranium_para8_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    uranium_para8_leg->Draw("Same");
    uranium_para8->Write();

    uranium_para9->cd();
    uranium_para9->SetLogx();
    uranium_para9->SetLogy();
    uranium_para9->SetGridx();
    uranium_para9->SetGridy();
    uranium_para9_gr->Draw("A");
    uranium_para9_gr->GetXaxis()->SetTitle("energy [MeV]");
    uranium_para9_gr->GetXaxis()->SetTitleOffset(1.1);
    uranium_para9_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    uranium_para9_leg->Draw("Same");
    uranium_para9->Write();

    uranium_para10->cd();
    uranium_para10->SetLogx();
    uranium_para10->SetLogy();
    uranium_para10->SetGridx();
    uranium_para10->SetGridy();
    uranium_para10_gr->Draw("A");
    uranium_para10_gr->GetXaxis()->SetTitle("energy [MeV]");
    uranium_para10_gr->GetXaxis()->SetTitleOffset(1.1);
    uranium_para10_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    uranium_para10_leg->Draw("Same");
    uranium_para10->Write();

    uranium_para11->cd();
    uranium_para11->SetLogx();
    uranium_para11->SetLogy();
    uranium_para11->SetGridx();
    uranium_para11->SetGridy();
    uranium_para11_gr->Draw("A");
    uranium_para11_gr->GetXaxis()->SetTitle("energy [MeV]");
    uranium_para11_gr->GetXaxis()->SetTitleOffset(1.1);
    uranium_para11_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    uranium_para11_leg->Draw("Same");
    uranium_para11->Write();

    uranium_para12->cd();
    uranium_para12->SetLogx();
    uranium_para12->SetLogy();
    uranium_para12->SetGridx();
    uranium_para12->SetGridy();
    uranium_para12_gr->Draw("A");
    uranium_para12_gr->GetXaxis()->SetTitle("energy [MeV]");
    uranium_para12_gr->GetXaxis()->SetTitleOffset(1.1);
    uranium_para12_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    uranium_para12_leg->Draw("Same");
    uranium_para12->Write();

    uranium_para13->cd();
    uranium_para13->SetLogx();
    uranium_para13->SetLogy();
    uranium_para13->SetGridx();
    uranium_para13->SetGridy();
    uranium_para13_gr->Draw("A");
    uranium_para13_gr->GetXaxis()->SetTitle("energy [MeV]");
    uranium_para13_gr->GetXaxis()->SetTitleOffset(1.1);
    uranium_para13_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    uranium_para13_leg->Draw("Same");
    uranium_para13->Write();

    uranium_para14->cd();
    uranium_para14->SetLogx();
    uranium_para14->SetLogy();
    uranium_para14->SetGridx();
    uranium_para14->SetGridy();
    uranium_para14_gr->Draw("A");
    uranium_para14_gr->GetXaxis()->SetTitle("energy [MeV]");
    uranium_para14_gr->GetXaxis()->SetTitleOffset(1.1);
    uranium_para14_gr->GetYaxis()->SetTitle("dN/dx [cm^{-1}]");
    uranium_para14_leg->Draw("Same");
    uranium_para14->Write();

    file->Close();

    return 0;
}
