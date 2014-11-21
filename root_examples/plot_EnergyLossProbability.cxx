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
#include "TProfile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TPaveStats.h"

class ProgressBar {
private:
    unsigned long _steps;
    unsigned long _currentCount;
    unsigned long _maxbarLength;
    unsigned long _nextStep;
    unsigned long _updateSteps;
    time_t _startTime;
    std::string stringTmpl;
    std::string arrow;

public:
    /// Initialize a ProgressBar with [steps] number of steps, updated at [updateSteps] intervalls
    ProgressBar(unsigned long steps = 0, unsigned long updateSteps = 100) :
            _steps(steps), _currentCount(0), _maxbarLength(10), _updateSteps(
                    updateSteps), _nextStep(1), _startTime(0) {
        if (_updateSteps > _steps)
            _updateSteps = _steps;
        arrow.append(">");
    }

    void start(const std::string &title) {
        _startTime = time(NULL);
        std::string s = ctime(&_startTime);
        s.erase(s.end() - 1, s.end());
        stringTmpl = "  Started ";
        stringTmpl.append(s);
        stringTmpl.append(" : [%-10s] %3i%%    %s: %02i:%02is %s\r");
        std::cout << title << std::endl;

    }
    /// update the progressbar
    /// should be called steps times in a loop
    void update() {
        _currentCount++;
        if (_currentCount == _nextStep || _currentCount == _steps) {
            _nextStep += long(_steps / float(_updateSteps));

            int percentage = int(100 * _currentCount / float(_steps));
            time_t currentTime = time(NULL);
            if (_currentCount < _steps) {
                int j = 0;
                if (arrow.size()
                        <= (_maxbarLength) * (_currentCount) / (_steps))
                    arrow.insert(0, "=");
                float tElapsed = currentTime - _startTime;
                float tToGo = (_steps - _currentCount) * tElapsed
                        / _currentCount;
                printf(stringTmpl.c_str(), arrow.c_str(), percentage,
                        "Finish in", int(tToGo / 60), int(tToGo) % 60, "");
                fflush(stdout);
            } else {
                float tElapsed = currentTime - _startTime;
                std::string s = " - Finished at ";
                s.append(ctime(&currentTime));
                char fs[255];
                sprintf(fs, "%c[%d;%dm Finished %c[%dm", 27, 1, 32, 27, 0);
                printf(stringTmpl.c_str(), fs, percentage, "Needed",
                        int(tElapsed / 60), int(tElapsed) % 60, s.c_str());
            }
        }
    }

    /// Mark the progressbar with an error
    void setError() {
        time_t currentTime = time(NULL);
        _currentCount++;
        float tElapsed = currentTime - _startTime;
        std::string s = " - Finished at ";
        s.append(ctime(&currentTime));
        char fs[255];
        sprintf(fs, "%c[%d;%dm  ERROR   %c[%dm", 27, 1, 31, 27, 0);
        printf(stringTmpl.c_str(), fs, _currentCount, "Needed",
                int(tElapsed / 60), int(tElapsed) % 60, s.c_str());
    }

};



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
        <<"or with 7 parameters using the first few parameters and\n"
        <<"giving log10(energy_min/MeV),log10(energy_max/MeV) and path. \n."
        <<"for the output file.\n"
        <<"------------------------------------------------------------------\n"
        <<endl;

    if(argc<5)
    {
        cerr<<"Not enough parameters(" << argc-1 << ") given! Expected 4 or 6 parameters! \n Aborting application!";
        exit(1);
    }
    else if(argc == 6|| argc > 8)
    {
        cerr<<"Number of parameters(" << argc-1 << ") not supported! Expected 4 or 6 parameters! \n Aborting application!";
        exit(1);
    }

    string path = "";



    stringstream ss;



    double ecut = atof(argv[1]);
    int seed = atoi(argv[2]);
    double vcut = atof(argv[3]);
    int statistic = atoi(argv[4]);
    double EminLog10 = -4;
    double EmaxLog10 = 14;
    if(argc>6)
    {
        EminLog10 = atof(argv[5]);
        EmaxLog10 = atof(argv[6]);
        path = argv[7];
    }
    string vcut_o;
    string ecut_o;

    ss<<path << "EnergyLosses.root";
    Medium *med = new Medium("ice",1.);
    EnergyCutSettings* cuts = new EnergyCutSettings(ecut,vcut);

    TFile *file     =   new TFile(ss.str().c_str(),"RECREATE");


    Propagator *pr = new Propagator(med,cuts,"mu","resources/tables",false,false,false,false);


    TH2D *ioniz = new TH2D("ioniz","Ionization;primary log10(energy/MeV);secondary log10(energy/MeV)",512,2,EmaxLog10,512,EminLog10,EmaxLog10);
    TH2D *brems = new TH2D("brems","Bremsstrahlung;primary log10(energy/MeV);secondary log10(energy/MeV)",512,2,EmaxLog10,512,EminLog10,EmaxLog10);
    TH2D *photo = new TH2D("photo","Photonuclear Effect;primary log10(energy/MeV);secondary log10(energy/MeV)",512,2,EmaxLog10,512,EminLog10,EmaxLog10);
    TH2D *epair = new TH2D("epair","e-Pair Production;primary log10(energy/MeV);secondary log10(energy/MeV)",512,2,EmaxLog10,512,EminLog10,EmaxLog10);

    double sec_energy;
    int ww;
    double energy;

    TTree *tree = new TTree("secondarys","secondarys");
    tree->SetDirectory(0);
    tree->Branch("sec_energy",&sec_energy,"sec_energy/D");
    tree->Branch("energy",&energy,"energy/D");
    tree->Branch("ww",&ww,"ww/I");


    pr->set_seed(seed);

    TH1D *hist_ioniz = new TH1D("hist_ioniz","Ionization;primary log10(energy/MeV);secondary log10(energy/MeV)",50,2,EmaxLog10);
    TH1D *hist_brems = new TH1D("hist_brems","Bremsstrahlung;primary log10(energy/MeV);secondary log10(energy/MeV)",50,2,EmaxLog10);
    TH1D *hist_photo = new TH1D("hist_photo","Photonuclear Effect;primary log10(energy/MeV);secondary log10(energy/MeV)",50,2,EmaxLog10);
    TH1D *hist_epair = new TH1D("hist_epair","e-Pair Production;primary log10(energy/MeV);secondary log10(energy/MeV)",50,2,EmaxLog10);

    ProgressBar* P = new ProgressBar(statistic,100);
    P->start("Propagation Starts!");
    for(int i =0;i<statistic;i++)
    {
        P->update();
        pr->GetParticle()->SetEnergy(pow(10,EmaxLog10));
        pr->GetParticle()->SetPropagatedDistance(0);
        pr->Propagate(1e20);
        vector<PROPOSALParticle*> bla (Output::getInstance().GetSecondarys());
        for(unsigned int k =0 ; k<bla.size();k++)
        {
            sec_energy  =   bla.at(k)->GetEnergy();
            energy  =   bla.at(k)->GetParentParticleEnergy();

            if(bla.at(k)->GetType()==-1002)
            {
                ww=1;
                ioniz->Fill(log10(energy),log10(sec_energy));
                hist_ioniz->Fill(log10(energy),sec_energy);

            }
            else if(bla.at(k)->GetType()==-1001)
            {
                ww=2;
                brems->Fill(log10(energy),log10(sec_energy));
                hist_brems->Fill(log10(energy),sec_energy);

            }
            else if(bla.at(k)->GetType()==-1004)
            {
                ww=3;
                photo->Fill(log10(energy),log10(sec_energy));
                hist_photo->Fill(log10(energy),sec_energy);

            }
            else if(bla.at(k)->GetType()==-1003)
            {
                ww=4;
                epair->Fill(log10(energy),log10(sec_energy));
                hist_epair->Fill(log10(energy),sec_energy);

            }

            tree->Fill();

        }
        Output::getInstance().ClearSecondaryVector();



    }

    TPaveStats *st;
    TCanvas *can = new TCanvas("stochatic","stochastic",1400,1050);

    can->Divide(2,2);
    can->cd(1);
    can->SetGridx();
    can->SetGridy();
    can->SetTickx();
    can->SetTicky();
    gPad->SetLogz(1);
    ioniz->Draw("colz");
    gPad->Update();
    st = (TPaveStats*)ioniz->FindObject("stats");
    st->SetOptStat(10);
    st->SetX1NDC(0.1);
    st->SetX2NDC(0.4);
    st->SetY1NDC(0.9);
    st->SetY2NDC(0.85);
    gPad->Update();

    can->cd(2);
    can->SetGridx();
    can->SetGridy();
    can->SetTickx();
    can->SetTicky();
    gPad->SetLogz(1);
    brems->Draw("colz");
    gPad->Update();
    st = (TPaveStats*)brems->FindObject("stats");
    st->SetOptStat(10);
    st->SetX1NDC(0.1);
    st->SetX2NDC(0.4);
    st->SetY1NDC(0.9);
    st->SetY2NDC(0.85);
    gPad->Update();

    can->cd(3);
    can->SetGridx();
    can->SetGridy();
    can->SetTickx();
    can->SetTicky();
    gPad->SetLogz(1);
    photo->Draw("colz");
    gPad->Update();
    st = (TPaveStats*)photo->FindObject("stats");
    st->SetOptStat(10);
    st->SetX1NDC(0.1);
    st->SetX2NDC(0.4);
    st->SetY1NDC(0.9);
    st->SetY2NDC(0.85);
    gPad->Update();

    can->cd(4);
    can->SetGridx();
    can->SetGridy();
    can->SetTickx();
    can->SetTicky();
    gPad->SetLogz(1);
    epair->Draw("colz");
    gPad->Update();
    st = (TPaveStats*)epair->FindObject("stats");
    st->SetOptStat(10);
    st->SetX1NDC(0.1);
    st->SetX2NDC(0.4);
    st->SetY1NDC(0.9);
    st->SetY2NDC(0.85);
    gPad->Update();

    can->Write();
    ioniz->Write();
    epair->Write();
    brems->Write();
    photo->Write();

    //tree->Write();



    TCanvas* can_tot = new TCanvas("total","Total Contribution",1);

    double content;
    for(int i = 1 ; i<hist_epair->GetNbinsX()+1;i++)
    {
        content = (hist_epair->GetBinContent(i)+1)/statistic;
        content = log10(content);
        hist_epair->SetBinContent(i,content);
    }
    hist_epair->SetStats(0);
    hist_epair->SetLineColor(4);
    hist_epair->SetLineWidth(2.);
    hist_epair->Draw("L");

    for(int i = 1 ; i<hist_ioniz->GetNbinsX()+1;i++)
    {
        content = (hist_ioniz->GetBinContent(i)+1)/statistic;
        content = log10(content);
        hist_ioniz->SetBinContent(i,content);
    }
    hist_ioniz->SetLineColor(1);
    hist_ioniz->SetLineWidth(2.);
    hist_ioniz->Draw("L Same");

    for(int i = 1 ; i<hist_brems->GetNbinsX()+1;i++)
    {
        content = (hist_brems->GetBinContent(i)+1)/statistic;
        content = log10(content);
        hist_brems->SetBinContent(i,content);
    }
    hist_brems->SetLineColor(2);
    hist_ioniz->SetLineWidth(2.);
    hist_brems->Draw("L same");

    for(int i = 1 ; i<hist_photo->GetNbinsX()+1;i++)
    {
        content = (hist_photo->GetBinContent(i)+1)/statistic;
        content = log10(content);
        hist_photo->SetBinContent(i,content);
    }
    hist_photo->SetLineColor(3);
    hist_ioniz->SetLineWidth(2.);
    hist_photo->Draw("L same");

    TLegend* leg = new TLegend(0.1,0.7,0.3,0.9);
    leg->SetFillColor(0);
    leg->AddEntry(hist_ioniz,"Ionization","L");
    leg->AddEntry(hist_brems,"Bremsstrahlung","L");
    leg->AddEntry(hist_photo,"Photonuc.","L");
    leg->AddEntry(hist_epair,"e-pair","L");

    leg->Draw();

    can_tot->Write();

        file->Close();
}


