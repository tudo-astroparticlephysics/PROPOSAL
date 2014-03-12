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
        <<"This is an example to plot the stochasticity of energy losses.\n"
        <<"This means how often an energy loss occurs at a certain primary\n"
        <<"energy and with which energy.\n"
        <<"A ROOT-file named Stochasticity.root will be generated,\n"
        <<"which contains all Particles and secondarys generated.\n"
        <<"You have to call the app. with 4 parameters: \n"
        <<"ecut seed vcut #particles \n"
        <<"e.g. 500 42 -1 1000\n"
        <<"Either ecut or vcut can be set to -1 then they are ignored.\n"
        <<"As a 5th parameter you can pass a path where the output\n"
        <<"will be saved to.\n"
        <<"------------------------------------------------------------------\n"
        <<endl;

    stringstream ss;
    //Propagator *pr = new Propagator("resources/configuration_IceOnly");
    double ecut = atof(argv[1]);
    int seed = atoi(argv[2]);
    double vcut = atof(argv[3]);
    int statistic = atoi(argv[4]);

    string path = "";
    if(argc == 6)
    {
        path = argv[5];
    }

    string vcut_o;
    string ecut_o;

    ss<< path << "Stochasticity.root";
    Medium *med = new Medium("ice",1.);
    EnergyCutSettings* cuts = new EnergyCutSettings(ecut,vcut);

    TFile *file     =   new TFile(ss.str().c_str(),"RECREATE");


    Propagator *pr = new Propagator(med,cuts,"mu","resources/tables",false,false,false,false,1);



    TH2D *ioniz = new TH2D("ioniz","ioniz",512,2,14,512,-4,14);
    TH2D *brems = new TH2D("brems","brems",512,2,14,512,-4,14);
    TH2D *photo = new TH2D("photo","photo",512,2,14,512,-4,14);
    TH2D *epair = new TH2D("epair","epair",512,2,14,512,-4,14);

    double sec_energy;
    double prim_energy;
    int ww;
    double energy;

    TTree *tree = new TTree("secondarys","secondarys");
    tree->SetDirectory(0);
    tree->Branch("sec_energy",&sec_energy,"sec_energy/D");
    tree->Branch("energy",&energy,"energy/D");
    tree->Branch("ww",&ww,"ww/I");


    pr->set_seed(seed);

    ProgressBar* P = new ProgressBar(statistic,100);
    P->start("Starting to build stochasticity plot.");
    for(int i =0;i<statistic;i++)
    {
        P->update();
        pr->GetParticle()->SetEnergy(pow(10,13));
        pr->GetParticle()->SetPropagatedDistance(0);
        pr->Propagate(1e20);
        vector<Particle*> bla (Output::getInstance().GetSecondarys());
        for(unsigned int k =0 ; k<bla.size();k++)
        {
            sec_energy  =   bla.at(k)->GetEnergy();
            energy  =   bla.at(k)->GetParentParticleEnergy();
            if(bla.at(k)->GetType()==-1002)
            {
                ww=1;
                ioniz->Fill(log10(energy),log10(sec_energy));

            }
            else if(bla.at(k)->GetType()==-1001)
            {
                ww=2;
                brems->Fill(log10(energy),log10(sec_energy));

            }
            else if(bla.at(k)->GetType()==-1004)
            {
                ww=3;
                photo->Fill(log10(energy),log10(sec_energy));

            }
            else if(bla.at(k)->GetType()==-1003)
            {
                ww=4;
                epair->Fill(log10(energy),log10(sec_energy));

            }


//            if(Output::getInstance().GetSecondarys().at(k)->GetName().compare("delta")==0)
//            {
//                ww=1;
//                ioniz->Fill(log10(energy),log10(sec_energy));

//            }
//            else if(Output::getInstance().GetSecondarys().at(k)->GetName().compare("brems")==0)
//            {
//                ww=2;
//                brems->Fill(log10(energy),log10(sec_energy));

//            }
//            else if(Output::getInstance().GetSecondarys().at(k)->GetName().compare("munu")==0)
//            {
//                ww=3;
//                photo->Fill(log10(energy),log10(sec_energy));

//            }
//            else if(Output::getInstance().GetSecondarys().at(k)->GetName().compare("epair")==0)
//            {
//                ww=4;
//                epair->Fill(log10(energy),log10(sec_energy));

//            }
            tree->Fill();

        }
        Output::getInstance().ClearSecondaryVector();



    }

    TCanvas *can = new TCanvas("stochatic","stochastic",1400,1050);
    can->Divide(2,2);
    can->cd(1);
    can->SetGridx();
    can->SetGridy();
    can->SetTickx();
    can->SetTicky();
    ioniz->Draw("colz");

    can->cd(2);
    can->SetGridx();
    can->SetGridy();
    can->SetTickx();
    can->SetTicky();
    brems->Draw("colz");

    can->cd(3);
    can->SetGridx();
    can->SetGridy();
    can->SetTickx();
    can->SetTicky();
    photo->Draw("colz");

    can->cd(4);
    can->SetGridx();
    can->SetGridy();
    can->SetTickx();
    can->SetTicky();
    epair->Draw("colz");

    can->Write();
    ioniz->Write();
    epair->Write();
    brems->Write();
    photo->Write();

    tree->Write();

    file->Close();
}
