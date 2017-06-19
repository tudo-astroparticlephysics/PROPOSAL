
// #include <iostream>
// #include <time.h>
// #include <fstream>

#include "TCanvas.h"
#include "TH2D.h"
#include "TTree.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TLegend.h"

#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/Output.h"
// #include "PROPOSAL/Bremsstrahlung.h"
// #include "PROPOSAL/Integral.h"
// #include "PROPOSAL/Medium.h"
// #include "PROPOSAL/Ionization.h"
// #include "PROPOSAL/Epairproduction.h"

using namespace PROPOSAL;

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
        <<"This is an example to plot the behavior of multi/dou muon events.\n"
        <<"An multi muon event contains 2 muons with half of the energy of the\n"
        <<"previous simulated muon.\n"
        <<"A ROOT-file named MultiMuon.root will be generated,\n"
        <<"which contains all Particles and secondarys generated.\n"
        <<"You have to call the app. with 4 parameters: \n"
        <<"ecut seed vcut #particles \n"
        <<"e.g. 500 42 -1 1000\n"
        <<"or with 6 parameters using the first few parameters and\n"
        <<"giving log10(energy_min/MeV) and log10(energy_max/MeV) in addition.\n"
        <<"------------------------------------------------------------------\n"
        <<endl;



    if(argc<5 && argc != 2)
    {
        cerr<<"Not enough parameters(" << argc-1 << ") given! Expected 4 or 6 parameters! \n Aborting application!";
        exit(1);
    }
    else if(argc == 6|| argc > 7)
    {
        cerr<<"Number of parameters(" << argc-1 << ") not supported! Expected 4 or 6 parameters! \n Aborting application!";
        exit(1);
    }

    stringstream ss;
    ss.precision(4);

    string OutputFile = "MultiMuon.root";
    double ecut = 500;
    int seed = 42;
    double vcut = -1;
    int statistic = 200;

    double EminLog10 = 5;
    double EmaxLog10 = 8;

    if(2 == argc)
    {
        statistic = atof(argv[1]);
    }

    if(argc==5)
    {
        ecut = atof(argv[1]);
        seed = atoi(argv[2]);
        vcut = atof(argv[3]);
        statistic = atoi(argv[4]);
    }

    if(argc>6)
    {
        EminLog10 = atof(argv[5]);
        EmaxLog10 = atof(argv[6]);
    }


    cerr << "Starting with: \t ecut=" << ecut;
    cerr << "\n \t\t seed=" << seed;
    cerr << "\n \t\t vcut=" << vcut;
    cerr << "\n \t\t statitic=" << statistic;
    cerr << "\n \t\t Emin=" << EminLog10;
    cerr << "\n \t\t Emax=" << EmaxLog10;
    cerr << "\n------------------------------";
    cerr << "\n------------------------------\n";

    ////////////////////////////////////////////////////////
    //Set up all the propagation variables
    //
    Propagator *pr = new Propagator("resources/configuration_IceOnly");
    Output::getInstance().EnableROOTOutput( OutputFile.c_str() );
    PROPOSALParticle* part;
    //
    ////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////
    //Set up variables for loop
    //
    double energy = pow(10,EmaxLog10);
    int ctr = (int)(EmaxLog10 - EminLog10)/(log10(2));
    energy *=2;
    //
    ////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////
    //Start the loop which generates N="statistic" particles with starting
    //energy and then again N="statistic" with half the energy and so forth.
    //
    for(;ctr > 0; ctr--)
    {
        energy /= 2;
        ss.str("");
        ss << "Propagation of particle with log10(energy/MeV)=" << log10(energy) << " starts!";
        ProgressBar* P = new ProgressBar(statistic,25);
        P->start(ss.str());
        for(int i =0;i<statistic;i++)
        {
            P->update();
            part = new PROPOSALParticle();
            part->SetProperties(ctr+ ctr*statistic,i,energy);
            pr->Propagate(part,SQRT2*1e2*1e3);
        }
        delete P;
    }

    Output::getInstance().Close();
    //
    ////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////
    //Start the loop which generates N="statistic" particles with starting
    //energy and then again N="statistic" with half the energy and so forth.
    //
    TFile f(OutputFile.c_str(),"update");
    TTree* TreeSec = (TTree*)f.Get("secondarys");
    TTree* TreePrim = (TTree*)f.Get("propagated_primarys");
    TTree* TreeInPrim = (TTree*)f.Get("primarys");

    double X,Y,Z;
    double CurrPartEnergySec;
    double energy_propagated;

    int ParentParticleSec,ParticleIdSec, ParentParticlePrim, ParticleIdPrim;

    double Distance,ELossHalf1, ELossHalf2;

    double EIn;
    TreeSec->SetBranchAddress("z",&Z);
    TreeSec->SetBranchAddress("y",&Y);
    TreeSec->SetBranchAddress("x",&X);
    TreeSec->SetBranchAddress("particle_id",&ParticleIdSec);
    TreeSec->SetBranchAddress("parent_particle_id",&ParentParticleSec);
    TreeSec->SetBranchAddress("current_primary_energy",&CurrPartEnergySec);

    TreePrim->SetBranchAddress("particle_id",&ParticleIdPrim);
    TreePrim->SetBranchAddress("parent_particle_id",&ParentParticlePrim);

    TreePrim->SetBranchAddress("energy",&energy_propagated);

    TreeInPrim->SetBranchAddress("energy",&EIn);
    long NMaxSec = TreeSec->GetEntries();
    long NMaxPrim = TreePrim->GetEntries();

    double MaxDistance = SQRT2*1e2*1e3;
    double HalfDistance = MaxDistance/2;

    cerr << "\n Secondarys: " << NMaxSec  << endl;

    TCanvas* can;
    TH1D* HistELossHalf1,*HistELossHalf2;
    vector<TH1D*> VecOfHistsHalf1;
    vector<TH1D*> VecOfHistsHalf2;
    vector<TCanvas*> VecOfCanvas;
    VecOfHistsHalf1.resize(0);
    VecOfHistsHalf2.resize(0);

    TreePrim->GetEntry(0);
    TreeInPrim->GetEntry(0);
    TreeSec->GetEntry(0);
    long ctrPrim = 0,    ctrSec = 0;
    double EInOld = -1;
    bool first = true;

    while(ctrPrim < NMaxPrim)
    {
        TreePrim->GetEntry(ctrPrim);
        TreeInPrim->GetEntry(ctrPrim);
        ctrPrim++;
//        cerr << "Particle" << endl;
        if(EInOld != EIn)
        {
            EInOld = EIn;
//            cerr << "Energy" << endl;
            if(!first)
            {
                HistELossHalf1->SetLineColor(kRed);
                HistELossHalf1->Draw();
                HistELossHalf1->SetStats(0);

                HistELossHalf2->SetLineColor(kBlue);
                HistELossHalf2->Draw("same");

                TLegend* Leg = new TLegend(0.2,0.7,0.4,0.9);
                Leg->AddEntry(HistELossHalf1,"First","L");
                Leg->AddEntry(HistELossHalf2,"Second","L");
                Leg->Draw("Same");


                can->Write();
                HistELossHalf1->Write();
                HistELossHalf2->Write();

                VecOfHistsHalf1.push_back(HistELossHalf1);
                VecOfHistsHalf2.push_back(HistELossHalf2);
                VecOfCanvas.push_back(can);

            }
            first=false;

            ss.str(""); ss << log10(EIn);
            can = new TCanvas("","");
            can->SetName(ss.str().c_str());
            ss.str(""); ss << "Muons with log10(EIn)=" << log10(EIn) << " and dist/2=" << HalfDistance;
            can->SetTitle(ss.str().c_str());
            ss.str(""); ss << log10(EIn)<<"First";
            HistELossHalf1 = new TH1D(ss.str().c_str(),"First",64,4,EmaxLog10);
            ss.str(""); ss << log10(EIn)<<"Second";
            HistELossHalf2 = new TH1D(ss.str().c_str(),"Second",64,4,EmaxLog10);
        }


        bool done=false;
        ELossHalf1=0;
        ELossHalf2=0;

        while(      ParentParticleSec == ParentParticlePrim
                &&  ParticleIdPrim == ParticleIdSec-1
                &&  NMaxSec > ctrSec)
        {
            Distance = sqrt( Z*Z + Y*Y + X*X);

            if(Distance>HalfDistance && !done)
            {
                ELossHalf1 = EIn - CurrPartEnergySec;
                ELossHalf2 = EIn - energy_propagated - ELossHalf1;
                done = true;
            }
            TreeSec->GetEntry(ctrSec);
            ctrSec++;

        }
        if(ELossHalf1<106)ELossHalf1=106;
        if(ELossHalf2<106)ELossHalf2=106;
        HistELossHalf1->Fill(log10(ELossHalf1));
        HistELossHalf2->Fill(log10(ELossHalf2));

    }

    vector<TH1D*> VecOfDouHistHalf1;
    vector<TH1D*> VecOfDouHistHalf2;

    vector<TH1D*> VecOfSingleHistHalf1;
    vector<TH1D*> VecOfSingleHistHalf2;

    vector<TCanvas*> VecOfFinalCanvas;
    double DouEnergyLoss;
    for(unsigned int k = 0; k < VecOfHistsHalf1.size() - 1; k++)
    {
        ss.str(""); ss << "Final_" << VecOfCanvas.at(k)->GetName();
        VecOfFinalCanvas.push_back( new TCanvas(ss.str().c_str(),ss.str().c_str()));

        VecOfDouHistHalf1.push_back( (TH1D*)VecOfHistsHalf1.at(k+1)->Clone("half1D") );
        VecOfDouHistHalf1.at(k)->Reset();
        VecOfDouHistHalf1.at(k)->SetMarkerStyle(20);
        VecOfDouHistHalf1.at(k)->SetMarkerColor(kRed);
        VecOfDouHistHalf1.at(k)->SetDrawOption("P");

        VecOfDouHistHalf2.push_back( (TH1D*)VecOfHistsHalf2.at(k+1)->Clone("half2D") );
        VecOfDouHistHalf2.at(k)->Reset();
        VecOfDouHistHalf2.at(k)->SetMarkerStyle(34);
        VecOfDouHistHalf2.at(k)->SetMarkerColor(kBlue);
        VecOfDouHistHalf2.at(k)->SetDrawOption("P");

        VecOfSingleHistHalf1.push_back( (TH1D*)VecOfHistsHalf1.at(k)->Clone("half1S") );
        VecOfSingleHistHalf1.at(k)->Reset();

        VecOfSingleHistHalf2.push_back( (TH1D*)VecOfHistsHalf2.at(k)->Clone("half1S") );
        VecOfSingleHistHalf2.at(k)->Reset();

        for(unsigned int i = 0; i< statistic/2; i++)
        {
            DouEnergyLoss = 0;
            DouEnergyLoss += pow(10,VecOfHistsHalf1.at(k+1)->GetRandom());
            DouEnergyLoss += pow(10,VecOfHistsHalf1.at(k+1)->GetRandom());
            DouEnergyLoss = log10(DouEnergyLoss);
            VecOfDouHistHalf1.at(k)->Fill(DouEnergyLoss);

            DouEnergyLoss = 0;
            DouEnergyLoss += pow(10,VecOfHistsHalf2.at(k+1)->GetRandom());
            DouEnergyLoss += pow(10,VecOfHistsHalf2.at(k+1)->GetRandom());
            DouEnergyLoss = log10(DouEnergyLoss);
            VecOfDouHistHalf2.at(k)->Fill(DouEnergyLoss);

            VecOfSingleHistHalf1.at(k)->Fill(VecOfHistsHalf1.at(k)->GetRandom());
            VecOfSingleHistHalf2.at(k)->Fill(VecOfHistsHalf2.at(k)->GetRandom());
        }
        VecOfDouHistHalf1.at(k)->Draw("P");
        VecOfDouHistHalf2.at(k)->Draw("SAME P");

        VecOfSingleHistHalf1.at(k)->Draw("SAME");
        VecOfSingleHistHalf2.at(k)->Draw("SAME");

        VecOfDouHistHalf1.at(k)->Write();
        VecOfDouHistHalf2.at(k)->Write();
        VecOfSingleHistHalf1.at(k)->Write();
        VecOfSingleHistHalf2.at(k)->Write();
        VecOfFinalCanvas.at(k)->Write();

    }

    f.Close();
}
