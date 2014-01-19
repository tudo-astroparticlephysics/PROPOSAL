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
        <<"or with 6 parameters using the first few parameters and\n"
        <<"giving log10(energy_min/MeV) and log10(energy_max/MeV) in addition.\n"
        <<"------------------------------------------------------------------\n"
        <<endl;

    if(argc<5)
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
    }
    string vcut_o;
    string ecut_o;


    Medium *med = new Medium("ice",1.);
    EnergyCutSettings* cuts = new EnergyCutSettings(ecut,vcut);


    Propagator *pr = new Propagator(med,cuts,"mu","resources/tables",false,false,false,false,1);

    Output::getInstance().EnableASCIIOutput("TestOutput.txt");

    ProgressBar* P = new ProgressBar(statistic,100);
    P->start("Propagation Starts!");
    for(int i =0;i<statistic;i++)
    {
        P->update();
        pr->GetParticle()->SetProperties(0,i,1e6);
        pr->Propagate(SQRT2*1e2*1e3);


    }

    Output::getInstance().Close();

}


