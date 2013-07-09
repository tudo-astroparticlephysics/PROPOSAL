#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TCanvas.h"
#include "PROPOSAL/Propagator.h"
#include "TFile.h"


using namespace std;

main()
{
    Propagator *pr = new Propagator("resources/configuration");

    double range;

    int number_of_particles =   1e4;

    TH1D *range_hist = new TH1D("rnage","Range distribution of muons with an energy of 9 TeV",64,0,1.2e6);
    TFile *file = new TFile("Range_distribution.root","RECREATE");


    for(int i =0 ;i< number_of_particles; i++)
    {
        Particle *p = new Particle("mu");
        p->SetEnergy(9e6);
        pr->SetParticle(p);

        range   =   pr->Propagate(p);

        range_hist->Fill(range);

    }
    range_hist->GetXaxis()->SetTitle("range [cm]");
    range_hist->GetYaxis()->SetTitle("#");
    range_hist->Write();
    file->Close();

	return 0;
}
