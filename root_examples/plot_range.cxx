#include "TH1D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "PROPOSAL/Propagator.h"
#include "TFile.h"


using namespace std;

int main()
{

    cout<<"This program will propagte 1e6 muons with an energy of 9 TeV."<<endl;
    Propagator *pr = new Propagator("resources/configuration");
    pr->set_seed(1234);
    pr->RandomDouble();
    double range;

    int number_of_particles =   1e5;
    int OnePercent = number_of_particles/100;
    TH1D *range_hist = new TH1D("Range","Range distribution of muons with an energy of 9 TeV",64,0,1.2e6);
    TFile *file = new TFile("Range_distribution.root","RECREATE");
    TTree *tree_range = new TTree("RangeTree","Range of the particles in a tree");
    tree_range->Branch("range",&range);


    for(int i =0 ;i< number_of_particles; i++)
    {
        if(i%OnePercent==0)cout<<"Progress: " <<100.*i/number_of_particles<<"%"<<endl;

        Particle *p = new Particle("mu");
        p->SetEnergy(9e6);
        pr->SetParticle(p);

        pr->Propagate(p);
        range = p->GetPropagatedDistance();

        range_hist->Fill(range);
        tree_range->Fill();

    }
    range_hist->GetXaxis()->SetTitle("range [cm]");
    range_hist->GetYaxis()->SetTitle("#");
    range_hist->Write();

    tree_range->Print();
    file->Write();
    file->Close();

	return 0;
}
