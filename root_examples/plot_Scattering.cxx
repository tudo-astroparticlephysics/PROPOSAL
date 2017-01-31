#include "TH1D.h"
#include "TCanvas.h"
#include "TBrowser.h"
#include "TTree.h"
#include "PROPOSAL/Propagator.h"
#include "TFile.h"
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include "TGraphErrors.h"


using namespace std;

int main()
{

    Medium* med = new Medium("ice",1.);
    EnergyCutSettings* ecut = new EnergyCutSettings(10,1e-3);
    Propagator* prop = new Propagator(med,ecut,"mu","resources/tables");//,false,true,true,true,1,12,1,1,1,1,false,2);

    TFile hfile("Decay_Reached_Scattering.root","RECREATE","Scattering_Info");

    double dist;
    double dev;
    double energy;
    double angle;
    double distanceMean=0,angleMean=0,devMean=0,energyMean=0;

    int N = (int)1e4;

    TTree *tree_decay = new TTree("Decay","An example of ROOT tree with a few branches");
    tree_decay->Branch("energy",&energy);
    tree_decay->Branch("distance",&dist);
    tree_decay->Branch("deviation",&dev);
    tree_decay->Branch("angle",&angle);

    TTree *tree_reached = new TTree("Reached","An example of ROOT tree with a few branches");
    tree_reached->Branch("energy",&energy);
    tree_reached->Branch("distance",&dist);
    tree_reached->Branch("deviation",&dev);
    tree_reached->Branch("angle",&angle);

    double propout;
    for(int i = 0; i<N ; i++)
    {
        prop->GetParticle()->SetEnergy(1e5);
        prop->GetParticle()->SetPhi(0);
        prop->GetParticle()->SetTheta(0);
        prop->GetParticle()->SetX(0);
        prop->GetParticle()->SetY(0);
        prop->GetParticle()->SetZ(0);
        prop->GetParticle()->SetT(0);
        prop->GetParticle()->SetPropagatedDistance(0);
        prop->GetParticle()->SetParticleId(0);

        propout = prop->Propagate(4.05e4);

        energy  =   prop->GetParticle()->GetEnergy();
        dist    =   prop->GetParticle()->GetPropagatedDistance();
        dev     =   sqrt(       prop->GetParticle()->GetX() * prop->GetParticle()->GetX()
                            +   prop->GetParticle()->GetY() * prop->GetParticle()->GetY() );
        angle   = asin(dev/dist)*180/PI;


        energyMean += energy;
        distanceMean += dist;
        devMean += dev;
        angleMean += angle;

        if(propout<0)
        {
            tree_decay->Fill();
        }
        else
        {
            tree_reached->Fill();
        }
    }
    tree_decay->Print();
    tree_reached->Print();
    hfile.Write();
    hfile.Close();

    distanceMean /= N;
    angleMean /= N;
    devMean /=N;
    energyMean /=N;
    cout.precision(16);
    cout << "Distance Mean: " << distanceMean << endl;
    cout << "Theta Mean: " << angleMean << endl;
    cout << "devMean: " << devMean << endl;
    cout << "energyMean: " << energyMean << endl;

	return 0;
}
