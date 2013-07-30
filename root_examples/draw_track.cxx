#include "TH3F.h"
#include "TCanvas.h"
#include "PROPOSAL/Propagator.h"
#include "TPolyMarker3D.h"
#include "TApplication.h"
#include "PROPOSAL/Output.h"
#include "TSystem.h"
#include "TStyle.h"

using namespace std;

int main(int argc, char **argv)
{
    gStyle->SetOptStat(0);

    TApplication theApp("Analysis", &argc, argv);
    double energy   =   9e6;
    int scale       =   1;

    double xmin =   1e20;
    double ymin =   1e20;
    double zmin =   1e20;
    double xmax =   -1e20;
    double ymax =   -1e20;
    double zmax =   -1e20;
    Particle *particle;
    Propagator *pr = new Propagator("resources/configuration");
    pr->set_seed(1234);

    vector<TPolyMarker3D*> marker;
    vector<Particle*> secondarys;

    string type;

    string s;


    cout<<"Enter the particle type (mu/tau/e)"<<endl;
    cin>>type;
    cout<<"Enter the particle energy [MeV]"<<endl;
    cin>>energy;
    particle = new Particle(type ,-200.,10.,-10.,90.,0.,0.,0.);

    particle->SetEnergy(energy);

    secondarys    =       pr->Propagate(particle);


    for(unsigned int i = 0; i<secondarys.size(); i++)
    {
        if(secondarys.at(i)->GetX() < xmin)
        {
            xmin    =   secondarys.at(i)->GetX();
        }
        if(secondarys.at(i)->GetX() > xmax)
        {
            xmax    =   secondarys.at(i)->GetX();
        }

        if(secondarys.at(i)->GetY() < ymin)
        {
            ymin    =   secondarys.at(i)->GetY();
        }
        if(secondarys.at(i)->GetY() > ymax)
        {
            ymax    =   secondarys.at(i)->GetY();
        }

        if(secondarys.at(i)->GetZ() < zmin)
        {
            zmin    =   secondarys.at(i)->GetZ();
        }
        if(secondarys.at(i)->GetZ() > zmax)
        {
            zmax    =   secondarys.at(i)->GetZ();
        }

    }

    if(xmin < 0)
        xmin    =   xmin +  0.1*xmin;
    else
        xmin    =   xmin -  0.1*xmin;

    if(ymin < 0)
        ymin    =   ymin +  0.1*ymin;
    else
        ymin    =   ymin -  0.1*ymin;

    if(zmin < 0)
        zmin    =   zmin +  0.1*zmin;
    else
        zmin    =   zmin -  0.1*zmin;

    if(xmax < 0)
        xmax    =   xmax -  0.1*xmax;
    else
        xmax    =   xmax +  0.1*xmax;

    if(ymax < 0)
        ymax    =   ymax -  0.1*ymax;
    else
        ymax    =   ymax +  0.1*ymax;

    if(zmax < 0)
        zmax    =   zmax -  0.1*zmax;
    else
        zmax    =   zmax +  0.1*zmax;

    if(xmin==0 && xmax ==0)
    {
        xmin    =   -1000;
        xmax    =   1000;
    }
    if(ymin==0 && ymax ==0)
    {
        ymin    =   -1000;
        ymax    =   1000;
    }
    if(zmin==0 && zmax ==0)
    {
        zmin    =   -1000;
        zmax    =   1000;
    }

    TH3F *h3 = new TH3F("track","track",10,xmin,xmax,10,ymin,ymax,10,zmin,zmax);

    for(unsigned int i = 0; i<secondarys.size(); i++)
    {
        if(secondarys.at(i)->GetName().compare("epair")==0)
        {
            TPolyMarker3D* m  = new TPolyMarker3D();
            m->SetPoint(0,secondarys.at(i)->GetX(),secondarys.at(i)->GetY(),secondarys.at(i)->GetZ());
            m->SetMarkerStyle(20);
            m->SetMarkerColor(kRed);
            m->SetMarkerSize(secondarys.at(i)->GetEnergy()/energy*400*scale);
            marker.push_back(m);
        }
        else if(secondarys.at(i)->GetName().compare("delta")==0)
        {
            TPolyMarker3D* m  = new TPolyMarker3D();
            m->SetPoint(0,secondarys.at(i)->GetX(),secondarys.at(i)->GetY(),secondarys.at(i)->GetZ());
            m->SetMarkerStyle(20);
            m->SetMarkerColor(kGreen);
            m->SetMarkerSize(secondarys.at(i)->GetEnergy()/energy*4000*scale);
            marker.push_back(m);
        }
        else if(secondarys.at(i)->GetName().compare("brems")==0)
        {
            TPolyMarker3D* m  = new TPolyMarker3D();
            m->SetPoint(0,secondarys.at(i)->GetX(),secondarys.at(i)->GetY(),secondarys.at(i)->GetZ());
            m->SetMarkerStyle(20);
            m->SetMarkerColor(kBlue);
            m->SetMarkerSize(secondarys.at(i)->GetEnergy()/energy*400*scale);
            marker.push_back(m);
        }
        else if(secondarys.at(i)->GetName().compare("munu")==0)
        {
            TPolyMarker3D* m  = new TPolyMarker3D();
            m->SetPoint(0,secondarys.at(i)->GetX(),secondarys.at(i)->GetY(),secondarys.at(i)->GetZ());
            m->SetMarkerStyle(20);
            m->SetMarkerColor(kBlack);
            m->SetMarkerSize(secondarys.at(i)->GetEnergy()/energy*4000*scale);
            marker.push_back(m);
        }
        else if(secondarys.at(i)->GetName().compare("e")==0)
        {
            TPolyMarker3D* m  = new TPolyMarker3D();
            m->SetPoint(0,secondarys.at(i)->GetX(),secondarys.at(i)->GetY(),secondarys.at(i)->GetZ());
            m->SetMarkerStyle(20);
            m->SetMarkerColor(kCyan);
            m->SetMarkerSize(secondarys.at(i)->GetEnergy()/energy*400000*scale);
            marker.push_back(m);
        }
    }
    h3->GetXaxis()->SetTitle("x [cm]");
    h3->GetYaxis()->SetTitle("y [cm]");
    h3->GetZaxis()->SetTitle("z [cm]");

    h3->GetXaxis()->SetTitleOffset(1.6);
    h3->GetYaxis()->SetTitleOffset(1.6);
    h3->GetZaxis()->SetTitleOffset(1.6);

    h3->Draw();

    for(unsigned int i = 0;i<marker.size();i++)
    {
        sleep(0.1);
        marker.at(i)->Draw("same");
        gPad->Modified();
        gPad->Update();
    }
    secondarys.clear();
    marker.clear();


    theApp.Run(kTRUE);

    return 0;
}
