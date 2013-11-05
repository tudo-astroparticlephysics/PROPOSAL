#include <iostream>
#include <fstream>
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/Interpolant.h"
#include "PROPOSAL/Ionization.h"
#include "PROPOSAL/Epairproduction.h"
#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/ContinuousRandomization.h"
#include "PROPOSAL/Geometry.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/StandardNormal.h"

//Stuff for LOG4CPLUS
#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
#include <log4cplus/configurator.h>
#include <log4cplus/layout.h>
#include <iomanip>

#include <time.h>


using namespace log4cplus;
using namespace std;



int main(int argc, char** argv)
{
    /*Medium *med = new Medium("ice",1);
        EnergyCutSettings *ecuts = new EnergyCutSettings();


        Propagator *prop = new Propagator("resources/configuration_IceOnly");

        Particle *part[100];

        prop->set_seed(142);

        ofstream out;
        out.precision(5);
        out.open("punktwolke.txt");

        for(int n = 0; n < 100; n++)
        {
            part[n] = new Particle("mu");
            part[n]->SetEnergy(10e6);
            part[n]->SetPropagatedDistance(5e5);

            prop->SetParticle(part[n]);
            prop->Propagate(part[n], out);

            //Trigger für ROOT einen neuen Graphen zu erstellen
            out << 99999 << "\n";

        }

        out.close();
        */


    // Initialisierung Teilchen und Scattering-Objekt
    Particle *part = new Particle("mu");
    part->SetEnergy(1e5);


    std::cout.precision(16);

    ScatteringMoliere *mol = new ScatteringMoliere();
    ScatteringFirstOrder *high = new ScatteringFirstOrder();

    //Initialisierung Zeitmesser (N Zeitmessungen)
    const int N = 50;
    clock_t start, end;
    double cpu_time_used[N];
    double mean, stdabw;


    //Antareswasser 10 m

    Medium *med = new Medium("antares_water",1.);

    //Moliere
    mean = 0;

    for(int n = 0;n<N;n++)
    {
        start = clock();

        for(int i = 0; i < 1e5; i++) mol->Scatter(1e3,part,med);

        end = clock();
        cpu_time_used[n] = ((double) (end - start)) / CLOCKS_PER_SEC;

        mean += cpu_time_used[n];
    }

    //Mittelwert und Fehler des Mittelwerts
    mean /=N;

    stdabw = 0;

    for(int i = 0; i<N; i++)
    {
        stdabw += (cpu_time_used[i]-mean)*(cpu_time_used[i]-mean);
    }

    stdabw = sqrt(1./(N*(N-1))*stdabw);

    std::cout << "Antareswasser Moliere : " << mean << "+-" << stdabw << endl;

    //Highland
    mean = 0;

    for(int n = 0;n<N;n++)
    {
        start = clock();

        for(int i = 0; i < 1e5; i++) high->Scatter(1e3,part,med);

        end = clock();
        cpu_time_used[n] = ((double) (end - start)) / CLOCKS_PER_SEC;

        mean += cpu_time_used[n];
    }

    //Mittelwert und Fehler des Mittelwerts
    mean /=N;

    stdabw = 0;

    for(int i = 0; i<N; i++)
    {
        stdabw += (cpu_time_used[i]-mean)*(cpu_time_used[i]-mean);
    }

    stdabw = sqrt(1./(N*(N-1))*stdabw);

    std::cout << "Antareswasser Highland : " << mean << "+-" << stdabw << endl;



    //Eis 10 m
    med = new Medium("ice",1.);

    //Moliere
    mean = 0;

        for(int n = 0;n<N;n++)
        {
            start = clock();

            for(int i = 0; i < 1e5; i++) mol->Scatter(1e3,part,med);

            end = clock();
            cpu_time_used[n] = ((double) (end - start)) / CLOCKS_PER_SEC;

            mean += cpu_time_used[n];
        }

        //Mittelwert und Fehler des Mittelwerts
        mean /=N;

        stdabw = 0;

        for(int i = 0; i<N; i++)
        {
            stdabw += (cpu_time_used[i]-mean)*(cpu_time_used[i]-mean);
        }

        stdabw = sqrt(1./(N*(N-1))*stdabw);

        std::cout << "Eis Moliere : " << mean << "+-" << stdabw << endl;

        //Highland
        mean = 0;

        for(int n = 0;n<N;n++)
        {
            start = clock();

            for(int i = 0; i < 1e5; i++) high->Scatter(1e3,part,med);

            end = clock();
            cpu_time_used[n] = ((double) (end - start)) / CLOCKS_PER_SEC;

            mean += cpu_time_used[n];
        }

        //Mittelwert und Fehler des Mittelwerts
        mean /=N;

        stdabw = 0;

        for(int i = 0; i<N; i++)
        {
            stdabw += (cpu_time_used[i]-mean)*(cpu_time_used[i]-mean);
        }

        stdabw = sqrt(1./(N*(N-1))*stdabw);

        std::cout << "Eis Highland : " << mean << "+-" << stdabw << endl;


}


