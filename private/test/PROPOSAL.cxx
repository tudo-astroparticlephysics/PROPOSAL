#include <iostream>
#include <fstream>

#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/Output.h"
#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"
#include "PROPOSAL/methods.h"


using namespace std;
namespace po	= boost::program_options;

boost::program_options::options_description CreateOptions();

string input_file;
string output_file;
string table_path;
string configuration_file;
string particle_type;

int NumberOfParticles;
double Estart,Eend, Gamma;

vector<double>ERange;


void ShowVersion();
void ShowHelp(po::options_description &all);

bool SaveArgsToCacheFile(int &argc, char** argv);
char** LoadArgsFromCacheFile(int &argc, char** argv);

int main(int argc, char** argv)
{
////    double x = (1-1e-16);
////    cerr << "erfInv(" << x << "): " << erfInv(x*2*(1-0.5)) << endl;
////    return 1;
//    Propagator* propa = new Propagator("resources/configuration");
//    Output::getInstance().WriteDescriptionFile();

//    if(argc != 1)
//    {
//        SaveArgsToCacheFile(argc,argv);
//    }
//    else
//    {
//        argv = LoadArgsFromCacheFile(argc,argv);
//    }


//    po::options_description all = CreateOptions();

//    //parse cmd line
//    po::variables_map vm;
//    po::store( po::command_line_parser(argc, argv).options(all).run(), vm);

//    //print help message if wanted
//    if(vm.count("help") || argc == 1)
//    {
//        ShowHelp(all);
//        exit(0);
//    }
//    //print version message if wanted
//    else if(vm.count("version"))
//    {
//        ShowVersion();
//        exit(0);
//    }
//    //notifies globalVar
//    try
//    {
//        //set the variables
//        vm.notify();
//    }
//    catch (po::invalid_command_line_syntax &e)
//    {
//        std::cerr<<"Error: "<<e.what()<<"\n";
//        exit(1);
//    }

    Propagator* propa = new Propagator("resources/configuration");

    Output::getInstance().EnableROOTOutput("TestOutput.root");
//Output::getInstance().SetLoggingConfigurationFile("resources/log4cplus.conf");
    for(int i = 0; i< (int)(1e4) ; i++)
    {
        PROPOSALParticle* part = new PROPOSALParticle(i,i,"mu",0,0,0,0,0,0,0,0);
        part->SetEnergy(1e6);

        propa->Propagate(part);
        log_error("asdsadada");
        delete part;
    }

    Output::getInstance().Close();
}




boost::program_options::options_description CreateOptions()
{

    po::options_description general("General options");
    general.add_options()
        ("help,h",		"Shows the help text. (this)")
        ("version,v",	"Shows the version of the program.");


    po::options_description propagator("Propagator options");
    propagator.add_options()
        ("input,i",       po::value<string>(&input_file)->default_value("")   ,   "Input file (ASCII/ROOT)")
        ("output,o",      po::value<string>(&output_file)->default_value("")  ,   "Output file (ASCII/.root)")
        ("tables,t",      po::value<string>(&table_path)->default_value("")  ,   "Absolute path where the interpolation tables are saved. If empty the path in the config file is used.")
        ("config,c",      po::value<string>(&configuration_file)->required()  ,   "Configuration file.");

    po::options_description quick_setup("No input file options");
    quick_setup.add_options()
        ("Nparticles,N",       po::value<int>(&NumberOfParticles)->default_value(0)   ,   "Number of Particles")
        ("Gamma,G",      po::value<double>(&Gamma)->default_value(0)  ,   "Spectral index for the generation E^gamma")
        ("Erange,E",     po::value< vector<double> >(&ERange)->multitoken()  ,   "Energy Range which will be generated in Log10(E/MeV). (Min Max)")
        ("Ptype,P",      po::value<string>(&particle_type)->default_value("mu")  ,   "Particle Type");


    po::options_description all("All options");
        all.add(general);
        all.add(propagator);
        all.add(quick_setup);


    return all;

}

void ShowVersion()
{
    std::cout << "PROPOSAL version " << "2.0.3" << endl;
}

void ShowHelp(po::options_description &all)
{

    std::cout << all;
}

bool SaveArgsToCacheFile(int &argc, char** argv)
{
    ofstream out;
    out.open(".PROPOSAL_CMDLINE_CACHE");
    if(out.is_open() != true)
    {
        cout << "WARNING:\tCould not write cache file!" << endl;
        return false;

    }

    out << argc << endl;
    for(int i=1; i<argc ; i++)
    {
        out << argv[i] << endl;
    }
    out.close();
    return true;
}

char** LoadArgsFromCacheFile(int &argc, char** argv)
{
    ifstream in;
    in.open(".PROPOSAL_CMDLINE_CACHE");

    if(in.is_open() != true)
    {
        cout << "WARNING:\tCould not open cache file!" << endl;
        return argv;
    }

    in >> argc;

    char** ReadINargv = new char*[argc];
    ReadINargv[0] = argv[0];


    string tmpString="sad";
    for(int i=1; i<argc ; i++)
    {
        in >> tmpString;
        ReadINargv[i] = new char[tmpString.length() + 1];
        std::strcpy(ReadINargv[i], tmpString.c_str());
    }
    argv = ReadINargv;
    in.close();

    return ReadINargv;
}
