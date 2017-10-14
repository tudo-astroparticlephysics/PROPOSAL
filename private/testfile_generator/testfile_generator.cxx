
#include <fstream>
#include <sys/stat.h>
#include <random>

#include <boost/bind.hpp>

#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Check if Directory/File exists
// ------------------------------------------------------------------------- //

bool is_dir(std::string pathname)
{
    struct stat sb;

    if (stat(pathname.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool is_file(std::string pathname)
{
    struct stat sb;

    if (stat(pathname.c_str(), &sb) == 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// ------------------------------------------------------------------------- //
// Simple commandline parser
// ------------------------------------------------------------------------- //

class InputParser{
    public:
        InputParser (int &argc, const char **argv){
            for (int i=1; i < argc; ++i)
                this->tokens.push_back(std::string(argv[i]));
        }

        const std::string& getCmdOption(const std::string &option) const{
            std::vector<std::string>::const_iterator itr;
            itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
            if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                return *itr;
            }
            static const std::string empty_string("");
            return empty_string;
        }

        bool cmdOptionExists(const std::string &option) const{
            return std::find(this->tokens.begin(), this->tokens.end(), option)
                   != this->tokens.end();
        }
    private:
        std::vector <std::string> tokens;
};

// ------------------------------------------------------------------------- //
// Enums
// ------------------------------------------------------------------------- //

namespace Interaction{
    enum Enum {Bremsstrahlung, Epairproduction, Ionization, Photonuclear};
}

namespace Interpolation{
    enum Enum {None, Both, dEdx, dNdx};
}

namespace EvalutateFunction{
    enum Enum {dEdx, dNdx, dNdxRnd, StochasticLoss};
}

// ------------------------------------------------------------------------- //
// Helper
// ------------------------------------------------------------------------- //

class RndFromFile{
    private:
        double rnd_;
        std::string Path_;
        std::ifstream in_;

    public:
        RndFromFile(std::string Path){
            Path_ = Path;
            in_.open(Path_.c_str(), std::ifstream::in);
            in_>>rnd_;
            if(!in_.good())printf("less than one rnd_number!");
        }

        double rnd(){
            in_>>rnd_;
            if(!in_.good())
            {
                in_.close();
                in_.clear();
                in_.open(Path_.c_str());
                in_>>rnd_;
            }
            return rnd_;
        }
};

// ------------------------------------------------------------------------- //
// Random numbers
// ------------------------------------------------------------------------- //

int Create_RandomNumbers(std::string filename)
{
    std::ofstream out;
    out.open(filename, std::ofstream::out);

    if (!out.good())
    {
        printf("Cannot open file %s for writing\n", filename.c_str());
        return -1;
    }

    std::cout << "Create tables for Random Numbers" << std::endl;
    std::cout << "================================" << std::endl;

    std::mt19937 generator;
    std::uniform_real_distribution<double> random_double(0.0,1.0);

    // You know 6264 is cool, right?
    for (int i = 0; i < 6264; ++i)
    {
        out << random_double(generator) << std::endl;
    }

    // Everthing ok
    return 0;
}


// ------------------------------------------------------------------------- //
// CrossSections
// ------------------------------------------------------------------------- //


int setting_loop(
                 Interaction::Enum interaction,
                 EvalutateFunction::Enum eval_func,
                 Interpolation::Enum interpol,
                 std::string& filename,
                 std::string& header,
                 std::string path_to_tables,
                 std::string path_to_rnd)
{
    std::ofstream out;
    out.open(filename, std::ofstream::out);

    if (!out.good())
    {
        printf("Cannot open file %s for writing\n", filename.c_str());
        return -1;
    }

    // Header
    switch (interaction)
    {
        case Interaction::Ionization:
            out << "ecut"
                << "\t"
                << "vcut"
                << "\t"
                << "energy"
                << "\t"
                << "medium"
                << "\t"
                << "particle"
                << "\t"
                << header
                << std::endl;
            break;
        case Interaction::Epairproduction:
            out << "ecut"
                << "\t"
                << "vcut"
                << "\t"
                << "lpm"
                << "\t"
                << "energy"
                << "\t"
                << "medium"
                << "\t"
                << "particle"
                << "\t"
                << header << std::endl;
            break;
        case Interaction::Bremsstrahlung:
            out << "para"
                << "\t"
                << "ecut"
                << "\t"
                << "vcut"
                << "\t"
                << "lpm"
                << "\t"
                << "energy"
                << "\t"
                << "medium"
                << "\t"
                << "particle"
                << "\t"
                << header
                << std::endl;
            break;
        case Interaction::Photonuclear:
            out << "para"
                << "\t"
                << "ecut"
                << "\t"
                << "vcut"
                << "\t"
                << "energy"
                << "\t"
                << "medium"
                << "\t"
                << "particle"
                << "\t"
                << header
                << std::endl;
        default: // Do nothing
            break;
    }

    RndFromFile* rand_gen1 = NULL;
    RndFromFile* rand_gen2 = NULL;

    switch (eval_func)
    {
        case EvalutateFunction::dNdxRnd:
            rand_gen1 = new RndFromFile(path_to_rnd);
            break;
        case EvalutateFunction::StochasticLoss:
            rand_gen1 = new RndFromFile(path_to_rnd);
            rand_gen2 = new RndFromFile(path_to_rnd);
            rand_gen2->rnd(); // To ensure randoms are different
            break;
        default: // No random generator needed
            break;
    }

    std::vector<int> param;

    switch (interaction)
    {
        case Interaction::Epairproduction:
            param.push_back(1); // Just to have one loop in param for loop
            break;
        case Interaction::Ionization:
            param.push_back(1); // Just to have one loop in param for loop
            break;
        case Interaction::Bremsstrahlung:
            param.push_back(11);
            param.push_back(12);
            param.push_back(13);
            param.push_back(14);
            break;
        case Interaction::Photonuclear:
            param.push_back(31);
            param.push_back(-31);
            param.push_back(32);
            param.push_back(-32);
            param.push_back(33);
            param.push_back(-33);
            param.push_back(34);
            param.push_back(-34);
            param.push_back(35);
            param.push_back(-35);
            param.push_back(36);
            param.push_back(-36);
            param.push_back(37);
            param.push_back(-37);
            param.push_back(38);
            param.push_back(-38);
            break;
        default:
            break;
    }

    std::vector<double> energy;
    energy.push_back(10000);
    for (int i = 0; i < 9; i++)
    {
        energy.push_back(energy.at(i) * 10);
    }

    std::vector<std::string> medium;
    medium.push_back("ice"); medium.push_back("uranium"); medium.push_back("air");

    std::vector<std::string> particle;
    particle.push_back("mu"), particle.push_back("tau"), particle.push_back("e");

    std::vector<std::pair<int, double> > cuts;
    cuts.push_back(std::pair<int, double>(-1, -1));
    cuts.push_back(std::pair<int, double>(500, -1));
    cuts.push_back(std::pair<int, double>(-1, 0.05));
    cuts.push_back(std::pair<int, double>(500, 0.05));

    double aux;

    Vector3D particle_position(1.,1.,1.);
    Vector3D particle_direction(0,0,0);
    particle_direction.SetSphericalCoordinates(1, .20/180*PI, 20/180*PI);

    out.precision(16);

    // cuts
    for (std::vector<std::pair<int, double> >::iterator it_cut = cuts.begin() ; it_cut != cuts.end(); ++it_cut)
    {
        EnergyCutSettings cuts(it_cut->first, it_cut->second);

        // medium
        for (std::vector<std::string>::iterator it_medium = medium.begin() ; it_medium != medium.end(); ++it_medium)
        {
            Medium* medium = MediumFactory::Get()->CreateMedium(*it_medium);

            // paricle
            for (std::vector<std::string>::iterator it_particle = particle.begin() ; it_particle != particle.end(); ++it_particle)
            {
                Particle particle(Particle::GetTypeFromName(*it_particle),particle_position,particle_direction,1e5,10);

                CrossSections* cross = NULL;
                switch (interaction)
                {
                    case Interaction::Bremsstrahlung:
                        cross = new Bremsstrahlung(&particle, medium, &cuts);
                        break;
                    case Interaction::Epairproduction:
                        cross = new Epairproduction(&particle, medium, &cuts);
                        break;
                    case Interaction::Ionization:
                        cross = new Ionization(&particle, medium, &cuts);
                        break;
                    case Interaction::Photonuclear:
                        cross = new Photonuclear(&particle, medium, &cuts);
                        break;
                    default:
                        break;
                }

                // params
                for (std::vector<int>::iterator it_param = param.begin() ; it_param != param.end(); ++it_param)
                {
                    // switch (interpol)
                    // {
                    //     case Interpolation::dNdx:
                    //         cross->DisableDNdxInterpolation();
                    //         break;
                    //     case Interpolation::dEdx:
                    //         cross->DisableDEdxInterpolation();
                    //         break;
                    //     case Interpolation::None:
                    //         break;
                    //     case Interpolation::Both:
                    //         cross->DisableDNdxInterpolation();
                    //         cross->DisableDEdxInterpolation();
                    //         break;
                    //     default: // do nothing here
                    //         break;
                    // }

                    switch (interaction)
                    {
                        case Interaction::Bremsstrahlung: // no break
                        case Interaction::Photonuclear:
                            cross->SetParametrization(static_cast<ParametrizationType::Enum>(*it_param));
                            break;
                        default: // Do nothing
                            break;
                    }

                    // lpm
                    int lpm_limit = 2;
                    if (interaction == Interaction::Ionization || interaction == Interaction::Photonuclear)
                    {
                        lpm_limit = 1;
                    }

                    for (int lpm = 0; lpm < lpm_limit; ++lpm)
                    {
                        cross->EnableLpmEffect(lpm);

                        switch (interpol)
                        {
                            case Interpolation::dNdx:
                                cross->EnableDNdxInterpolation(path_to_tables);
                                break;
                            case Interpolation::dEdx:
                                cross->EnableDEdxInterpolation(path_to_tables);
                                break;
                            case Interpolation::Both:
                                cross->EnableDNdxInterpolation(path_to_tables);
                                cross->EnableDEdxInterpolation(path_to_tables);
                                break;
                            case Interpolation::None:
                                break;
                            default: // do nothing here
                                break;
                        }


                        // energy
                        for (std::vector<double>::iterator it_energy = energy.begin() ; it_energy != energy.end(); ++it_energy)
                        {
                            particle.SetEnergy(*it_energy);

                            switch (eval_func)
                            {
                                case EvalutateFunction::dEdx:
                                    aux = cross->CalculatedEdx();
                                    break;
                                case EvalutateFunction::dNdx:
                                    aux = cross->CalculatedNdx();
                                    break;
                                case EvalutateFunction::dNdxRnd:
                                    {
                                        double rnd1 = rand_gen1->rnd();
                                        aux = cross->CalculatedNdx(rnd1);
                                        break;
                                    }
                                case EvalutateFunction::StochasticLoss:
                                    {
                                        double rnd1 = rand_gen1->rnd();
                                        double rnd2 = rand_gen2->rnd();
                                        aux = cross->CalculateStochasticLoss(rnd1, rnd2);
                                        break;
                                    }
                                    break;
                                default:
                                    break;
                            }

                            switch (interaction)
                            {
                                case Interaction::Bremsstrahlung:
                                    out << *it_param
                                        << "\t"
                                        << it_cut->first
                                        << "\t"
                                        << it_cut->second
                                        << "\t"
                                        << lpm
                                        << "\t"
                                        << *it_energy
                                        << "\t"
                                        << *it_medium
                                        << "\t"
                                        << *it_particle
                                        << "\t"
                                        << aux
                                        << std::endl;
                                    break;
                                case Interaction::Epairproduction:
                                    out << it_cut->first
                                        << "\t"
                                        << it_cut->second
                                        << "\t"
                                        << lpm
                                        << "\t"
                                        << *it_energy
                                        << "\t"
                                        << *it_medium
                                        << "\t"
                                        << *it_particle
                                        << "\t"
                                        << aux
                                        << std::endl;
                                    break;
                                case Interaction::Ionization:
                                    out << it_cut->first
                                        << "\t"
                                        << it_cut->second
                                        << "\t"
                                        << *it_energy
                                        << "\t"
                                        << *it_medium
                                        << "\t"
                                        << *it_particle
                                        << "\t"
                                        << aux
                                        << std::endl;
                                    break;
                                case Interaction::Photonuclear:
                                    out << *it_param
                                        << "\t"
                                        << it_cut->first
                                        << "\t"
                                        << it_cut->second
                                        << "\t"
                                        << *it_energy
                                        << "\t"
                                        << *it_medium
                                        << "\t"
                                        << *it_particle
                                        << "\t"
                                        << aux
                                        << std::endl;
                                    break;
                                default:
                                    out << *it_param
                                        << "\t"
                                        << it_cut->first
                                        << "\t"
                                        << it_cut->second
                                        << "\t"
                                        << lpm
                                        << "\t"
                                        << *it_energy
                                        << "\t"
                                        << *it_medium
                                        << "\t"
                                        << *it_particle
                                        << "\t"
                                        << aux
                                        << std::endl;
                                    break;
                            }
                        }

                        switch (interpol)
                        {
                            case Interpolation::dNdx:
                                cross->DisableDNdxInterpolation();
                                if (cross->GetType() == ParticleType::NuclInt)
                                {
                                    static_cast<Photonuclear*>(cross)->DisablePhotoInterpolation();
                                }
                                if (cross->GetType() == ParticleType::EPair)
                                {
                                    static_cast<Epairproduction*>(cross)->DisableEpairInterpolation();
                                }
                                break;
                            case Interpolation::dEdx:
                                cross->DisableDEdxInterpolation();
                                if (cross->GetType() == ParticleType::NuclInt)
                                {
                                    static_cast<Photonuclear*>(cross)->DisablePhotoInterpolation();
                                }
                                if (cross->GetType() == ParticleType::EPair)
                                {
                                    static_cast<Epairproduction*>(cross)->DisableEpairInterpolation();
                                }
                                break;
                            case Interpolation::None:
                                break;
                            case Interpolation::Both:
                                cross->DisableDNdxInterpolation();
                                cross->DisableDEdxInterpolation();
                                if (cross->GetType() == ParticleType::NuclInt)
                                {
                                    static_cast<Photonuclear*>(cross)->DisablePhotoInterpolation();
                                }
                                if (cross->GetType() == ParticleType::EPair)
                                {
                                    static_cast<Epairproduction*>(cross)->DisableEpairInterpolation();
                                }
                                break;
                            default: // do nothing here
                                break;
                        }
                    }
                }

                delete cross;
                cross = NULL;
            }
        }
    }


    if (rand_gen1 != NULL)
    {
        delete rand_gen1;
    }
    if (rand_gen2 != NULL)
    {
        delete rand_gen2;
    }

    // Evering ok
    return 0;
}

// ------------------------------------------------------------------------- //
int Bremsstrahlung_Test_of_dEdx(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dEdx tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Bremsstrahlung, EvalutateFunction::dEdx, Interpolation::None, filename, header, path_to_tables, path_to_rnd);
}


// ------------------------------------------------------------------------- //
int Bremsstrahlung_Test_of_dNdx(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dNdx tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Bremsstrahlung, EvalutateFunction::dNdx, Interpolation::None, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Bremsstrahlung_Test_of_dNdxrnd(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dNdx with random number tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Bremsstrahlung, EvalutateFunction::dNdxRnd, Interpolation::None, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Bremsstrahlung_Test_of_e(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of StochasticLoss tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Bremsstrahlung, EvalutateFunction::StochasticLoss, Interpolation::None, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Bremsstrahlung_Test_of_dEdx_Interpolant(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dEdx Interpolant tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Bremsstrahlung, EvalutateFunction::dEdx, Interpolation::dEdx, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Bremsstrahlung_Test_of_dNdx_Interpolant(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dNdx Interpolant tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Bremsstrahlung, EvalutateFunction::dNdx, Interpolation::dNdx, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Bremsstrahlung_Test_of_dNdxrnd_Interpolant(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dNdxRnd Interpolant tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Bremsstrahlung, EvalutateFunction::dNdxRnd, Interpolation::dNdx, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Bremsstrahlung_Test_of_e_Interpolant(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of StochasticLoss tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Bremsstrahlung, EvalutateFunction::StochasticLoss, Interpolation::dNdx, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
// Epair
// ------------------------------------------------------------------------- //


int Epairproduction_Test_of_dEdx(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dEdx tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Epairproduction, EvalutateFunction::dEdx, Interpolation::None, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Epairproduction_Test_of_dNdx(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dNdx tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Epairproduction, EvalutateFunction::dNdx, Interpolation::None, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Epairproduction_Test_of_dNdxrnd(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dNdx with random number tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Epairproduction, EvalutateFunction::dNdxRnd, Interpolation::None, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Epairproduction_Test_of_e(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of StochasticLoss tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Epairproduction, EvalutateFunction::StochasticLoss, Interpolation::None, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Epairproduction_Test_of_dEdx_Interpolant(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dEdx Interpolant tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Epairproduction, EvalutateFunction::dEdx, Interpolation::dEdx, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Epairproduction_Test_of_dNdx_Interpolant(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dNdx Interpolant tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Epairproduction, EvalutateFunction::dNdx, Interpolation::dNdx, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Epairproduction_Test_of_dNdxrnd_Interpolant(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dNdxRnd Interpolant tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Epairproduction, EvalutateFunction::dNdxRnd, Interpolation::dNdx, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Epairproduction_Test_of_e_Interpolant(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of StochasticLoss tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Epairproduction, EvalutateFunction::StochasticLoss, Interpolation::Both, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
// Ionization
// ------------------------------------------------------------------------- //


int Ionization_Test_of_dEdx(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dEdx tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Ionization, EvalutateFunction::dEdx, Interpolation::None, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Ionization_Test_of_dNdx(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dNdx tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Ionization, EvalutateFunction::dNdx, Interpolation::None, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Ionization_Test_of_dNdxrnd(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dNdx with random number tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Ionization, EvalutateFunction::dNdxRnd, Interpolation::None, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Ionization_Test_of_e(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of StochasticLoss tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Ionization, EvalutateFunction::StochasticLoss, Interpolation::None, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Ionization_Test_of_dEdx_Interpolant(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dEdx Interpolant tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Ionization, EvalutateFunction::dEdx, Interpolation::dEdx, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Ionization_Test_of_dNdx_Interpolant(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dNdx Interpolant tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Ionization, EvalutateFunction::dNdx, Interpolation::dNdx, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Ionization_Test_of_dNdxrnd_Interpolant(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dNdxRnd Interpolant tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Ionization, EvalutateFunction::dNdxRnd, Interpolation::dNdx, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Ionization_Test_of_e_Interpolant(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of StochasticLoss tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Ionization, EvalutateFunction::StochasticLoss, Interpolation::dEdx, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
// Photonuclear
// ------------------------------------------------------------------------- //


int Photonuclear_Test_of_dEdx(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dEdx tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Photonuclear, EvalutateFunction::dEdx, Interpolation::None, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Photonuclear_Test_of_dNdx(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dNdx tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Photonuclear, EvalutateFunction::dNdx, Interpolation::None, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Photonuclear_Test_of_dNdxrnd(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dNdx with random number tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Photonuclear, EvalutateFunction::dNdxRnd, Interpolation::None, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Photonuclear_Test_of_e(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of StochasticLoss tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Photonuclear, EvalutateFunction::StochasticLoss, Interpolation::None, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Photonuclear_Test_of_dEdx_Interpolant(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dEdx Interpolant tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Photonuclear, EvalutateFunction::dEdx, Interpolation::dEdx, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Photonuclear_Test_of_dNdx_Interpolant(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dNdx Interpolant tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Photonuclear, EvalutateFunction::dNdx, Interpolation::dNdx, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Photonuclear_Test_of_dNdxrnd_Interpolant(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of dNdxRnd Interpolant tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Photonuclear, EvalutateFunction::dNdxRnd, Interpolation::dNdx, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
int Photonuclear_Test_of_e_Interpolant(std::string filename, std::string header, std::string path_to_tables, std::string path_to_rnd)
{
    printf("Generate Test of StochasticLoss tables in %s\n", filename.c_str());
    return setting_loop(Interaction::Photonuclear, EvalutateFunction::StochasticLoss, Interpolation::dNdx, filename, header, path_to_tables, path_to_rnd);
}

// ------------------------------------------------------------------------- //
// Continous Randomization
// ------------------------------------------------------------------------- //

int setting_loop_contrand(std::string& filename)
{
    std::ofstream out;
    out.open(filename, std::ofstream::out);

    if (!out.good())
    {
        printf("Cannot open file %s for writing\n", filename.c_str());
        return -1;
    }

    printf("Generate Test of ContinousRandomization Interpoltion tables in %s\n", filename.c_str());

    // Header
    out << "rnd"
        << "\t"
        << "paricle"
        << "\t"
        << "medium"
        << "\t"
        << "ecut"
        << "\t"
        << "vcut"
        << "\t"
        << "initial_energy"
        << "\t"
        << "final_energy"
        << "\t"
        << "randomized_energy"
        << "\t"
        << ""
        << std::endl;

    std::mt19937 generator;
    std::uniform_real_distribution<double> random_double(0.0,1.0);

    std::vector<double> energy;
    energy.push_back(10000);
    for (int i = 0; i < 9; i++)
    {
        energy.push_back(energy.at(i) * 10);
    }

    std::vector<std::string> medium;
    medium.push_back("ice"); medium.push_back("hydrogen"); medium.push_back("uranium");

    std::vector<std::string> particle;
    particle.push_back("mu"), particle.push_back("tau"), particle.push_back("e");

    std::vector<std::pair<int, double> > cuts;
    cuts.push_back(std::pair<int, double>(-1, -1));
    cuts.push_back(std::pair<int, double>(500, -1));
    cuts.push_back(std::pair<int, double>(1000, -1));
    cuts.push_back(std::pair<int, double>(-1, 0.05));
    cuts.push_back(std::pair<int, double>(500, 0.05));
    cuts.push_back(std::pair<int, double>(1000, 0.05));
    cuts.push_back(std::pair<int, double>(-1, 0.01));
    cuts.push_back(std::pair<int, double>(1000, 0.01));

    out.precision(16);

    Vector3D particle_position(1.,1.,1.);
    Vector3D particle_direction(0,0,0);
    particle_direction.SetSphericalCoordinates(1, .20/180*PI, 20/180*PI);

    // paricle
    for (std::vector<std::string>::iterator it_particle = particle.begin() ; it_particle != particle.end(); ++it_particle)
    {
        Particle particle(Particle::GetTypeFromName(*it_particle),particle_position,particle_direction,1e5,10);

        // medium
        for (std::vector<std::string>::iterator it_medium = medium.begin() ; it_medium != medium.end(); ++it_medium)
        {
            Medium* medium = MediumFactory::Get()->CreateMedium(*it_medium);

            // cuts
            for (std::vector<std::pair<int, double> >::iterator it_cut = cuts.begin() ; it_cut != cuts.end(); ++it_cut)
            {
                EnergyCutSettings cuts(it_cut->first, it_cut->second);

                // Propagator prop(&medium, &cuts, particle.GetType(), "");

                std::vector<CrossSections*> crosssections;

                crosssections.resize(4);

                // Ionization ioniz(&particle, &medium, &cuts);
                // Bremsstrahlung brems(&particle, &medium, &cuts);
                // Photonuclear photo(&particle, &medium, &cuts);
                // Epairproduction epair(&particle, &medium, &cuts);
                //
                // crosssections.push_back(&ioniz);
                // crosssections.push_back(&brems);
                // crosssections.push_back(&photo);
                // crosssections.push_back(&photo);

                crosssections.at(0) = new Ionization(&particle, medium, &cuts);
                crosssections.at(1) = new Bremsstrahlung(&particle, medium, &cuts);
                crosssections.at(2) = new Photonuclear(&particle, medium, &cuts);
                crosssections.at(3) = new Epairproduction(&particle, medium, &cuts);

                ContinuousRandomization cont(&particle, medium, crosssections);

                for(unsigned int i=0 ; i<crosssections.size();i++)
                {
                    crosssections.at(i)->EnableDEdxInterpolation();
                    crosssections.at(i)->EnableDNdxInterpolation();
                }
                cont.EnableDE2dxInterpolation();
                cont.EnableDE2deInterpolation();

                double rnd;

                // energy
                for (std::vector<double>::iterator it_energy = energy.begin() ; it_energy != energy.end(); ++it_energy)
                {

                    rnd = random_double(generator);

                    out << rnd
                        << "\t"
                        << *it_particle
                        << "\t"
                        << *it_medium
                        << "\t"
                        << it_cut->first
                        << "\t"
                        << it_cut->second
                        << "\t"
                        << *it_energy
                        << "\t"
                        << *it_energy - 1000
                        << "\t"
                        << cont.Randomize(*it_energy, *it_energy - 1000, rnd)
                        << std::endl;
                }

                for(unsigned int i = 0 ; i < crosssections.size() ; i++){

                    delete crosssections.at(i);
                }
            }
        }
    }

    // Evering ok
    return 0;
}

// ------------------------------------------------------------------------- //
// ProcessCollection
// ------------------------------------------------------------------------- //

int ProcColl_Stochastics(std::string filename, std::string path_to_tables)
{
    std::ofstream out;
    out.open(filename, std::ofstream::out);

    if (!out.good())
    {
        printf("Cannot open file %s for writing\n", filename.c_str());
        return -1;
    }

    printf("Generate Test of ProcessCollection Stochastics tables in %s\n", filename.c_str());

    // Header
    out << "ecut"
        << "\t"
        << "vcut"
        << "\t"
        << "lpm"
        << "\t"
        << "energy"
        << "\t"
        << "medium"
        << "\t"
        << "paricle"
        << "\t"
        << "IonizEvents"
        << "\t"
        << "BremsEvents"
        << "\t"
        << "PhotoEvents"
        << "\t"
        << "EpaiEvents"
        << std::endl;


    std::vector<double> energy;
    energy.push_back(10000);
    for (int i = 0; i < 9; i++)
    {
        energy.push_back(energy.at(i) * 10);
    }

    std::vector<std::string> medium;
    medium.push_back("ice"); medium.push_back("uranium"); medium.push_back("air");

    std::vector<std::string> particle;
    particle.push_back("mu"), particle.push_back("tau"), particle.push_back("e");

    std::vector<std::pair<int, double> > cuts;
    cuts.push_back(std::pair<int, double>(-1, -1));
    cuts.push_back(std::pair<int, double>(500, -1));
    cuts.push_back(std::pair<int, double>(-1, 0.05));
    cuts.push_back(std::pair<int, double>(500, 0.05));

    int lpm = 0;
    int brems_count = 0;
    int ioniz_count = 0;
    int epair_count = 0;
    int photo_count = 0;
    int statistics = 100000;
    std::pair<double, ParticleType::Enum> loss_return;

    out.precision(16);

    Vector3D particle_position(1.,1.,1.);
    Vector3D particle_direction(0,0,0);
    particle_direction.SetSphericalCoordinates(1, .20/180*PI, 20/180*PI);

    // cuts
    for (std::vector<std::pair<int, double> >::iterator it_cut = cuts.begin() ; it_cut != cuts.end(); ++it_cut)
    {
        EnergyCutSettings cuts(it_cut->first, it_cut->second);
        // medium
        for (std::vector<std::string>::iterator it_medium = medium.begin() ; it_medium != medium.end(); ++it_medium)
        {
            Medium* medium = MediumFactory::Get()->CreateMedium(*it_medium);

            // paricle
            for (std::vector<std::string>::iterator it_particle = particle.begin() ; it_particle != particle.end(); ++it_particle)
            {
                Particle particle(Particle::GetTypeFromName(*it_particle),particle_position,particle_direction,1e5,10);
                // Propagator prop(&medium, &cuts, particle.GetType(), "");
                ProcessCollection proc_col(&particle, medium, &cuts);
                proc_col.EnableInterpolation(path_to_tables);

                // energy
                for (std::vector<double>::iterator it_energy = energy.begin() ; it_energy != energy.end(); ++it_energy)
                {
                    brems_count = 0;
                    ioniz_count = 0;
                    epair_count = 0;
                    photo_count = 0;

                    for (int i = 0; i < statistics; ++i)
                    {
                        proc_col.GetParticle()->SetEnergy(*it_energy);
                        loss_return = proc_col.MakeStochasticLoss();
                        switch (loss_return.second)
                        {
                            case ParticleType::Brems:
                                brems_count++;
                                break;
                            case ParticleType::DeltaE:
                                ioniz_count++;
                                break;
                            case ParticleType::EPair:
                                epair_count++;
                                break;
                            case ParticleType::NuclInt:
                                photo_count++;
                                break;
                            default: // Do nothing
                                break;
                        }
                    }

                    out << it_cut->first
                        << "\t"
                        << it_cut->second
                        << "\t"
                        << lpm
                        << "\t"
                        << *it_energy
                        << "\t"
                        << *it_medium
                        << "\t"
                        << *it_particle
                        << "\t"
                        << ioniz_count
                        << "\t"
                        << brems_count
                        << "\t"
                        << photo_count
                        << "\t"
                        << epair_count
                        << std::endl;
                }
            }
        }
    }

    // Evering ok
    return 0;
}


int ProcColl_Displacement(std::string filename)
{
    std::ofstream out;
    out.open(filename, std::ofstream::out);

    if (!out.good())
    {
        printf("Cannot open file %s for writing\n", filename.c_str());
        return -1;
    }

    printf("Generate Test of ProcessCollection Displacement tables in %s\n", filename.c_str());

    // Header
    out << "ecut"
        << "\t"
        << "vcut"
        << "\t"
        << "lpm"
        << "\t"
        << "medium"
        << "\t"
        << "paricle"
        << "\t"
        << "energy"
        << "\t"
        << "ef"
        << "\t"
        << "dx"
        << std::endl;


    // std::vector<double> energy;
    // energy.push_back(10000);
    // for (int i = 0; i < 9; i++)
    // {
    //     energy.push_back(energy.at(i) * 10);
    // }

    std::vector<std::pair<double, double> > energy;
    energy.push_back(std::pair<double, double>(100000,         10000));
    energy.push_back(std::pair<double, double>(10000000,       10000));
    energy.push_back(std::pair<double, double>(1000000000,     10000));
    energy.push_back(std::pair<double, double>(100000000000,   10000));
    energy.push_back(std::pair<double, double>(10000000000000, 10000));
    energy.push_back(std::pair<double, double>(10000000,       1000000));
    energy.push_back(std::pair<double, double>(1000000000,     1000000));
    energy.push_back(std::pair<double, double>(100000000000,   1000000));
    energy.push_back(std::pair<double, double>(10000000000000, 1000000));
    energy.push_back(std::pair<double, double>(1000000000,     100000000));
    energy.push_back(std::pair<double, double>(100000000000,   100000000));
    energy.push_back(std::pair<double, double>(10000000000000, 100000000));
    energy.push_back(std::pair<double, double>(100000000000,   10000000000));
    energy.push_back(std::pair<double, double>(10000000000000, 10000000000));
    energy.push_back(std::pair<double, double>(10000000000000, 1000000000000));

    std::vector<std::string> medium;
    medium.push_back("ice"); medium.push_back("uranium"); medium.push_back("air");

    std::vector<std::string> particle;
    particle.push_back("mu"), particle.push_back("tau"), particle.push_back("e");

    std::vector<std::pair<int, double> > cuts;
    cuts.push_back(std::pair<int, double>(-1, -1));
    cuts.push_back(std::pair<int, double>(500, -1));
    cuts.push_back(std::pair<int, double>(-1, 0.05));
    cuts.push_back(std::pair<int, double>(500, 0.05));

    int lpm = 0;
    double dx = 0;
    double dist = 1.;

    out.precision(16);

    Vector3D particle_position(1.,1.,1.);
    Vector3D particle_direction(0,0,0);
    particle_direction.SetSphericalCoordinates(1, .20/180*PI, 20/180*PI);

    // cuts
    for (std::vector<std::pair<int, double> >::iterator it_cut = cuts.begin() ; it_cut != cuts.end(); ++it_cut)
    {
        EnergyCutSettings cuts(it_cut->first, it_cut->second);
        // medium
        for (std::vector<std::string>::iterator it_medium = medium.begin() ; it_medium != medium.end(); ++it_medium)
        {
            Medium* medium = MediumFactory::Get()->CreateMedium(*it_medium);

            // paricle
            for (std::vector<std::string>::iterator it_particle = particle.begin() ; it_particle != particle.end(); ++it_particle)
            {
                Particle particle(Particle::GetTypeFromName(*it_particle),particle_position,particle_direction,1e5,10);
                // Propagator prop(&medium, &cuts, particle.GetType(), "");
                ProcessCollection proc_col(&particle, medium, &cuts);
                proc_col.EnableInterpolation("../resources/tables");

                // energy
                for (std::vector<std::pair<double, double> >::iterator it_energy = energy.begin() ; it_energy != energy.end(); ++it_energy)
                {
                    proc_col.GetParticle()->SetEnergy((*it_energy).first);
                    dx = proc_col.CalculateDisplacement((*it_energy).first, (*it_energy).second, dist);

                    out << it_cut->first
                        << "\t"
                        << it_cut->second
                        << "\t"
                        << lpm
                        << "\t"
                        << *it_medium
                        << "\t"
                        << *it_particle
                        << "\t"
                        << (*it_energy).first
                        << "\t"
                        << (*it_energy).second
                        << "\t"
                        << dx
                        << std::endl;
                }
            }
        }
    }

    // Evering ok
    return 0;
}


int ProcColl_FinalEnergyDist(std::string filename, std::string path_to_tables)
{
    std::ofstream out;
    out.open(filename, std::ofstream::out);

    if (!out.good())
    {
        printf("Cannot open file %s for writing\n", filename.c_str());
        return -1;
    }

    printf("Generate Test of ProcessCollection Final Energy Displacement tables in %s\n", filename.c_str());

    // Header
    out << "ecut"
        << "\t"
        << "vcut"
        << "\t"
        << "lpm"
        << "\t"
        << "medium"
        << "\t"
        << "paricle"
        << "\t"
        << "dist"
        << "\t"
        << "energy"
        << "\t"
        << "ef"
        << "\t"
        << "FinalEnergy"
        << std::endl;

    std::vector<double> dist;
    dist.push_back(1);
    dist.push_back(10);
    dist.push_back(100);
    dist.push_back(1000);
    dist.push_back(10000);
    dist.push_back(100000);

    std::vector<std::pair<double, double> > energy;
    energy.push_back(std::pair<double, double>(100000,         10000));
    energy.push_back(std::pair<double, double>(10000000,       10000));
    energy.push_back(std::pair<double, double>(1000000000,     10000));
    energy.push_back(std::pair<double, double>(100000000000,   10000));
    energy.push_back(std::pair<double, double>(10000000000000, 10000));
    energy.push_back(std::pair<double, double>(10000000,       1000000));
    energy.push_back(std::pair<double, double>(1000000000,     1000000));
    energy.push_back(std::pair<double, double>(100000000000,   1000000));
    energy.push_back(std::pair<double, double>(10000000000000, 1000000));
    energy.push_back(std::pair<double, double>(1000000000,     100000000));
    energy.push_back(std::pair<double, double>(100000000000,   100000000));
    energy.push_back(std::pair<double, double>(10000000000000, 100000000));
    energy.push_back(std::pair<double, double>(100000000000,   10000000000));
    energy.push_back(std::pair<double, double>(10000000000000, 10000000000));
    energy.push_back(std::pair<double, double>(10000000000000, 1000000000000));

    std::vector<std::string> medium;
    medium.push_back("ice"); medium.push_back("uranium"); medium.push_back("air");

    std::vector<std::string> particle;
    particle.push_back("mu"), particle.push_back("tau"), particle.push_back("e");

    std::vector<std::pair<int, double> > cuts;
    cuts.push_back(std::pair<int, double>(-1, -1));
    cuts.push_back(std::pair<int, double>(500, -1));
    cuts.push_back(std::pair<int, double>(-1, 0.05));
    cuts.push_back(std::pair<int, double>(500, 0.05));

    int lpm = 0;
    double final_energy = 0;

    out.precision(16);

    Vector3D particle_position(1.,1.,1.);
    Vector3D particle_direction(0,0,0);
    particle_direction.SetSphericalCoordinates(1, .20/180*PI, 20/180*PI);

    // cuts
    for (std::vector<std::pair<int, double> >::iterator it_cut = cuts.begin() ; it_cut != cuts.end(); ++it_cut)
    {
        EnergyCutSettings cuts(it_cut->first, it_cut->second);
        // medium
        for (std::vector<std::string>::iterator it_medium = medium.begin() ; it_medium != medium.end(); ++it_medium)
        {
            Medium* medium = MediumFactory::Get()->CreateMedium(*it_medium);

            // paricle
            for (std::vector<std::string>::iterator it_particle = particle.begin() ; it_particle != particle.end(); ++it_particle)
            {
                Particle particle(Particle::GetTypeFromName(*it_particle),particle_position,particle_direction,1e5,10);
                // Propagator prop(&medium, &cuts, particle.GetType(), "");
                ProcessCollection proc_col(&particle, medium, &cuts);
                proc_col.EnableInterpolation(path_to_tables);

                for (std::vector<double>::iterator it_dist = dist.begin(); it_dist != dist.end(); ++it_dist)
                {
                    // energy
                    for (std::vector<std::pair<double, double> >::iterator it_energy = energy.begin() ; it_energy != energy.end(); ++it_energy)
                    {
                        proc_col.GetParticle()->SetEnergy((*it_energy).first);
                        proc_col.CalculateDisplacement((*it_energy).first, (*it_energy).second, *it_dist);
                        final_energy = proc_col.CalculateFinalEnergy((*it_energy).first, *it_dist);

                        out << it_cut->first
                            << "\t"
                            << it_cut->second
                            << "\t"
                            << lpm
                            << "\t"
                            << *it_medium
                            << "\t"
                            << *it_particle
                            << "\t"
                            << *it_dist
                            << "\t"
                            << (*it_energy).first
                            << "\t"
                            << (*it_energy).second
                            << "\t"
                            << final_energy
                            << std::endl;
                    }
                }
            }
        }
    }

    // Evering ok
    return 0;
}


int ProcColl_FinalEnergyParticleInteraction(std::string filename, std::string path_to_tables)
{
    std::ofstream out;
    out.open(filename, std::ofstream::out);

    if (!out.good())
    {
        printf("Cannot open file %s for writing\n", filename.c_str());
        return -1;
    }

    printf("Generate Test of ProcessCollection Final Energy Particle Interaction tables in %s\n", filename.c_str());

    // Header
    out << "ecut"
        << "\t"
        << "vcut"
        << "\t"
        << "lpm"
        << "\t"
        << "medium"
        << "\t"
        << "paricle"
        << "\t"
        << "energy"
        << "\t"
        << "particle_interaction"
        << "\t"
        << "rnd"
        << "\t"
        << "FinalEnergy"
        << std::endl;


    std::vector<double> energy;
    energy.push_back(10000);
    for (int i = 0; i < 10; i++)
    {
        energy.push_back(energy.at(i) * 10);
    }

    std::vector<std::string> medium;
    medium.push_back("ice"); medium.push_back("uranium"); medium.push_back("air");

    std::vector<std::string> particle;
    particle.push_back("mu"), particle.push_back("tau"), particle.push_back("e");

    std::vector<std::pair<int, double> > cuts;
    cuts.push_back(std::pair<int, double>(-1, -1));
    cuts.push_back(std::pair<int, double>(500, -1));
    cuts.push_back(std::pair<int, double>(-1, 0.05));
    cuts.push_back(std::pair<int, double>(500, 0.05));

    int lpm = 0;
    double rnd = 0.0;
    double final_energy = 0.0;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    out.precision(16);

    Vector3D particle_position(1.,1.,1.);
    Vector3D particle_direction(0,0,0);
    particle_direction.SetSphericalCoordinates(1, .20/180*PI, 20/180*PI);

    // cuts
    for (std::vector<std::pair<int, double> >::iterator it_cut = cuts.begin() ; it_cut != cuts.end(); ++it_cut)
    {
        EnergyCutSettings cuts(it_cut->first, it_cut->second);
        // medium
        for (std::vector<std::string>::iterator it_medium = medium.begin() ; it_medium != medium.end(); ++it_medium)
        {
            Medium* medium = MediumFactory::Get()->CreateMedium(*it_medium);

            // paricle
            for (std::vector<std::string>::iterator it_particle = particle.begin() ; it_particle != particle.end(); ++it_particle)
            {
                Particle particle(Particle::GetTypeFromName(*it_particle),particle_position,particle_direction,1e5,10);
                // Propagator prop(&medium, &cuts, particle.GetType(), "");
                ProcessCollection proc_col(&particle, medium, &cuts);
                proc_col.EnableInterpolation(path_to_tables);

                // energy
                for (std::vector<double>::iterator it_energy = energy.begin() ; it_energy != energy.end(); ++it_energy)
                {
                    // particle interaction
                    for (int particle_interaction = 1; particle_interaction > -1; --particle_interaction)
                    {
                        for (int i = 0; i < 5; ++i)
                        {
                            rnd = distribution(generator);
                            proc_col.GetParticle()->SetEnergy(*it_energy);
                            final_energy = proc_col.CalculateTrackingIntegal(*it_energy, rnd, particle_interaction);
                            final_energy = proc_col.CalculateFinalEnergy(*it_energy, rnd, particle_interaction);

                            out << it_cut->first
                                << "\t"
                                << it_cut->second
                                << "\t"
                                << lpm
                                << "\t"
                                << *it_medium
                                << "\t"
                                << *it_particle
                                << "\t"
                                << *it_energy
                                << "\t"
                                << particle_interaction
                                << "\t"
                                << rnd
                                << "\t"
                                << final_energy
                                << std::endl;
                        }
                    }
                }
            }
        }
    }

    // Evering ok
    return 0;
}


int ProcColl_Tracking(std::string filename, std::string path_to_tables)
{
    std::ofstream out;
    out.open(filename, std::ofstream::out);

    if (!out.good())
    {
        printf("Cannot open file %s for writing\n", filename.c_str());
        return -1;
    }

    printf("Generate Test of ProcessCollection Tracking tables in %s\n", filename.c_str());

    // Header
    out << "ecut"
        << "\t"
        << "vcut"
        << "\t"
        << "lpm"
        << "\t"
        << "medium"
        << "\t"
        << "paricle"
        << "\t"
        << "particle_interaction"
        << "\t"
        << "energy"
        << "\t"
        << "tracking"
        << std::endl;


    std::vector<double> energy;
    energy.push_back(10000);
    for (int i = 0; i < 7; i++)
    {
        energy.push_back(energy.at(i) * 10);
    }

    std::vector<std::string> medium;
    medium.push_back("ice"); medium.push_back("uranium"); medium.push_back("air");

    std::vector<std::string> particle;
    particle.push_back("mu"), particle.push_back("tau"), particle.push_back("e");

    std::vector<std::pair<int, double> > cuts;
    cuts.push_back(std::pair<int, double>(-1, -1));
    cuts.push_back(std::pair<int, double>(500, -1));
    cuts.push_back(std::pair<int, double>(-1, 0.05));
    cuts.push_back(std::pair<int, double>(500, 0.05));

    int lpm = 0;
    double rnd = 0.1;
    double tracking = 0;

    // std::default_random_engine generator;
    // std::uniform_real_distribution<double> distribution(0.0,1.0);

    out.precision(16);

    Vector3D particle_position(1.,1.,1.);
    Vector3D particle_direction(0,0,0);
    particle_direction.SetSphericalCoordinates(1, .20/180*PI, 20/180*PI);

    // cuts
    for (std::vector<std::pair<int, double> >::iterator it_cut = cuts.begin() ; it_cut != cuts.end(); ++it_cut)
    {
        EnergyCutSettings cuts(it_cut->first, it_cut->second);
        // medium
        for (std::vector<std::string>::iterator it_medium = medium.begin() ; it_medium != medium.end(); ++it_medium)
        {
            Medium* medium = MediumFactory::Get()->CreateMedium(*it_medium);

            // paricle
            for (std::vector<std::string>::iterator it_particle = particle.begin() ; it_particle != particle.end(); ++it_particle)
            {
                Particle particle(Particle::GetTypeFromName(*it_particle),particle_position,particle_direction,1e5,10);
                // Propagator prop(&medium, &cuts, particle.GetType(), "");
                ProcessCollection proc_col(&particle, medium, &cuts);
                proc_col.EnableInterpolation(path_to_tables);

                // energy
                for (std::vector<double>::iterator it_energy = energy.begin() ; it_energy != energy.end(); ++it_energy)
                {
                    // particle interaction
                    for (int particle_interaction = 1; particle_interaction > -1; --particle_interaction)
                    {
                        proc_col.GetParticle()->SetEnergy(*it_energy);
                        tracking = proc_col.CalculateTrackingIntegal(*it_energy, rnd, particle_interaction);

                        out << it_cut->first
                            << "\t"
                            << it_cut->second
                            << "\t"
                            << lpm
                            << "\t"
                            << *it_medium
                            << "\t"
                            << *it_particle
                            << "\t"
                            << particle_interaction
                            << "\t"
                            << *it_energy
                            << "\t"
                            << tracking
                            << std::endl;
                    }
                }
            }
        }
    }

    // Evering ok
    return 0;
}


int ProcColl_MakeDecay(std::string filename, std::string path_to_tables)
{
    std::ofstream out;
    out.open(filename, std::ofstream::out);

    if (!out.good())
    {
        printf("Cannot open file %s for writing\n", filename.c_str());
        return -1;
    }

    printf("Generate Test of ProcessCollection MakeDecay tables in %s\n", filename.c_str());

    // Header
    out << "ecut"
        << "\t"
        << "vcut"
        << "\t"
        << "lpm"
        << "\t"
        << "medium"
        << "\t"
        << "paricle"
        << "\t"
        << "rnd1"
        << "\t"
        << "rnd2"
        << "\t"
        << "rnd3"
        << "\t"
        << "energy"
        << "\t"
        << "product_energy"
        << "\t"
        << "out"
        << std::endl;


    std::vector<double> energy;
    energy.push_back(10000);
    for (int i = 0; i < 10; i++)
    {
        energy.push_back(energy.at(i) * 10);
    }

    std::vector<std::string> medium;
    medium.push_back("ice"); medium.push_back("uranium"); medium.push_back("air");

    std::vector<std::string> particle;
    particle.push_back("mu"), particle.push_back("tau"), particle.push_back("e");

    std::vector<std::pair<int, double> > cuts;
    cuts.push_back(std::pair<int, double>(-1, -1));
    cuts.push_back(std::pair<int, double>(500, -1));
    cuts.push_back(std::pair<int, double>(-1, 0.05));
    cuts.push_back(std::pair<int, double>(500, 0.05));

    int lpm = 0;

    double rnd1 = 0.0;
    double rnd2 = 0.0;
    double rnd3 = 0.0;

    std::pair<double, ParticleType::Enum> decay_out;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    out.precision(16);

    Vector3D particle_position(1.,1.,1.);
    Vector3D particle_direction(0,0,0);
    particle_direction.SetSphericalCoordinates(1, .20/180*PI, 20/180*PI);

    // cuts
    for (std::vector<std::pair<int, double> >::iterator it_cut = cuts.begin() ; it_cut != cuts.end(); ++it_cut)
    {
        EnergyCutSettings cuts(it_cut->first, it_cut->second);
        // medium
        for (std::vector<std::string>::iterator it_medium = medium.begin() ; it_medium != medium.end(); ++it_medium)
        {
            Medium* medium = MediumFactory::Get()->CreateMedium(*it_medium);

            // paricle
            for (std::vector<std::string>::iterator it_particle = particle.begin() ; it_particle != particle.end(); ++it_particle)
            {
                Particle particle(Particle::GetTypeFromName(*it_particle),particle_position,particle_direction,1e5,10);
                // Propagator prop(&medium, &cuts, particle.GetType(), "");
                ProcessCollection proc_col(&particle, medium, &cuts);
                proc_col.EnableInterpolation(path_to_tables);

                // energy
                for (std::vector<double>::iterator it_energy = energy.begin() ; it_energy != energy.end(); ++it_energy)
                {
                    // particle interaction
                    for (int i = 0; i < 5; ++i)
                    {
                        rnd1 = distribution(generator);
                        rnd2 = distribution(generator);
                        rnd3 = distribution(generator);

                        proc_col.GetParticle()->SetEnergy(*it_energy);
                        decay_out = proc_col.MakeDecay(rnd1,rnd2,rnd3);

                        out << it_cut->first
                            << "\t"
                            << it_cut->second
                            << "\t"
                            << lpm
                            << "\t"
                            << *it_medium
                            << "\t"
                            << *it_particle
                            << "\t"
                            << rnd1
                            << "\t"
                            << rnd2
                            << "\t"
                            << rnd3
                            << "\t"
                            << *it_energy
                            << "\t"
                            << decay_out.first
                            << "\t"
                            << Particle::GetName(decay_out.second)
                            << std::endl;
                    }
                }
            }
        }
    }

    // Evering ok
    return 0;
}


void usage(const char *argv[])
{
    std::cout << "Usage: "
        << argv[0] << " -t <TEST> [-h | -f <PATH_TO_TABLES> | -s <PATH_TO_SAVE_FILES>] \n"
        << "\t"
        << "-h\thelp\n"
        << "\t"
        << "-f\tPath to interpolation tables\n"
        << "\t"
        << "-s\tPath where test tables should be stored\n"
        << "\t"
        << "-t\t<TEST> specifies which test tables should be created.\n\n"
        << "\t"
        << "Choose <TEST> from:\n"
        << "\t\t brems\t Bremsstrahlung\n"
        << "\t\t epair\t Epair Production\n"
        << "\t\t ioniz\t Ionization\n"
        << "\t\t photo\t Photonuclear\n"
        << "\t\t pcol\t Process Collection\n"
        << "\t\t cont\t Continous Randomization\n"
        << "\t\t rnd\t Random Numbers\n"
        << "\t\t all\t All\n\n";
}

// ------------------------------------------------------------------------- //
// Main
// ------------------------------------------------------------------------- //

int main(int argc, const char *argv[])
{

    InputParser input(argc, argv);

    // Help
    if (input.cmdOptionExists("-h"))
    {
        usage(argv);
        exit(0);
    }

    // Path to tables
    std::string path_to_tables = input.getCmdOption("-f");

    if (path_to_tables.empty())
    {
        // std::cout << "Warning: No path given! Path to tables is set to " << "\"" << path_to_tables << "\"" << std::endl;
        std::cout << "Warning: No path given! Interplolation tables will be stored in memory." << std::endl;
    }
    else
    {
        if (is_dir(path_to_tables) == false)
        {
            std::cout << "Warning: Given path \""<< path_to_tables << "\" is not a directory! Interpolation tables will be stored in memory." << std::endl;
            path_to_tables = "";
        }
    }

    // Path to save tables
    std::string path_to_save = input.getCmdOption("-s");

    if (path_to_save.empty())
    {
        path_to_save = "";
        std::cout << "Warning: No destination path given! Test tables will be save to " << "\"" << "." << "\"" << std::endl;
    }
    else
    {
        std::string::iterator it = path_to_save.end();
        if (strcmp(&*it, "/") != 0)
        {
            path_to_save.append("/");
        }

        if (is_dir(path_to_save) == false)
        {
            std::cout << "Warning: Given destination path \""<< path_to_save << "\" is not a directory! Test tables will be saved in \".\" ." << std::endl;
            path_to_save = "";
        }
        else
        {
            std::cout << "Test tables will be save to " << "\"" << path_to_save << "\"" << std::endl;
        }
    }


    // Test decision
    std::string opt_brems = "brems";
    std::string opt_epair = "epair";
    std::string opt_ioniz = "ioniz";
    std::string opt_photo = "photo";
    std::string opt_pcol = "pcol";
    std::string opt_cont = "cont";
    std::string opt_rnd = "rnd";
    std::string opt_all = "all";

    std::vector<std::string> tests;
    tests.push_back(opt_brems);
    tests.push_back(opt_ioniz);
    tests.push_back(opt_epair);
    tests.push_back(opt_photo);
    tests.push_back(opt_pcol);
    tests.push_back(opt_cont);
    tests.push_back(opt_rnd);
    tests.push_back(opt_all);

    std::string chosen_test = input.getCmdOption("-t");
    if (chosen_test.empty())
    {
        std::cout << "Warning: No test given! Abort! " << std::endl;
        usage(argv);
        exit(0);
    }
    else if (std::find(tests.begin(), tests.end(), chosen_test) == tests.end())
    {
        std::cout << "Warning: " << chosen_test << "is no valid test argument! Abort! " << std::endl;
        usage(argv);
        exit(0);
    }
    else
    {
        std::cout << "Create tables for " << chosen_test << std::endl;
        std::cout << "==================================" << std::endl;
    }

    // File with random numbers. Needed for many tests.
    std::string path_to_rnd = path_to_save + "rnd.txt";

    // --------------------------------------------------------------------- //
    // Tests
    // --------------------------------------------------------------------- //

    if ((chosen_test == opt_rnd) || (chosen_test == opt_all))
    {
        Create_RandomNumbers(path_to_rnd);
    }

    if ((chosen_test == opt_brems) || (chosen_test == opt_all))
    {
        std::cout << "Create tables for Bremsstrahlung" << std::endl;
        std::cout << "================================" << std::endl;

        if (is_file(path_to_rnd) == false)
        {
            Create_RandomNumbers(path_to_rnd);
        }

        Bremsstrahlung_Test_of_dEdx(path_to_save + "Brems_dEdx.txt", "dEdx", "", "");
        Bremsstrahlung_Test_of_dNdx(path_to_save + "Brems_dNdx.txt", "dNdx", "", "");
        Bremsstrahlung_Test_of_dNdxrnd(path_to_save + "Brems_dNdxrnd.txt", "dNdxRnd", "", path_to_rnd);
        Bremsstrahlung_Test_of_e(path_to_save + "Brems_e.txt", "StochasticLoss", "", path_to_rnd);
        Bremsstrahlung_Test_of_dEdx_Interpolant(path_to_save + "Brems_dEdx_interpol.txt", "dEdx", path_to_tables, "");
        Bremsstrahlung_Test_of_dNdx_Interpolant(path_to_save + "Brems_dNdx_interpol.txt", "dNdx", path_to_tables, "");
        Bremsstrahlung_Test_of_dNdxrnd_Interpolant(path_to_save + "Bremsstrahlung_dNdxrnd_interpol.txt", "dNdxRnd", path_to_tables, path_to_rnd);
        Bremsstrahlung_Test_of_e_Interpolant(path_to_save + "Brems_e_interpol.txt", "StochasticLoss", path_to_tables, path_to_rnd);
    }

    if ((chosen_test == opt_epair) || (chosen_test == opt_all))
    {
        std::cout << "Create tables for Epair Production" << std::endl;
        std::cout << "==================================" << std::endl;

        if (is_file(path_to_rnd) == false)
        {
            std::cout << path_to_rnd << std::endl;
            Create_RandomNumbers(path_to_rnd);
        }
        std::cout << "out " << path_to_rnd << std::endl;

        Epairproduction_Test_of_dEdx(path_to_save + "Epair_dEdx.txt", "dEdx", "", "");
        Epairproduction_Test_of_dNdx(path_to_save + "Epair_dNdx.txt", "dNdx", "", "");
        Epairproduction_Test_of_dNdxrnd(path_to_save + "Epair_dNdxrnd.txt", "dNdxRnd", "", path_to_rnd);
        Epairproduction_Test_of_e(path_to_save + "Epair_e.txt", "StochasticLoss", "", path_to_rnd);
        Epairproduction_Test_of_dEdx_Interpolant(path_to_save + "Epair_dEdx_interpol.txt", "dEdx", path_to_tables, "");
        Epairproduction_Test_of_dNdx_Interpolant(path_to_save + "Epair_dNdx_interpol.txt", "dNdx", path_to_tables, "");
        Epairproduction_Test_of_dNdxrnd_Interpolant(path_to_save + "Epair_dNdxrnd_interpol.txt", "dNdxRnd", path_to_tables, path_to_rnd);
        Epairproduction_Test_of_e_Interpolant(path_to_save + "Epair_e_interpol.txt", "StochasticLoss", path_to_tables, path_to_rnd);
    }

    if ((chosen_test == opt_ioniz) || (chosen_test == opt_all))
    {
        std::cout << "Create tables for Ionization" << std::endl;
        std::cout << "============================" << std::endl;

        if (is_file(path_to_rnd) == false)
        {
            Create_RandomNumbers(path_to_rnd);
        }

        Ionization_Test_of_dEdx(path_to_save + "Ioniz_dEdx.txt", "dEdx", "", "");
        Ionization_Test_of_dNdx(path_to_save + "Ioniz_dNdx.txt", "dNdx", "", "");
        Ionization_Test_of_dNdxrnd(path_to_save + "Ioniz_dNdxrnd.txt", "dNdxRnd", "", path_to_rnd);
        Ionization_Test_of_e(path_to_save + "Ioniz_e.txt", "StochasticLoss", "", path_to_rnd);
        Ionization_Test_of_dEdx_Interpolant(path_to_save + "Ioniz_dEdx_interpol.txt", "dEdx", path_to_tables, "");
        Ionization_Test_of_dNdx_Interpolant(path_to_save + "Ioniz_dNdx_interpol.txt", "dNdx", path_to_tables, "");
        Ionization_Test_of_dNdxrnd_Interpolant(path_to_save + "Ionization_dNdxrnd_interpol.txt", "dNdxRnd", path_to_tables, path_to_rnd);
        Ionization_Test_of_e_Interpolant(path_to_save + "Ioniz_e_interpol.txt", "StochasticLoss", path_to_tables, path_to_rnd);
    }

    if ((chosen_test == opt_photo) || (chosen_test == opt_all))
    {
        std::cout << "Create tables for Photonuclear" << std::endl;
        std::cout << "==============================" << std::endl;

        if (is_file(path_to_rnd) == false)
        {
            Create_RandomNumbers(path_to_rnd);
        }

        Photonuclear_Test_of_dEdx(path_to_save + "Photo_dEdx.txt", "dEdx", "", "");
        Photonuclear_Test_of_dNdx(path_to_save + "Photo_dNdx.txt", "dNdx", "", "");
        Photonuclear_Test_of_dNdxrnd(path_to_save + "Photo_dNdxrnd.txt", "dNdxRnd", "", path_to_rnd);
        Photonuclear_Test_of_e(path_to_save + "Photo_e.txt", "StochasticLoss", "", path_to_rnd);
        Photonuclear_Test_of_dEdx_Interpolant(path_to_save + "Photo_dEdx_interpol.txt", "dEdx", path_to_tables, "");
        Photonuclear_Test_of_dNdx_Interpolant(path_to_save + "Photo_dNdx_interpol.txt", "dNdx", path_to_tables, "");
        Photonuclear_Test_of_dNdxrnd_Interpolant(path_to_save + "Photo_dNdxrnd_interpol.txt", "dNdxRnd", path_to_tables, path_to_rnd);
        Photonuclear_Test_of_e_Interpolant(path_to_save + "Photo_e_interpol.txt", "StochasticLoss", path_to_tables, path_to_rnd);
    }

    if ((chosen_test == opt_cont) || (chosen_test == opt_all))
    {
        std::cout << "Create tables for Continous Randomization" << std::endl;
        std::cout << "=========================================" << std::endl;

        std::string a = path_to_save + "ContinousRandomization_interpol.txt";
        setting_loop_contrand(a);
    }

    if ((chosen_test == opt_pcol) || (chosen_test == opt_all))
    {
        std::cout << "Create tables for Process Collection" << std::endl;
        std::cout << "====================================" << std::endl;

        ProcColl_Stochastics(path_to_save + "ProcColl_Stoch.txt", path_to_tables);
        ProcColl_Displacement(path_to_save + "ProcColl_Disp.txt");
        ProcColl_FinalEnergyDist(path_to_save + "ProcColl_FinalEnergyDist.txt", path_to_tables);
        ProcColl_FinalEnergyParticleInteraction(path_to_save + "ProcColl_FinalEnergyParticleInteraction.txt", path_to_tables);
        ProcColl_Tracking(path_to_save + "ProcColl_Tracking.txt", path_to_tables);
        ProcColl_MakeDecay(path_to_save + "ProcColl_MakeDecay.txt", path_to_tables);
    }

    printf("Done...\n");
    return 0;
}
