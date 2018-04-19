
#include <fstream>
#include <iostream>

#include "PROPOSAL/PROPOSAL.h"

using namespace PROPOSAL;

int main(int argc, const char* argv[])
{
    int statistics = 1;

    if (argc >= 2)
    {
        int result;
        std::istringstream ss(argv[1]);
        ss.imbue(std::locale::classic());
        ss >> result;

        statistics = result;
        std::cout << "propagate " << statistics << " particles" << std::endl;
    }

    if (argc >= 3)
    {
        std::cerr << "reading random generator seed from" << argv[2] << std::endl;
        std::ifstream myfile(argv[2]);
        PROPOSAL::RandomGenerator::Get().Deserialize(myfile);
        myfile.close();
    }

    RandomGenerator::Get().SetSeed(1234);

    // ----[ Sector ]---------------------------------------- //

    Sector::Definition sec_def;

    sec_def.utility_def.epair_def.lpm_effect = true;
    sec_def.utility_def.brems_def.lpm_effect = true;

    sec_def.location = Sector::ParticleLocation::InsideDetector;

    sec_def.do_continuous_randomization = true;
    sec_def.do_exact_time_calculation   = true;

    sec_def.scattering_model = ScatteringFactory::HighlandIntegral;

    sec_def.cut_settings.SetEcut(500);
    sec_def.cut_settings.SetVcut(-1);

    sec_def.SetGeometry(Sphere(Vector3D(), 1e18, 0));
    sec_def.SetMedium(Ice());

    std::vector<Sector::Definition> sector_defintions;
    sector_defintions.push_back(sec_def);

    Propagator prop(MuMinusDef::Get(), sector_defintions, Sphere(Vector3D(), 1e18, 0));

    Particle& particle_mu = prop.GetParticle();

    particle_mu.SetEnergy(1e9); // [MeV]
    particle_mu.SetPropagatedDistance(0);
    particle_mu.SetPosition(Vector3D(0, 0, 0));
    particle_mu.SetDirection(Vector3D(0, 0, -1));

    Particle backup_mu = Particle(particle_mu);

    // std::ofstream myfile;
    // std::cout << '\n';
    for (int i = 0; i < statistics; ++i)
    {
        // std::cerr << "loop: " << i << std::endl;
        // std::cout << "RandomGenerator state:\n";
        // PROPOSAL::RandomGenerator::Get().Serialize(std::cout);
        // std::cout << '\n';

        // PROPOSAL::RandomGenerator::Get().Serialize(std::cout);
        // std::cout << '\n' << std::endl;

        particle_mu.InjectState(backup_mu);

        prop.Propagate(100);

        // myfile.open ("example.txt", std::ios::app);
        // myfile << "Writing this to a file.\n";
        // myfile.close();
    }
    std::cout << "concentration on difficult energy" << std::endl;

    double energy = 146.768; // 0.274374 + 4.63716e-10
    // double energy = 142.58; // 0.287524 - 1.01164e-8
    ParticleDef particle_def = MuMinusDef::Get();
    Ice medium;
    EnergyCutSettings ecuts;
    double multiplier = 1.;
    Ionization Ioniz(particle_def, medium, ecuts, multiplier);
    IonizIntegral Ioniz_Int(Ioniz);
    double dEdx_new = Ioniz_Int.CalculatedEdx(energy);

    Output::getInstance().ClearSecondaryVector();
}
