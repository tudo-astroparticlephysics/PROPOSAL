
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/CrossSections.h"
#include <iostream>
#include <string>

using namespace PROPOSAL;

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

int ConvertOldToNewBremsstrahlungParametrization(int old_param)
{
    switch (old_param)
    {
        case 1:
            return ParametrizationType::BremsKelnerKokoulinPetrukhin;
            break;
        case 2:
            return ParametrizationType::BremsAndreevBezrukovBugaev;
            break;
        case 3:
            return ParametrizationType::BremsPetrukhinShestakov;
            break;
        case 4:
            return ParametrizationType::BremsCompleteScreeningCase;
            break;
        default:
            printf("Wrong Bremsstrahlung parametrization type '%i'. Default to BremsKelnerKokoulinPetrukhin."
                , old_param);
            return ParametrizationType::BremsKelnerKokoulinPetrukhin;
    }
}

// ------------------------------------------------------------------------- //
int Bremsstrahlung_Test_of_dEdx(std::string filename)
{
    std::ofstream out;
    out.open(filename, std::ofstream::out);

    if (!out.good())
    {
        printf("Cannot open file %s for writing\n", filename.c_str());
        return -1;
    }

    // Header
    out << "para" << "\t" << "ecut" << "\t" << "vcut" << "\t" << "lpm" << "\t" << "energy" << "\t" << "medium" << "\t" << "particle" << "\t" << "dEdx" << std::endl;

    std::vector<int> param;
    param.push_back(1);
    param.push_back(2);
    param.push_back(3);
    param.push_back(4);

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

    double dEdx_new;

    out.precision(16);

    // cuts
    for (std::vector<std::pair<int, double> >::iterator it_cut = cuts.begin() ; it_cut != cuts.end(); ++it_cut)
    {
        EnergyCutSettings cuts(it_cut->first, it_cut->second);

        // medium
        for (std::vector<std::string>::iterator it_medium = medium.begin() ; it_medium != medium.end(); ++it_medium)
        {
            Medium medium(Medium::GetTypeFromName(*it_medium), 1.);

            // paricle
            for (std::vector<std::string>::iterator it_particle = particle.begin() ; it_particle != particle.end(); ++it_particle)
            {
                PROPOSALParticle particle(PROPOSALParticle::GetTypeFromName(*it_particle),1.,1.,1,.20,20,1e5,10);
                Bremsstrahlung brems(&particle, &medium, &cuts);

                // params
                for (std::vector<int>::iterator it_param = param.begin() ; it_param != param.end(); ++it_param)
                {
                    brems.SetParametrization(static_cast<ParametrizationType::Enum>(ConvertOldToNewBremsstrahlungParametrization(*it_param)));

                    // lpm
                    for (int lpm = 0; lpm < 2; ++lpm)
                    {
                        brems.EnableLpmEffect(lpm);

                        // energy
                        for (std::vector<double>::iterator it_energy = energy.begin() ; it_energy != energy.end(); ++it_energy)
                        {
                            particle.SetEnergy(*it_energy);

                            dEdx_new=brems.CalculatedEdx();

                            out << *it_param << "\t" << it_cut->first << "\t" << it_cut->second << "\t" << lpm << "\t" << *it_energy << "\t" << *it_medium << "\t" << *it_particle << "\t" << dEdx_new << std::endl;
                        }
                    }
                }
            }
        }
    }

    // Evering ok
    return 0;
}


// ------------------------------------------------------------------------- //
int Bremsstrahlung_Test_of_dNdx(std::string filename)
{
    std::ofstream out;
    out.open(filename, std::ofstream::out);

    if (!out.good())
    {
        printf("Cannot open file %s for writing\n", filename.c_str());
        return -1;
    }

    // Header
    out << "para" << "\t" << "ecut" << "\t" << "vcut" << "\t" << "lpm" << "\t" << "energy" << "\t" << "medium" << "\t" << "particle" << "\t" << "dNdx" << std::endl;

    std::vector<int> param;
    param.push_back(1);
    param.push_back(2);
    param.push_back(3);
    param.push_back(4);

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

    double dNdx_new;

    out.precision(16);

    // cuts
    for (std::vector<std::pair<int, double> >::iterator it_cut = cuts.begin() ; it_cut != cuts.end(); ++it_cut)
    {
        EnergyCutSettings cuts(it_cut->first, it_cut->second);

        // medium
        for (std::vector<std::string>::iterator it_medium = medium.begin() ; it_medium != medium.end(); ++it_medium)
        {
            Medium medium(Medium::GetTypeFromName(*it_medium), 1.);

            // paricle
            for (std::vector<std::string>::iterator it_particle = particle.begin() ; it_particle != particle.end(); ++it_particle)
            {
                PROPOSALParticle particle(PROPOSALParticle::GetTypeFromName(*it_particle),1.,1.,1,.20,20,1e5,10);
                Bremsstrahlung brems(&particle, &medium, &cuts);

                // params
                for (std::vector<int>::iterator it_param = param.begin() ; it_param != param.end(); ++it_param)
                {
                    brems.SetParametrization(static_cast<ParametrizationType::Enum>(ConvertOldToNewBremsstrahlungParametrization(*it_param)));

                    // lpm
                    for (int lpm = 0; lpm < 2; ++lpm)
                    {
                        brems.EnableLpmEffect(lpm);

                        // energy
                        for (std::vector<double>::iterator it_energy = energy.begin() ; it_energy != energy.end(); ++it_energy)
                        {
                            particle.SetEnergy(*it_energy);

                            dNdx_new=brems.CalculatedNdx();

                            out << *it_param << "\t" << it_cut->first << "\t" << it_cut->second << "\t" << lpm << "\t" << *it_energy << "\t" << *it_medium << "\t" << *it_particle << "\t" << dNdx_new << std::endl;
                        }
                    }
                }
            }
        }
    }

    // Evering ok
    return 0;
}

// ------------------------------------------------------------------------- //
int Bremsstrahlung_Test_of_dNdxrnd(std::string filename)
{
    std::ofstream out;
    out.open(filename, std::ofstream::out);

    if (!out.good())
    {
        printf("Cannot open file %s for writing\n", filename.c_str());
        return -1;
    }

    // Header
    out << "para" << "\t" << "ecut" << "\t" << "vcut" << "\t" << "lpm" << "\t" << "energy" << "\t" << "medium" << "\t" << "particle" << "\t" << "dNdxrnd" << std::endl;

    RndFromFile Rand("bin/TestFiles/rnd.txt");

    std::vector<int> param;
    param.push_back(1);
    param.push_back(2);
    param.push_back(3);
    param.push_back(4);

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

    double dNdx_new;

    out.precision(16);

    // cuts
    for (std::vector<std::pair<int, double> >::iterator it_cut = cuts.begin() ; it_cut != cuts.end(); ++it_cut)
    {
        EnergyCutSettings cuts(it_cut->first, it_cut->second);

        // medium
        for (std::vector<std::string>::iterator it_medium = medium.begin() ; it_medium != medium.end(); ++it_medium)
        {
            Medium medium(Medium::GetTypeFromName(*it_medium), 1.);

            // paricle
            for (std::vector<std::string>::iterator it_particle = particle.begin() ; it_particle != particle.end(); ++it_particle)
            {
                PROPOSALParticle particle(PROPOSALParticle::GetTypeFromName(*it_particle),1.,1.,1,.20,20,1e5,10);
                Bremsstrahlung brems(&particle, &medium, &cuts);

                // params
                for (std::vector<int>::iterator it_param = param.begin() ; it_param != param.end(); ++it_param)
                {
                    brems.SetParametrization(static_cast<ParametrizationType::Enum>(ConvertOldToNewBremsstrahlungParametrization(*it_param)));

                    // lpm
                    for (int lpm = 0; lpm < 2; ++lpm)
                    {
                        brems.EnableLpmEffect(lpm);

                        // energy
                        for (std::vector<double>::iterator it_energy = energy.begin() ; it_energy != energy.end(); ++it_energy)
                        {
                            particle.SetEnergy(*it_energy);

                            dNdx_new=brems.CalculatedNdx(Rand.rnd());

                            out << *it_param << "\t" << it_cut->first << "\t" << it_cut->second << "\t" << lpm << "\t" << *it_energy << "\t" << *it_medium << "\t" << *it_particle << "\t" << dNdx_new << std::endl;
                        }
                    }
                }
            }
        }
    }

    // Evering ok
    return 0;
}

// ------------------------------------------------------------------------- //
int Bremsstrahlung_Test_of_e(std::string filename)
{
    std::ofstream out;
    out.open(filename, std::ofstream::out);

    if (!out.good())
    {
        printf("Cannot open file %s for writing\n", filename.c_str());
        return -1;
    }

    // Header
    out << "para" << "\t" << "ecut" << "\t" << "vcut" << "\t" << "lpm" << "\t" << "energy" << "\t" << "medium" << "\t" << "particle" << "\t" << "dNdxrnd" << std::endl;

    RndFromFile Rand("bin/TestFiles/rnd.txt");
    RndFromFile Rand2("bin/TestFiles/rnd.txt");

    std::vector<int> param;
    param.push_back(1);
    param.push_back(2);
    param.push_back(3);
    param.push_back(4);

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

    double e_new;
    double rnd1;
    double rnd2;

    out.precision(16);

    // cuts
    for (std::vector<std::pair<int, double> >::iterator it_cut = cuts.begin() ; it_cut != cuts.end(); ++it_cut)
    {
        EnergyCutSettings cuts(it_cut->first, it_cut->second);

        // medium
        for (std::vector<std::string>::iterator it_medium = medium.begin() ; it_medium != medium.end(); ++it_medium)
        {
            Medium medium(Medium::GetTypeFromName(*it_medium), 1.);

            // paricle
            for (std::vector<std::string>::iterator it_particle = particle.begin() ; it_particle != particle.end(); ++it_particle)
            {
                PROPOSALParticle particle(PROPOSALParticle::GetTypeFromName(*it_particle),1.,1.,1,.20,20,1e5,10);
                Bremsstrahlung brems(&particle, &medium, &cuts);

                // params
                for (std::vector<int>::iterator it_param = param.begin() ; it_param != param.end(); ++it_param)
                {
                    brems.SetParametrization(static_cast<ParametrizationType::Enum>(ConvertOldToNewBremsstrahlungParametrization(*it_param)));

                    // lpm
                    for (int lpm = 0; lpm < 2; ++lpm)
                    {
                        brems.EnableLpmEffect(lpm);

                        // energy
                        for (std::vector<double>::iterator it_energy = energy.begin() ; it_energy != energy.end(); ++it_energy)
                        {
                            particle.SetEnergy(*it_energy);

                            rnd1 = Rand.rnd();
                            rnd2 = Rand2.rnd();

                            e_new=brems.CalculateStochasticLoss(rnd1, rnd2);

                            out << *it_param << "\t" << it_cut->first << "\t" << it_cut->second << "\t" << lpm << "\t" << *it_energy << "\t" << *it_medium << "\t" << *it_particle << "\t" << e_new << std::endl;
                        }
                    }
                }
            }
        }
    }

    // Evering ok
    return 0;
}

// ------------------------------------------------------------------------- //
int Bremsstrahlung_Test_of_dEdx_Interpolant(std::string filename)
{
    std::ofstream out;
    out.open(filename, std::ofstream::out);

    if (!out.good())
    {
        printf("Cannot open file %s for writing\n", filename.c_str());
        return -1;
    }

    // Header
    out << "para" << "\t" << "ecut" << "\t" << "vcut" << "\t" << "lpm" << "\t" << "energy" << "\t" << "medium" << "\t" << "particle" << "\t" << "dEdx" << std::endl;

    std::vector<int> param;
    param.push_back(1);
    param.push_back(2);
    param.push_back(3);
    param.push_back(4);

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

    double dEdx_new;

    out.precision(16);

    // cuts
    for (std::vector<std::pair<int, double> >::iterator it_cut = cuts.begin() ; it_cut != cuts.end(); ++it_cut)
    {
        EnergyCutSettings cuts(it_cut->first, it_cut->second);

        // medium
        for (std::vector<std::string>::iterator it_medium = medium.begin() ; it_medium != medium.end(); ++it_medium)
        {
            Medium medium(Medium::GetTypeFromName(*it_medium), 1.);

            // paricle
            for (std::vector<std::string>::iterator it_particle = particle.begin() ; it_particle != particle.end(); ++it_particle)
            {
                PROPOSALParticle particle(PROPOSALParticle::GetTypeFromName(*it_particle),1.,1.,1,.20,20,1e5,10);
                Bremsstrahlung brems(&particle, &medium, &cuts);

                // params
                for (std::vector<int>::iterator it_param = param.begin() ; it_param != param.end(); ++it_param)
                {
                    brems.SetParametrization(static_cast<ParametrizationType::Enum>(ConvertOldToNewBremsstrahlungParametrization(*it_param)));

                    // lpm
                    for (int lpm = 0; lpm < 2; ++lpm)
                    {
                        brems.EnableLpmEffect(lpm);
                        brems.EnableDEdxInterpolation();

                        // energy
                        for (std::vector<double>::iterator it_energy = energy.begin() ; it_energy != energy.end(); ++it_energy)
                        {
                            particle.SetEnergy(*it_energy);

                            dEdx_new=brems.CalculatedEdx();

                            out << *it_param << "\t" << it_cut->first << "\t" << it_cut->second << "\t" << lpm << "\t" << *it_energy << "\t" << *it_medium << "\t" << *it_particle << "\t" << dEdx_new << std::endl;
                        }
                    }
                }
            }
        }
    }

    // Evering ok
    return 0;
}

// ------------------------------------------------------------------------- //
int Bremsstrahlung_Test_of_dNdx_Interpolant(std::string filename)
{
    std::ofstream out;
    out.open(filename, std::ofstream::out);

    if (!out.good())
    {
        printf("Cannot open file %s for writing\n", filename.c_str());
        return -1;
    }

    // Header
    out << "para" << "\t" << "ecut" << "\t" << "vcut" << "\t" << "lpm" << "\t" << "energy" << "\t" << "medium" << "\t" << "particle" << "\t" << "dNdx" << std::endl;

    std::vector<int> param;
    param.push_back(1);
    param.push_back(2);
    param.push_back(3);
    param.push_back(4);

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

    double dNdx_new;

    out.precision(16);

    // cuts
    for (std::vector<std::pair<int, double> >::iterator it_cut = cuts.begin() ; it_cut != cuts.end(); ++it_cut)
    {
        EnergyCutSettings cuts(it_cut->first, it_cut->second);

        // medium
        for (std::vector<std::string>::iterator it_medium = medium.begin() ; it_medium != medium.end(); ++it_medium)
        {
            Medium medium(Medium::GetTypeFromName(*it_medium), 1.);

            // paricle
            for (std::vector<std::string>::iterator it_particle = particle.begin() ; it_particle != particle.end(); ++it_particle)
            {
                PROPOSALParticle particle(PROPOSALParticle::GetTypeFromName(*it_particle),1.,1.,1,.20,20,1e5,10);
                Bremsstrahlung brems(&particle, &medium, &cuts);

                // params
                for (std::vector<int>::iterator it_param = param.begin() ; it_param != param.end(); ++it_param)
                {
                    brems.SetParametrization(static_cast<ParametrizationType::Enum>(ConvertOldToNewBremsstrahlungParametrization(*it_param)));

                    // lpm
                    for (int lpm = 0; lpm < 2; ++lpm)
                    {
                        brems.EnableLpmEffect(lpm);
                        brems.EnableDNdxInterpolation();

                        // energy
                        for (std::vector<double>::iterator it_energy = energy.begin() ; it_energy != energy.end(); ++it_energy)
                        {
                            particle.SetEnergy(*it_energy);

                            dNdx_new=brems.CalculatedNdx();

                            out << *it_param << "\t" << it_cut->first << "\t" << it_cut->second << "\t" << lpm << "\t" << *it_energy << "\t" << *it_medium << "\t" << *it_particle << "\t" << dNdx_new << std::endl;
                        }
                    }
                }
            }
        }
    }

    // Evering ok
    return 0;
}

// ------------------------------------------------------------------------- //
int Bremsstrahlung_Test_of_e_Interpolant(std::string filename)
{
    std::ofstream out;
    out.open(filename, std::ofstream::out);

    if (!out.good())
    {
        printf("Cannot open file %s for writing\n", filename.c_str());
        return -1;
    }

    // Header
    out << "para" << "\t" << "ecut" << "\t" << "vcut" << "\t" << "lpm" << "\t" << "energy" << "\t" << "medium" << "\t" << "particle" << "\t" << "dNdxrnd" << std::endl;

    RndFromFile Rand("bin/TestFiles/rnd.txt");
    RndFromFile Rand2("bin/TestFiles/rnd.txt");

    std::vector<int> param;
    param.push_back(1);
    param.push_back(2);
    param.push_back(3);
    param.push_back(4);

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

    double e_new;
    double rnd1;
    double rnd2;

    out.precision(16);

    // cuts
    for (std::vector<std::pair<int, double> >::iterator it_cut = cuts.begin() ; it_cut != cuts.end(); ++it_cut)
    {
        EnergyCutSettings cuts(it_cut->first, it_cut->second);

        // medium
        for (std::vector<std::string>::iterator it_medium = medium.begin() ; it_medium != medium.end(); ++it_medium)
        {
            Medium medium(Medium::GetTypeFromName(*it_medium), 1.);

            // paricle
            for (std::vector<std::string>::iterator it_particle = particle.begin() ; it_particle != particle.end(); ++it_particle)
            {
                PROPOSALParticle particle(PROPOSALParticle::GetTypeFromName(*it_particle),1.,1.,1,.20,20,1e5,10);
                Bremsstrahlung brems(&particle, &medium, &cuts);

                // params
                for (std::vector<int>::iterator it_param = param.begin() ; it_param != param.end(); ++it_param)
                {
                    brems.SetParametrization(static_cast<ParametrizationType::Enum>(ConvertOldToNewBremsstrahlungParametrization(*it_param)));

                    // lpm
                    for (int lpm = 0; lpm < 2; ++lpm)
                    {
                        brems.EnableLpmEffect(lpm);
                        brems.EnableDNdxInterpolation();

                        // energy
                        for (std::vector<double>::iterator it_energy = energy.begin() ; it_energy != energy.end(); ++it_energy)
                        {
                            particle.SetEnergy(*it_energy);

                            rnd1 = Rand.rnd();
                            rnd2 = Rand2.rnd();

                            e_new=brems.CalculateStochasticLoss(rnd1, rnd2);

                            out << *it_param << "\t" << it_cut->first << "\t" << it_cut->second << "\t" << lpm << "\t" << *it_energy << "\t" << *it_medium << "\t" << *it_particle << "\t" << e_new << std::endl;
                        }
                    }
                }
            }
        }
    }

    // Evering ok
    return 0;
}

// ------------------------------------------------------------------------- //
// Main
// ------------------------------------------------------------------------- //

int main(int argc, const char *argv[])
{
    Bremsstrahlung_Test_of_dEdx("Brems_dEdx_test.txt");
    Bremsstrahlung_Test_of_dNdx("Brems_dNdx_test.txt");
    Bremsstrahlung_Test_of_dNdxrnd("Brems_Test_of_dNdxrnd.txt");
    Bremsstrahlung_Test_of_e("Brems_e.txt");
    Bremsstrahlung_Test_of_dEdx_Interpolant("Brems_dEdx_interpol.txt");
    Bremsstrahlung_Test_of_dNdx_Interpolant("Brems_dNdx_interpol.txt");
    Bremsstrahlung_Test_of_e_Interpolant("Brems_e_interpol.txt");

    return 0;
}
