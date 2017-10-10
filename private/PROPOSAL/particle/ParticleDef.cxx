/*! \file   ParticleDef.cxx
*   \brief  Source file for definition of the ParticleDef class object.
*
*   For more details see the class documentation.
*
*   \date   Sat Aug  5 14:47:16 CEST 2017
*   \author Mario Dunsch
*/

#include <string>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/decay/LeptonicDecayChannel.h"
#include "PROPOSAL/decay/TwoBodyPhaseSpace.h"
#include "PROPOSAL/decay/StableChannel.h"

#include "PROPOSAL/methods.h"

#define PARTICLE_IMP(cls, MASS, LIFETIME, CHARGE)                                                                      \
    cls##Def::cls##Def()                                                                                               \
        : ParticleDef()                                                                                                \
    {                                                                                                                  \
        name     = #cls;                                                                                               \
        mass     = MASS;                                                                                               \
        lifetime = LIFETIME;                                                                                           \
        charge   = CHARGE;                                                                                             \
        if (lifetime < 0.0)                                                                                            \
        {                                                                                                              \
            decay_table.SetStable();                                                                                   \
        }                                                                                                              \
    }                                                                                                                  \
                                                                                                                       \
    cls##Def::~cls##Def() {}

using namespace PROPOSAL;

/******************************************************************************
*                                HardBBTables                                *
******************************************************************************/

const double HardBBTables::muon[8][7] = {
    { 7.174409e-4, 1.7132e-3, 4.082304e-3, 8.628455e-3, 0.01244159, 0.02204591, 0.03228755 },
    { -0.2436045, -0.5756682, -1.553973, -3.251305, -5.976818, -9.495636, -13.92918 },
    { -0.2942209, -0.68615, -2.004218, -3.999623, -6.855045, -10.05705, -14.37232 },
    { -0.1658391, -0.3825223, -1.207777, -2.33175, -3.88775, -5.636636, -8.418409 },
    { -0.05227727, -0.1196482, -0.4033373, -0.7614046, -1.270677, -1.883845, -2.948277 },
    { -9.328318e-3, -0.02124577, -0.07555636, -0.1402496, -0.2370768, -0.3614146, -0.5819409 },
    { -8.751909e-4, -1.987841e-3, -7.399682e-3, -0.01354059, -0.02325118, -0.03629659, -0.059275 },
    { -3.343145e-5, -7.584046e-5, -2.943396e-4, -5.3155e-4, -9.265136e-4, -1.473118e-3, -2.419946e-3 }
};

const double HardBBTables::tau[8][7] = {
    { -1.269205e-4, -2.843877e-4, -5.761546e-4, -1.195445e-3, -1.317386e-3, -9.689228e-15, -6.4595e-15 },
    { -0.01563032, -0.03589573, -0.07768545, -0.157375, -0.2720009, -0.4186136, -0.8045046 },
    { 0.04693954, 0.1162945, 0.3064255, 0.7041273, 1.440518, 2.533355, 3.217832 },
    { 0.05338546, 0.130975, 0.3410341, 0.7529364, 1.425927, 2.284968, 2.5487 },
    { 0.02240132, 0.05496, 0.144945, 0.3119032, 0.5576727, 0.8360727, 0.8085682 },
    { 4.658909e-3, 0.01146659, 0.03090286, 0.06514455, 0.1109868, 0.1589677, 0.1344223 },
    { 4.822364e-4, 1.193018e-3, 3.302773e-3, 6.843364e-3, 0.011191, 0.015614, 0.01173827 },
    { 1.9837e-5, 4.940182e-5, 1.409573e-4, 2.877909e-4, 4.544877e-4, 6.280818e-4, 4.281932e-4 }
};

HardBBTables::VecType HardBBTables::getBBVector(const double table[8][7])
{
    HardBBTables::VecType var;

    for (int i = 0; i < 8; ++i)
    {
        var.push_back(std::vector<double>(table[i], table[i] + sizeof(table[i]) / sizeof(double)));
    }

    return var;
}

const HardBBTables::VecType HardBBTables::MuonTable = HardBBTables::getBBVector(HardBBTables::muon);
const HardBBTables::VecType HardBBTables::TauTable = HardBBTables::getBBVector(HardBBTables::tau);


/******************************************************************************
*                                ParticleDef                                 *
******************************************************************************/

namespace PROPOSAL
{

std::ostream& operator<<(std::ostream& os, ParticleDef const& def)
{
    std::stringstream ss;
    ss << " ParticleDef (" << &def << ") ";

    os << Helper::Centered(60, ss.str()) << '\n';

    os << def.name << '\n';
    os << "Mass:"
       << "\t\t" << def.mass << '\n';
    os << "Low:"
       << "\t\t" << def.low << '\n';
    os << "Lifetime:"
       << "\t" << def.lifetime << '\n';
    os << "Charge:"
       << "\t\t" << def.charge << '\n';

    os << Helper::Centered(60, "");
    return os;
}

} /* PROPOSAL */


ParticleDef::ParticleDef()
    : name("")
    , mass(0.0)
    , low(0.0)
    , lifetime(0.0)
    , charge(0.0)
    , hardbb_table(NULL)
    , decay_table()
{
}

ParticleDef::ParticleDef(std::string name,
            double mass,
            double lifetime,
            double charge,
            const HardBBTables::VecType* table)
    : name(name)
    , mass(mass)
    , low(mass)
    , lifetime(lifetime)
    , charge(charge)
    , hardbb_table(table)
    , decay_table()
{
}

ParticleDef::~ParticleDef()
{
}

ParticleDef::ParticleDef(const ParticleDef& def)
{
    name = def.name;
    mass = def.mass;
    low = def.low;
    lifetime = def.lifetime;
    charge = def.charge;
    hardbb_table = def.hardbb_table;
    decay_table = def.decay_table;
}

ParticleDef& ParticleDef::operator=(const ParticleDef& def)
{
    if (this != &def)
    {
      ParticleDef tmp(def);
      swap(tmp);
    }
    return *this;
}

void ParticleDef::swap(ParticleDef& def)
{
    using std::swap;

    swap(name, def.name);
    swap(mass, def.mass);
    swap(low, def.low);
    swap(lifetime, def.lifetime);
    swap(charge, def.charge);
    swap(hardbb_table, def.hardbb_table);

    decay_table.swap(def.decay_table);
}

bool ParticleDef::operator==(const ParticleDef& def) const
{
    if (name != def.name)
    {
        return false;
    }
    else if (mass != def.mass)
    {
        return false;
    }
    else if (low != def.low)
    {
        return false;
    }
    else if (lifetime != def.lifetime)
    {
        return false;
    }
    else if (charge != def.charge)
    {
        return false;
    }
    else if (hardbb_table != def.hardbb_table)
    {
        return false;
    }
    else if (decay_table != def.decay_table)
    {
        return false;
    }
    else
    {
        return true;
    }
}

bool ParticleDef::operator!=(const ParticleDef& def) const
{
    return !(*this == def);
}

std::size_t PROPOSAL::hash_value(ParticleDef const& particle_def) {
    std::size_t seed = 0;
    boost::hash_combine(seed, particle_def.mass);
    boost::hash_combine(seed, particle_def.lifetime);
    boost::hash_combine(seed, particle_def.charge);
    return seed;
}

// ------------------------------------------------------------------------- //
// ParticleDef definitions
// ------------------------------------------------------------------------- //

MuMinusDef::MuMinusDef()
    : ParticleDef()
{
    name = "MuMinus";
    mass = MMU;
    low = MMU;
    lifetime = LMU;
    charge = -1.0;
    hardbb_table = &HardBBTables::MuonTable;

    // Decay modes
    LeptonicDecayChannel mode;
    decay_table.addChannel(1.0, mode);
}

MuMinusDef::~MuMinusDef()
{
}

MuPlusDef::MuPlusDef()
    : ParticleDef()
{
    name = "MuPlus";
    mass = MMU;
    low = MMU;
    lifetime = LMU;
    charge = -1.0;
    hardbb_table = &HardBBTables::MuonTable;

    // Decay modes
    LeptonicDecayChannel mode;
    decay_table.addChannel(1.0, mode);
}

MuPlusDef::~MuPlusDef()
{
}

TauMinusDef::TauMinusDef()
    : ParticleDef()
{
    name = "TauMinus";
    mass = MTAU;
    low = MTAU;
    lifetime = LTAU;
    charge = -1.0;
    hardbb_table = &HardBBTables::TauTable;

    // Decay modes
    DecayChannel* mode = new LeptonicDecayChannel();

    //TODO(mario): Different leptonic modes Mon 2017/08/21
    decay_table.addChannel(0.1737, *mode);
    decay_table.addChannel(0.1783, *mode);
    delete mode;

    mode = new TwoBodyPhaseSpace(MPI, 0.0);
    decay_table.addChannel(0.1109, *mode);
    delete mode;

    mode = new TwoBodyPhaseSpace(MRH, 0.0);
    decay_table.addChannel(0.2540, *mode);
    delete mode;

    mode = new TwoBodyPhaseSpace(MA1, 0.0);
    decay_table.addChannel(0.1826, *mode);
    delete mode;

    mode = new TwoBodyPhaseSpace(MRS, 0.0);
    decay_table.addChannel(1.0, *mode); // Else
    delete mode;
}

TauMinusDef::~TauMinusDef()
{
}

TauPlusDef::TauPlusDef()
    : ParticleDef()
{
    name = "TauPlus";
    mass = MMU;
    low = MMU;
    lifetime = LMU;
    charge = -1.0;
    hardbb_table = &HardBBTables::TauTable;

    // Decay modes
    DecayChannel* mode = new LeptonicDecayChannel();

    //TODO(mario): Different leptonic modes Mon 2017/08/21
    decay_table.addChannel(0.1737, *mode);
    decay_table.addChannel(0.1783, *mode);
    delete mode;

    mode = new TwoBodyPhaseSpace(MPI, 0.0);
    decay_table.addChannel(0.1109, *mode);
    delete mode;

    mode = new TwoBodyPhaseSpace(MRH, 0.0);
    decay_table.addChannel(0.2540, *mode);
    delete mode;

    mode = new TwoBodyPhaseSpace(MA1, 0.0);
    decay_table.addChannel(0.1826, *mode);
    delete mode;

    mode = new TwoBodyPhaseSpace(MRS, 0.0);
    decay_table.addChannel(1.0, *mode); // Else
    delete mode;
}

TauPlusDef::~TauPlusDef()
{
}

// ------------------------------------------------------------------------- //
// Signature for following macro definitions:
//
// class name, mass, lifetime, charge
//
// ------------------------------------------------------------------------- //

PARTICLE_IMP(EMinus, ME, STABLE_PARTICLE, -1.0)
PARTICLE_IMP(EPlus, ME, STABLE_PARTICLE, 1.0)

PARTICLE_IMP(StauMinus, MSTAU, STABLE_PARTICLE, -1.0)
PARTICLE_IMP(StauPlus, MSTAU, STABLE_PARTICLE, 1.0)

PARTICLE_IMP(P0, MPI0, LPI0, 0.0)
PARTICLE_IMP(PiMinus, MPI, LPI, -1.0)
PARTICLE_IMP(PiPlus, MPI, LPI, 1.0)

PARTICLE_IMP(KMinus, MKAON, LKAON, -1.0)
PARTICLE_IMP(KPlus, MKAON, LKAON, 1.0)

PARTICLE_IMP(PMinus, MP, STABLE_PARTICLE, -1.0)
PARTICLE_IMP(PPlus, MP, STABLE_PARTICLE, 1.0)

PARTICLE_IMP(NuE, 0.0, STABLE_PARTICLE, 0.0)
PARTICLE_IMP(NuEBar, 0.0, STABLE_PARTICLE, 0.0)

PARTICLE_IMP(NuMu, 0.0, STABLE_PARTICLE, 0.0)
PARTICLE_IMP(NuMuBar, 0.0, STABLE_PARTICLE, 0.0)

PARTICLE_IMP(NuTau, 0.0, STABLE_PARTICLE, 0.0)
PARTICLE_IMP(NuTauBar, 0.0, STABLE_PARTICLE, 0.0)

PARTICLE_IMP(DeltaE, 0.0, 0.0, 0.0)
PARTICLE_IMP(Brems, 0.0, 0.0, 0.0)
PARTICLE_IMP(NuclInt, 0.0, 0.0, 0.0)
PARTICLE_IMP(Hadrons, 0.0, 0.0, 0.0)
PARTICLE_IMP(Epair, 0.0, 0.0, 0.0)
PARTICLE_IMP(Mupair, 0.0, 0.0, 0.0)

PARTICLE_IMP(ContinuousEnergyLoss, 0.0, 0.0, 0.0)
PARTICLE_IMP(Monopole, MMON, STABLE_PARTICLE, CMON)
PARTICLE_IMP(Gamma, 0.0, STABLE_PARTICLE, 0.0)
PARTICLE_IMP(StableMassiveParticle, MSMP, STABLE_PARTICLE, -1.0)
