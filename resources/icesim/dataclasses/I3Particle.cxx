#include <boost/functional/hash/hash.hpp>
#include <dataclasses/I3Constants.h>
#include <dataclasses/physics/I3Particle.h>
#include <icetray/I3Units.h>
#include <icetray/serialization.h>

#include <boost/assign/list_of.hpp>
#include <limits>
#include <map>

#include <boost/lexical_cast.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/seq/for_each.hpp>

#if BOOST_VERSION >= 105300
#define USE_ATOMICS
#endif

#ifdef USE_ATOMICS
#include <boost/atomic.hpp>
#include <boost/thread/lock_guard.hpp>
#include <boost/thread/mutex.hpp>
#endif

#include <cstdlib>
#include <unistd.h>

#ifdef USE_ATOMICS
namespace {
boost::mutex global_id_lock;

boost::atomic<int32_t> global_last_pid_(0);
boost::atomic<int32_t> global_minor_id_(0);
boost::atomic<uint64_t> global_major_id_(0);
} // namespace
#else
static int32_t global_last_pid_  = 0;
static int32_t global_minor_id_  = 0;
static uint64_t global_major_id_ = 0;
#endif

void I3Particle::generateID()
{
    int this_pid = getpid();

#ifdef USE_ATOMICS // thread-safe version
    int32_t last_pid = global_last_pid_.load(boost::memory_order_relaxed);
    boost::atomic_thread_fence(boost::memory_order_acquire); // keep memory ops from wandering
    if (this_pid != last_pid)
    {
        boost::lock_guard<boost::mutex> lg(global_id_lock); // acquire the lock
        // check whether another thread already updated this
        last_pid = global_last_pid_.load(boost::memory_order_relaxed);
        if (this_pid != last_pid)
        {
            log_debug("PID has changed from %i to %i. regenerating I3Particle::majorID.", last_pid, this_pid);
            boost::atomic_thread_fence(boost::memory_order_release);
            global_last_pid_.store(this_pid, boost::memory_order_relaxed);
            global_major_id_.store(0, boost::memory_order_relaxed); // this will cause a new major ID to be generated
            global_minor_id_.store(0, boost::memory_order_relaxed); // reset the minor ID, too
        }
    }

    uint64_t old_major_id = global_major_id_.load(boost::memory_order_relaxed);
    boost::atomic_thread_fence(boost::memory_order_acquire); // keep memory ops from wandering
    if (old_major_id == 0)
    {
        boost::lock_guard<boost::mutex> lg(global_id_lock); // acquire the lock
        // check whether another thread already updated this
        old_major_id = global_major_id_.load(boost::memory_order_relaxed);
        if (old_major_id == 0)
        {
            boost::hash<std::string> string_hash;
            std::stringstream s;
            s << time(0) << this_pid << gethostid();
            boost::atomic_thread_fence(boost::memory_order_release);
            global_major_id_.store(string_hash(s.str()), boost::memory_order_relaxed);
        }
    }
    ID_.majorID = global_major_id_.load();
    ID_.minorID = global_minor_id_.fetch_add(1);
#else // unsafe version if atomics aren't available
    if (this_pid != global_last_pid_)
    {
        log_debug("PID has changed from %i to %i. regenerating I3Particle::majorID.", global_last_pid_, this_pid);
        global_last_pid_ = this_pid;
        global_major_id_ = 0; // this will cause a new major ID to be generated
        global_minor_id_ = 0; // reset the minor ID, too
    }
    if (global_major_id_ == 0)
    {
        boost::hash<std::string> string_hash;
        std::stringstream s;
        s << time(0) << this_pid << gethostid();
        global_major_id_ = string_hash(s.str());
    }
    ID_.majorID = global_major_id_;
    ID_.minorID = global_minor_id_++;
#endif
}

I3Particle::~I3Particle() {}
I3Particle::I3Particle(ParticleShape shape, ParticleType type)
    : pdgEncoding_(type)
    , shape_(shape)
    , status_(NotSet)
    , pos_()
    , dir_()
    , time_(NAN)
    , energy_(NAN)
    , length_(NAN)
    , speed_(I3Constants::c)
    , locationType_(Anywhere)
{
    generateID();

    log_trace_stream("Calling I3Particle::I3Particle(ParticleShape " << shape << ", ParticleType " << type << ").");
}

I3Particle::I3Particle(const I3Position pos,
                       const I3Direction dir,
                       const double vertextime,
                       ParticleShape shape,
                       double length)
    : pdgEncoding_(unknown)
    , shape_(shape)
    , status_(NotSet)
    , pos_(pos)
    , dir_(dir)
    , time_(vertextime)
    , energy_(NAN)
    , length_(length)
    , speed_(I3Constants::c)
    , locationType_(Anywhere)
{
    generateID();
}

I3Particle::I3Particle(const uint64_t major, const int32_t minor)
    : pdgEncoding_(unknown)
    , shape_(Null)
    , status_(NotSet)
    , pos_()
    , dir_()
    , time_(NAN)
    , energy_(NAN)
    , length_(NAN)
    , speed_(I3Constants::c)
    , locationType_(Anywhere)
{
    ID_.majorID = major;
    ID_.minorID = minor;
}

I3Particle I3Particle::Clone() const
{
    I3Particle p;
    p.SetPdgEncoding(pdgEncoding_);
    p.SetShape(shape_);
    p.SetFitStatus(status_);
    p.SetPos(pos_);
    p.SetDir(dir_);
    p.SetTime(time_);
    p.SetEnergy(energy_);
    p.SetLength(length_);
    p.SetSpeed(speed_);
    p.SetLocationType(locationType_);
    return p;
}

// the masses are taken from geant4.9.6.p03 except for WPlus, WMinus and Z0,
// which were taken from PDG Booklet July 2012
typedef std::map<I3Particle::ParticleType, double> particle_type_mass_conversion_t;
static const particle_type_mass_conversion_t fromParticleTypeMassTable =
    boost::assign::list_of<std::pair<I3Particle::ParticleType, double> >(I3Particle::Gamma, 0.0 * I3Units::GeV)(
        I3Particle::EPlus,
        0.00051099891 * I3Units::GeV)(I3Particle::EMinus, 0.00051099891 * I3Units::GeV)(
        I3Particle::MuPlus,
        0.105658367 * I3Units::GeV)(I3Particle::MuMinus, 0.105658367 * I3Units::GeV)(
        I3Particle::Pi0,
        0.1349766 * I3Units::GeV)(I3Particle::PiPlus, 0.1395701 * I3Units::GeV)(I3Particle::PiMinus,
                                                                                0.1395701 * I3Units::GeV)(
        I3Particle::K0_Long,
        0.497614 * I3Units::GeV)(I3Particle::KPlus, 0.493677 * I3Units::GeV)(I3Particle::KMinus,
                                                                             0.493677 * I3Units::GeV)(
        I3Particle::Neutron,
        0.93956536 * I3Units::GeV)(I3Particle::PPlus, 0.938272013 * I3Units::GeV)(I3Particle::PMinus,
                                                                                  0.938272013 * I3Units::GeV)(
        I3Particle::K0_Short,
        0.497614 * I3Units::GeV)(I3Particle::Eta, 0.547853 * I3Units::GeV)(I3Particle::Lambda, 1.115683 * I3Units::GeV)(
        I3Particle::SigmaPlus,
        1.18937 * I3Units::GeV)(I3Particle::Sigma0, 1.192642 * I3Units::GeV)(
        I3Particle::SigmaMinus,
        1.18937 * I3Units::GeV)(I3Particle::Xi0, 1.31486 * I3Units::GeV)(I3Particle::XiMinus, 1.32171 * I3Units::GeV)(
        I3Particle::OmegaMinus,
        1.67245 * I3Units::GeV)(I3Particle::NeutronBar, 0.93956536 * I3Units::GeV)(I3Particle::LambdaBar,
                                                                                   1.115683 * I3Units::GeV)(
        I3Particle::SigmaMinusBar,
        1.18937 * I3Units::GeV)(I3Particle::Sigma0Bar, 1.192642 * I3Units::GeV)(I3Particle::SigmaPlusBar,
                                                                                1.18937 * I3Units::GeV)(
        I3Particle::Xi0Bar,
        1.31486 * I3Units::GeV)(I3Particle::XiPlusBar, 1.32171 * I3Units::GeV)(I3Particle::OmegaPlusBar,
                                                                               1.67245 * I3Units::GeV)(
        I3Particle::DPlus,
        1.86957 * I3Units::GeV)(I3Particle::DMinus, 1.86957 * I3Units::GeV)(I3Particle::D0, 1.8648 * I3Units::GeV)(
        I3Particle::D0Bar,
        1.8648 * I3Units::GeV)(I3Particle::DsPlus, 1.96845 * I3Units::GeV)(I3Particle::DsMinusBar,
                                                                           1.96845 * I3Units::GeV)(
        I3Particle::LambdacPlus,
        2.28646 * I3Units::GeV)(I3Particle::WPlus, 80.385 * I3Units::GeV)(I3Particle::WMinus, 80.385 * I3Units::GeV)(
        I3Particle::Z0,
        91.1876 * I3Units::GeV)(I3Particle::NuE, 0.0 * I3Units::GeV)(I3Particle::NuEBar, 0.0 * I3Units::GeV)(
        I3Particle::NuMu,
        0.0 * I3Units::GeV)(I3Particle::NuMuBar, 0.0 * I3Units::GeV)(I3Particle::TauPlus, 1.77682 * I3Units::GeV)(
        I3Particle::TauMinus,
        1.77682 * I3Units::GeV)(I3Particle::NuTau, 0.0 * I3Units::GeV)(I3Particle::NuTauBar, 0.0 * I3Units::GeV)

    // Nuclei
    (I3Particle::He3Nucleus, 2.808391 * I3Units::GeV)(I3Particle::He4Nucleus, 3.727379 * I3Units::GeV)(
        I3Particle::Li6Nucleus,
        5.60151816372 * I3Units::GeV)(I3Particle::Li7Nucleus, 6.53383353972 * I3Units::GeV)(
        I3Particle::Be9Nucleus,
        8.39275030104 * I3Units::GeV)(I3Particle::B10Nucleus, 9.32443669262 * I3Units::GeV)(
        I3Particle::B11Nucleus,
        10.2525479206 * I3Units::GeV)(I3Particle::C12Nucleus, 11.174863388 * I3Units::GeV)(
        I3Particle::C13Nucleus,
        12.1094824273 * I3Units::GeV)(I3Particle::N14Nucleus, 13.0402043278 * I3Units::GeV)(
        I3Particle::N15Nucleus,
        13.9689363768 * I3Units::GeV)(I3Particle::O16Nucleus, 14.8950815346 * I3Units::GeV)(
        I3Particle::O17Nucleus,
        15.830503751 * I3Units::GeV)(I3Particle::O18Nucleus, 16.76202507 * I3Units::GeV)(
        I3Particle::F19Nucleus,
        17.6923029112 * I3Units::GeV)(I3Particle::Ne20Nucleus, 18.6177321841 * I3Units::GeV)(
        I3Particle::Ne21Nucleus,
        19.5505363674 * I3Units::GeV)(I3Particle::Ne22Nucleus, 20.4797374564 * I3Units::GeV)(
        I3Particle::Na23Nucleus,
        21.4092162538 * I3Units::GeV)(I3Particle::Mg24Nucleus, 22.3357965987 * I3Units::GeV)(
        I3Particle::Mg25Nucleus,
        23.2680313677 * I3Units::GeV)(I3Particle::Mg26Nucleus, 24.1965036397 * I3Units::GeV)(
        I3Particle::Al26Nucleus,
        24.1999980695 * I3Units::GeV)(I3Particle::Al27Nucleus, 25.1265057485 * I3Units::GeV)(
        I3Particle::Si28Nucleus,
        26.0531939251 * I3Units::GeV)(I3Particle::Si29Nucleus, 26.9842857039 * I3Units::GeV)(
        I3Particle::Si30Nucleus,
        27.9132418499 * I3Units::GeV)(I3Particle::Si31Nucleus, 28.8462197999 * I3Units::GeV)(
        I3Particle::Si32Nucleus,
        29.7765819269 * I3Units::GeV)(I3Particle::P31Nucleus, 28.8442183429 * I3Units::GeV)(
        I3Particle::P32Nucleus,
        29.7758480379 * I3Units::GeV)(I3Particle::P33Nucleus, 30.7053097979 * I3Units::GeV)(
        I3Particle::S32Nucleus,
        29.773628119 * I3Units::GeV)(I3Particle::S33Nucleus, 30.70455185 * I3Units::GeV)(I3Particle::S34Nucleus,
                                                                                         31.632700084 * I3Units::GeV)(
        I3Particle::S35Nucleus,
        32.565279544 * I3Units::GeV)(I3Particle::S36Nucleus, 33.494955853 * I3Units::GeV)(I3Particle::Cl35Nucleus,
                                                                                          32.5646030619 * I3Units::GeV)(
        I3Particle::Cl36Nucleus,
        33.4955887729 * I3Units::GeV)(I3Particle::Cl37Nucleus, 34.4248431259 * I3Units::GeV)(
        I3Particle::Ar36Nucleus,
        33.4943699372 * I3Units::GeV)(I3Particle::Ar37Nucleus, 34.4251478462 * I3Units::GeV)(
        I3Particle::Ar38Nucleus,
        35.3528749822 * I3Units::GeV)(I3Particle::Ar39Nucleus, 36.2858415502 * I3Units::GeV)(
        I3Particle::Ar40Nucleus,
        37.2155376931 * I3Units::GeV)(I3Particle::Ar41Nucleus, 38.1490041502 * I3Units::GeV)(
        I3Particle::Ar42Nucleus,
        39.0791429702 * I3Units::GeV)(I3Particle::K39Nucleus, 36.284767546 * I3Units::GeV)(
        I3Particle::K40Nucleus,
        37.21653338 * I3Units::GeV)(I3Particle::K41Nucleus, 38.146003539 * I3Units::GeV)(I3Particle::Ca40Nucleus,
                                                                                         37.2147134578 * I3Units::GeV)(
        I3Particle::Ca41Nucleus,
        38.1459160018 * I3Units::GeV)(I3Particle::Ca42Nucleus, 39.0740007168 * I3Units::GeV)(
        I3Particle::Ca43Nucleus,
        40.0056331778 * I3Units::GeV)(I3Particle::Ca44Nucleus, 40.9340673658 * I3Units::GeV)(
        I3Particle::Ca45Nucleus,
        41.8662179228 * I3Units::GeV)(I3Particle::Ca46Nucleus, 42.7953888238 * I3Units::GeV)(
        I3Particle::Ca47Nucleus,
        43.7276778058 * I3Units::GeV)(I3Particle::Ca48Nucleus, 44.6572978278 * I3Units::GeV)(
        I3Particle::Sc44Nucleus,
        40.9372110547 * I3Units::GeV)(I3Particle::Sc45Nucleus, 41.8654533837 * I3Units::GeV)(
        I3Particle::Sc46Nucleus,
        42.7962580887 * I3Units::GeV)(I3Particle::Sc47Nucleus, 43.7251771107 * I3Units::GeV)(
        I3Particle::Sc48Nucleus,
        44.6565071587 * I3Units::GeV)(I3Particle::Ti44Nucleus, 40.9369701498 * I3Units::GeV)(
        I3Particle::Ti45Nucleus,
        41.8670068998 * I3Units::GeV)(I3Particle::Ti46Nucleus, 42.7933832428 * I3Units::GeV)(
        I3Particle::Ti47Nucleus,
        43.7240682988 * I3Units::GeV)(I3Particle::Ti48Nucleus, 44.6520069938 * I3Units::GeV)(
        I3Particle::Ti49Nucleus,
        45.5834299498 * I3Units::GeV)(I3Particle::Ti50Nucleus, 46.5120561048 * I3Units::GeV)(
        I3Particle::V48Nucleus,
        44.6555109582 * I3Units::GeV)(I3Particle::V49Nucleus, 45.5835234282 * I3Units::GeV)(
        I3Particle::V50Nucleus,
        46.5137528452 * I3Units::GeV)(I3Particle::V51Nucleus, 47.4422670442 * I3Units::GeV)(
        I3Particle::Cr50Nucleus,
        46.5122066868 * I3Units::GeV)(I3Particle::Cr51Nucleus, 47.4425114068 * I3Units::GeV)(
        I3Particle::Cr52Nucleus,
        48.3700373088 * I3Units::GeV)(I3Particle::Cr53Nucleus, 49.3016635288 * I3Units::GeV)(
        I3Particle::Cr54Nucleus,
        50.2315097528 * I3Units::GeV)(I3Particle::Mn52Nucleus, 48.3742407516 * I3Units::GeV)(
        I3Particle::Mn53Nucleus,
        49.3017523196 * I3Units::GeV)(I3Particle::Mn54Nucleus, 50.2323788816 * I3Units::GeV)(
        I3Particle::Mn55Nucleus,
        51.1617176996 * I3Units::GeV)(I3Particle::Fe54Nucleus, 50.2311739195 * I3Units::GeV)(
        I3Particle::Fe55Nucleus,
        51.1614410355 * I3Units::GeV)(I3Particle::Fe56Nucleus, 52.0898090795 * I3Units::GeV)(
        I3Particle::Fe57Nucleus,
        53.0217283295 * I3Units::GeV)(I3Particle::Fe58Nucleus, 53.9512490695 * I3Units::GeV);

bool I3Particle::HasMass() const
{
    ParticleType type                                  = ParticleType(pdgEncoding_);
    particle_type_mass_conversion_t::const_iterator it = fromParticleTypeMassTable.find(type);
    if (it == fromParticleTypeMassTable.end())
    {
        return false;
    } else
    {
        return true;
    }
}

double I3Particle::GetMass() const
{
    ParticleType type = ParticleType(pdgEncoding_);
    if (this->HasMass())
    {
        return fromParticleTypeMassTable.at(type);
    } else
    {
        log_fatal("\"%s\" has no mass implemented.", (this->GetTypeString()).c_str());
    }
}

double I3Particle::GetMassForType(ParticleType type)
{
    particle_type_mass_conversion_t::const_iterator it = fromParticleTypeMassTable.find(type);
    if (it == fromParticleTypeMassTable.end())
    {
        log_fatal("\"%d\" has no mass implemented.", type);
    } else
    {
        return it->second;
    }
}

double I3Particle::GetTotalEnergy() const
{
    ParticleType type = ParticleType(pdgEncoding_);
    if (this->HasMass())
    {
        double mass = fromParticleTypeMassTable.at(type);
        return energy_ + mass;
    } else
    {
        log_fatal("\"%s\" has no mass implemented. Cannot get total energy.", (this->GetTypeString()).c_str());
    }
}

void I3Particle::SetTotalEnergy(double total_energy)
{
    ParticleType type = ParticleType(pdgEncoding_);
    double energy;
    if (this->HasMass())
    {
        energy = total_energy - fromParticleTypeMassTable.at(type);
    } else
    {
        log_fatal("\"%s\" has no mass implemented. Cannot set total energy.", (this->GetTypeString()).c_str());
    }
    if (energy < 0.)
    {
        log_fatal("Particle must not have negative energy.");
    }
    this->SetEnergy(energy);
}

// using the magic of the preprocessor, expand
// the existing list of enum entries into a case
// line converting from enum to string
#define MAKE_ENUM_TO_STRING_CASE_LINE(r, data, t)                                                                      \
    case t:                                                                                                            \
        return BOOST_PP_STRINGIZE(t);

std::string I3Particle::GetTypeString() const
{
    switch (pdgEncoding_)
    {
        BOOST_PP_SEQ_FOR_EACH(MAKE_ENUM_TO_STRING_CASE_LINE, ~, I3PARTICLE_H_I3Particle_ParticleType)
    }
    return (boost::lexical_cast<std::string>(pdgEncoding_));
}

std::string I3Particle::GetShapeString() const
{
    switch (shape_)
    {
        BOOST_PP_SEQ_FOR_EACH(MAKE_ENUM_TO_STRING_CASE_LINE, ~, I3PARTICLE_H_I3Particle_ParticleShape)
    }
    return (boost::lexical_cast<std::string>(shape_));
}

std::string I3Particle::GetFitStatusString() const
{
    switch (status_)
    {
        BOOST_PP_SEQ_FOR_EACH(MAKE_ENUM_TO_STRING_CASE_LINE, ~, I3PARTICLE_H_I3Particle_FitStatus)
    }
    return (boost::lexical_cast<std::string>(status_));
}

std::string I3Particle::GetLocationTypeString() const
{
    switch (locationType_)
    {
        BOOST_PP_SEQ_FOR_EACH(MAKE_ENUM_TO_STRING_CASE_LINE, ~, I3PARTICLE_H_I3Particle_LocationType)
    }
    return (boost::lexical_cast<std::string>(locationType_));
}

#define MAKE_STRING_TO_ENUM_IF_LINE(r, data, t)                                                                        \
    else if (str == BOOST_PP_STRINGIZE(t)) { data = t; }

void I3Particle::SetTypeString(const std::string& str)
{
    ParticleType type;

    if (false)
    {
    }
    BOOST_PP_SEQ_FOR_EACH(MAKE_STRING_TO_ENUM_IF_LINE, type, I3PARTICLE_H_I3Particle_ParticleType)
    else
    {
        try
        {
            type = static_cast<ParticleType>(boost::lexical_cast<int>(str));
        } catch (boost::bad_lexical_cast& bad)
        {
            log_fatal("\"%s\" is not a valid ParticleType.", str.c_str());
        }
    }

    SetType(type);
}

void I3Particle::SetShapeString(const std::string& str)
{
    if (false)
    {
    }
    BOOST_PP_SEQ_FOR_EACH(MAKE_STRING_TO_ENUM_IF_LINE, shape_, I3PARTICLE_H_I3Particle_ParticleShape)
    else
    {
        try
        {
            shape_ = static_cast<ParticleShape>(boost::lexical_cast<int>(str));
        } catch (boost::bad_lexical_cast& bad)
        {
            log_fatal("\"%s\" is not a valid ParticleShape.", str.c_str());
        }
    }
}

void I3Particle::SetFitStatusString(const std::string& str)
{
    if (false)
    {
    }
    BOOST_PP_SEQ_FOR_EACH(MAKE_STRING_TO_ENUM_IF_LINE, status_, I3PARTICLE_H_I3Particle_FitStatus)
    else
    {
        try
        {
            status_ = static_cast<FitStatus>(boost::lexical_cast<int>(str));
        } catch (boost::bad_lexical_cast& bad)
        {
            log_fatal("\"%s\"is not a valid FitStatus.", str.c_str());
        }
    }
}

void I3Particle::SetLocationTypeString(const std::string& str)
{
    if (false)
    {
    }
    BOOST_PP_SEQ_FOR_EACH(MAKE_STRING_TO_ENUM_IF_LINE, locationType_, I3PARTICLE_H_I3Particle_LocationType)
    else
    {
        try
        {
            locationType_ = static_cast<LocationType>(boost::lexical_cast<int>(str));
        } catch (boost::bad_lexical_cast& bad)
        {
            log_fatal("\"%s\" is not a valid LocationType.", str.c_str());
        }
    }
}

bool I3Particle::IsNucleus() const
{
    return (abs(pdgEncoding_) >= 1000000000 && abs(pdgEncoding_) <= 1099999999);
}

bool I3Particle::IsTrack() const
{
    if (shape_ == InfiniteTrack || shape_ == StartingTrack || shape_ == StoppingTrack || shape_ == ContainedTrack ||
        pdgEncoding_ == MuPlus || pdgEncoding_ == MuMinus || pdgEncoding_ == TauPlus || pdgEncoding_ == TauMinus ||
        pdgEncoding_ == STauPlus || pdgEncoding_ == STauMinus || pdgEncoding_ == SMP || pdgEncoding_ == Monopole ||
        (shape_ == Primary &&
         (pdgEncoding_ == PPlus || pdgEncoding_ == PMinus || IsNucleus() || pdgEncoding_ == Gamma)))
        return true;
    else
        return false;
}

bool I3Particle::IsCascade() const
{
    const int32_t type = pdgEncoding_;

    if (shape_ == Cascade || shape_ == CascadeSegment || type == EPlus || type == EMinus || type == Brems ||
        type == DeltaE || type == PairProd || type == NuclInt || type == Hadrons || type == Pi0 || type == PiPlus ||
        type == PiMinus || (shape_ != Primary && (type == PPlus || type == PMinus || IsNucleus() || type == Gamma)))
        return true;
    else
        return false;
}

bool I3Particle::IsNeutrino() const
{
    const ParticleType type = ParticleType(pdgEncoding_);

    if (type == NuE || type == NuEBar || type == NuMu || type == NuMuBar || type == NuTau || type == NuTauBar ||
        type == Nu)
        return true;
    else
        return false;
}

bool I3Particle::HasPosition() const
{
    if (std::isnan(GetX()) || std::isnan(GetY()) || std::isnan(GetZ()))
        return false;
    else
        return true;
}

bool I3Particle::HasDirection() const
{
    if (std::isnan(GetZenith()) || std::isnan(GetAzimuth()))
        return false;
    else
        return true;
}

bool I3Particle::HasEnergy() const
{
    if (std::isnan(energy_))
        return false;
    else
        return true;
}

bool I3Particle::IsTopShower() const
{
    if (shape_ == TopShower)
        return true;
    else
        return false;
}

I3Position I3Particle::ShiftTimeTrack(const double time) const
{
    if (std::isfinite(speed_) && speed_ > 0.)
    {
        return pos_ + (speed_ * time) * dir_;
    } else
    {
        log_error("ShiftTimeTrack needs the particle to have a valid speed");
        return I3Position();
    }
}

I3Position I3Particle::GetStartPos() const
{
    if (shape_ == StartingTrack || shape_ == ContainedTrack || shape_ == MCTrack)
        return pos_;
    else
    {
        log_warn("GetStartPos undefined for a particle that is neither starting "
                 "nor contained.");
        return I3Position();
    }
}

double I3Particle::GetStartTime() const
{
    if (shape_ == StartingTrack || shape_ == ContainedTrack || shape_ == MCTrack)
        return time_;
    else
    {
        log_warn("GetStartTime undefined for a particle that is neither starting "
                 "nor contained.");
        return NAN;
    }
}

I3Position I3Particle::GetStopPos() const
{
    if (shape_ == StoppingTrack)
        return pos_;
    else if (shape_ == ContainedTrack || shape_ == MCTrack)
    {
        return pos_ + length_ * dir_;
    } else
    {
        log_warn("GetStopPos undefined for a particle that is neither stopping "
                 "nor contained.");
        I3Position nullpos;
        return nullpos;
    }
}

double I3Particle::GetStopTime() const
{
    if (shape_ == StoppingTrack)
        return time_;
    else if (shape_ == ContainedTrack || shape_ == MCTrack)
    {
        return time_ + length_ / speed_;
    } else
    {
        log_warn("GetStopTime undefined for a particle that is neither stopping "
                 "nor contained.");
        return NAN;
    }
}

// XXX: These don't need to be bimaps. The look is only one way and they could
// be regular maps easily enough.

// Old IceCube numbering conventions
typedef std::map<int, I3Particle::ParticleType> particle_type_conversion_t;
static const particle_type_conversion_t fromOldI3ParticleTable =
    boost::assign::list_of<std::pair<int, I3Particle::ParticleType> >(0, I3Particle::unknown)(1, I3Particle::Gamma)(
        2,
        I3Particle::EPlus)(3, I3Particle::EMinus)(5, I3Particle::MuPlus)(6, I3Particle::MuMinus)(7, I3Particle::Pi0)(
        8,
        I3Particle::PiPlus)(9, I3Particle::PiMinus)(10, I3Particle::K0_Long)(11, I3Particle::KPlus)(
        12,
        I3Particle::KMinus)(13, I3Particle::Neutron)(14, I3Particle::PPlus)(15, I3Particle::PMinus)(
        16,
        I3Particle::K0_Short)(17, I3Particle::Eta)(18, I3Particle::Lambda)(19, I3Particle::SigmaPlus)(
        20,
        I3Particle::Sigma0)(21, I3Particle::SigmaMinus)(22, I3Particle::Xi0)(23, I3Particle::XiMinus)(
        24,
        I3Particle::OmegaMinus)(25, I3Particle::NeutronBar)(26, I3Particle::LambdaBar)(27, I3Particle::SigmaMinusBar)(
        28,
        I3Particle::Sigma0Bar)(29, I3Particle::SigmaPlusBar)(30, I3Particle::Xi0Bar)(31, I3Particle::XiPlusBar)(
        32,
        I3Particle::OmegaPlusBar)(35, I3Particle::DPlus)(36, I3Particle::DMinus)(37, I3Particle::D0)(
        38,
        I3Particle::D0Bar)(39, I3Particle::DsPlus)(40, I3Particle::DsMinusBar)(41, I3Particle::LambdacPlus)(
        42,
        I3Particle::WPlus)(43, I3Particle::WMinus)(44, I3Particle::Z0)(66, I3Particle::NuE)(67, I3Particle::NuEBar)(
        68,
        I3Particle::NuMu)(69, I3Particle::NuMuBar)(131, I3Particle::TauPlus)(132, I3Particle::TauMinus)(
        133,
        I3Particle::NuTau)(134, I3Particle::NuTauBar)(302, I3Particle::He3Nucleus)(402, I3Particle::He4Nucleus)(
        603,
        I3Particle::Li6Nucleus)(703, I3Particle::Li7Nucleus)(904, I3Particle::Be9Nucleus)(1005, I3Particle::B10Nucleus)(
        1105,
        I3Particle::B11Nucleus)(1206, I3Particle::C12Nucleus)(1306, I3Particle::C13Nucleus)(
        1407,
        I3Particle::N14Nucleus)(1507, I3Particle::N15Nucleus)(1608, I3Particle::O16Nucleus)(
        1708,
        I3Particle::O17Nucleus)(1808, I3Particle::O18Nucleus)(1909, I3Particle::F19Nucleus)(
        2010,
        I3Particle::Ne20Nucleus)(2110, I3Particle::Ne21Nucleus)(2210, I3Particle::Ne22Nucleus)(
        2311,
        I3Particle::Na23Nucleus)(2412, I3Particle::Mg24Nucleus)(2512, I3Particle::Mg25Nucleus)(
        2612,
        I3Particle::Mg26Nucleus)(2613, I3Particle::Al26Nucleus)(2713, I3Particle::Al27Nucleus)(
        2814,
        I3Particle::Si28Nucleus)(2914, I3Particle::Si29Nucleus)(3014, I3Particle::Si30Nucleus)(
        3114,
        I3Particle::Si31Nucleus)(3214, I3Particle::Si32Nucleus)(3115, I3Particle::P31Nucleus)(
        3215,
        I3Particle::P32Nucleus)(3315, I3Particle::P33Nucleus)(3216, I3Particle::S32Nucleus)(
        3316,
        I3Particle::S33Nucleus)(3416, I3Particle::S34Nucleus)(3516, I3Particle::S35Nucleus)(
        3616,
        I3Particle::S36Nucleus)(3517, I3Particle::Cl35Nucleus)(3617, I3Particle::Cl36Nucleus)(
        3717,
        I3Particle::Cl37Nucleus)(3618, I3Particle::Ar36Nucleus)(3718, I3Particle::Ar37Nucleus)(
        3818,
        I3Particle::Ar38Nucleus)(3918, I3Particle::Ar39Nucleus)(4018, I3Particle::Ar40Nucleus)(
        4118,
        I3Particle::Ar41Nucleus)(4218, I3Particle::Ar42Nucleus)(3919, I3Particle::K39Nucleus)(
        4019,
        I3Particle::K40Nucleus)(4119, I3Particle::K41Nucleus)(4020, I3Particle::Ca40Nucleus)(
        4120,
        I3Particle::Ca41Nucleus)(4220, I3Particle::Ca42Nucleus)(4320, I3Particle::Ca43Nucleus)(
        4420,
        I3Particle::Ca44Nucleus)(4520, I3Particle::Ca45Nucleus)(4620, I3Particle::Ca46Nucleus)(
        4720,
        I3Particle::Ca47Nucleus)(4820, I3Particle::Ca48Nucleus)(4421, I3Particle::Sc44Nucleus)(
        4521,
        I3Particle::Sc45Nucleus)(4621, I3Particle::Sc46Nucleus)(4721, I3Particle::Sc47Nucleus)(
        4821,
        I3Particle::Sc48Nucleus)(4422, I3Particle::Ti44Nucleus)(4522, I3Particle::Ti45Nucleus)(
        4622,
        I3Particle::Ti46Nucleus)(4722, I3Particle::Ti47Nucleus)(4822, I3Particle::Ti48Nucleus)(
        4922,
        I3Particle::Ti49Nucleus)(5022, I3Particle::Ti50Nucleus)(4823, I3Particle::V48Nucleus)(
        4923,
        I3Particle::V49Nucleus)(5023, I3Particle::V50Nucleus)(5123, I3Particle::V51Nucleus)(
        5024,
        I3Particle::Cr50Nucleus)(5124, I3Particle::Cr51Nucleus)(5224, I3Particle::Cr52Nucleus)(
        5324,
        I3Particle::Cr53Nucleus)(5424, I3Particle::Cr54Nucleus)(5225, I3Particle::Mn52Nucleus)(
        5325,
        I3Particle::Mn53Nucleus)(5425, I3Particle::Mn54Nucleus)(5525, I3Particle::Mn55Nucleus)(
        5426,
        I3Particle::Fe54Nucleus)(5526, I3Particle::Fe55Nucleus)(5626, I3Particle::Fe56Nucleus)(
        5726,
        I3Particle::Fe57Nucleus)(5826, I3Particle::Fe58Nucleus)(9900, I3Particle::CherenkovPhoton)(-4, I3Particle::Nu)(
        -41,
        I3Particle::Monopole)(-1001, I3Particle::Brems)(-1002, I3Particle::DeltaE)(-1003, I3Particle::PairProd)(
        -1004,
        I3Particle::NuclInt)(-1005, I3Particle::MuPair)(-1006, I3Particle::Hadrons)(
        -1111,
        I3Particle::ContinuousEnergyLoss)(-2100, I3Particle::FiberLaser)(-2101, I3Particle::N2Laser)(
        -2201,
        I3Particle::YAGLaser)(-9131, I3Particle::STauPlus)(-9132, I3Particle::STauMinus)(-9500, I3Particle::SMP);

static const particle_type_conversion_t fromRDMCTable =
    boost::assign::list_of<std::pair<int, I3Particle::ParticleType> >(-100, I3Particle::unknown)(1, I3Particle::Gamma)(
        2,
        I3Particle::EPlus)(3, I3Particle::EMinus)(4, I3Particle::Nu)(5, I3Particle::MuPlus)(6, I3Particle::MuMinus)(
        7,
        I3Particle::Pi0)(8, I3Particle::PiPlus)(9, I3Particle::PiMinus)(11, I3Particle::KPlus)(12, I3Particle::KMinus)(
        14,
        I3Particle::PPlus)(15, I3Particle::PMinus)(33, I3Particle::TauPlus)(34, I3Particle::TauMinus)(
        41,
        I3Particle::Monopole)(201, I3Particle::NuE)(202, I3Particle::NuMu)(203, I3Particle::NuTau)(
        204,
        I3Particle::NuEBar)(205, I3Particle::NuMuBar)(206, I3Particle::NuTauBar)(1001, I3Particle::Brems)(
        1002,
        I3Particle::DeltaE)(1003, I3Particle::PairProd)(1004, I3Particle::NuclInt)(1005, I3Particle::MuPair)(
        1006,
        I3Particle::Hadrons)(2100, I3Particle::FiberLaser)(2101, I3Particle::N2Laser)(2201, I3Particle::YAGLaser)(
        3002,
        I3Particle::He4Nucleus)(3003, I3Particle::Li7Nucleus)(3004, I3Particle::Be9Nucleus)(
        3005,
        I3Particle::B11Nucleus)(3006, I3Particle::C12Nucleus)(3007, I3Particle::N14Nucleus)(
        3008,
        I3Particle::O16Nucleus)(3009, I3Particle::F19Nucleus)(3010, I3Particle::Ne20Nucleus)(
        3011,
        I3Particle::Na23Nucleus)(3012, I3Particle::Mg24Nucleus)(3013, I3Particle::Al27Nucleus)(
        3014,
        I3Particle::Si28Nucleus)(3015, I3Particle::P31Nucleus)(3016, I3Particle::S32Nucleus)(
        3017,
        I3Particle::Cl35Nucleus)(3018, I3Particle::Ar40Nucleus)(3019, I3Particle::K39Nucleus)(
        3020,
        I3Particle::Ca40Nucleus)(3021, I3Particle::Sc45Nucleus)(3022, I3Particle::Ti48Nucleus)(
        3023,
        I3Particle::V51Nucleus)(3024, I3Particle::Cr52Nucleus)(3025, I3Particle::Mn55Nucleus)(3026,
                                                                                              I3Particle::Fe56Nucleus);

#ifndef __CINT__
I3Particle::I3Particle(const boost::optional<I3Particle>& p)
    : ID_(p->ID_)
    , pdgEncoding_(p->pdgEncoding_)
    , shape_(p->shape_)
    , status_(p->status_)
    , pos_(p->pos_)
    , dir_(p->dir_)
    , time_(p->time_)
    , energy_(p->energy_)
    , length_(p->length_)
    , speed_(p->speed_)
    , locationType_(p->locationType_)
{
}
#endif

template<class Archive>
void I3Particle::save(Archive& ar, unsigned version) const
{
    ar& make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar& make_nvp("ID", ID_.minorID);
    ar& make_nvp("MajorID", ID_.majorID);
    ar& make_nvp("pdgEncoding", pdgEncoding_);
    ar& make_nvp("shape", shape_);
    ar& make_nvp("fitStatus", status_);

    ar& make_nvp("pos", pos_);
    ar& make_nvp("dir", dir_);
    ar& make_nvp("time", time_);

    ar& make_nvp("energy", energy_);
    ar& make_nvp("length", length_);
    ar& make_nvp("speed", speed_);
    ar& make_nvp("LocationType", locationType_);
}

template<class Archive>
void I3Particle::load(Archive& ar, unsigned version)
{
    if (version > i3particle_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3Particle class.",
                  version,
                  i3particle_version_);

    ar& make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar& make_nvp("ID", ID_.minorID);
    if (version > 1)
    {
        ar& make_nvp("MajorID", ID_.majorID);
    }
    if (version == 0)
    {
        int32_t junk;
        ar& make_nvp("parentID", junk);
        ar& make_nvp("primaryID", junk);
    }
    if (version <= 2)
    {
        int t;
        ar& make_nvp("type", t);
        particle_type_conversion_t::const_iterator it = fromRDMCTable.find(t);

        if (it == fromRDMCTable.end())
        {
            log_warn("unknown RDMC code \"%i\" cannot be converted to a I3Particle::ParticleType. It will appear as "
                     "\"unknown\".",
                     t);
            pdgEncoding_ = unknown;
        } else
        {
            pdgEncoding_ = it->second;
        }
    } else if (version <= 4)
    {
        int t;
        ar& make_nvp("type", t);
        particle_type_conversion_t::const_iterator it = fromOldI3ParticleTable.find(t);

        if (it == fromOldI3ParticleTable.end())
        {
            log_warn(
                "unknown code \"%i\" cannot be converted to a I3Particle::ParticleType. It will appear as \"unknown\".",
                t);
            pdgEncoding_ = unknown;
        } else
        {
            pdgEncoding_ = it->second;
        }
    } else
    { // version >= 5
        ar& make_nvp("pdgEncoding", pdgEncoding_);
    }
    ar& make_nvp("shape", shape_);
    ar& make_nvp("fitStatus", status_);

    ar& make_nvp("pos", pos_);
    ar& make_nvp("dir", dir_);

    ar& make_nvp("time", time_);
    ar& make_nvp("energy", energy_);
    ar& make_nvp("length", length_);
    ar& make_nvp("speed", speed_);

    if (version == 0)
    {
        std::vector<I3Particle> junk;
        ar& make_nvp("composite", junk);
    }
    if (version > 0)
        ar& make_nvp("LocationType", locationType_);
    if (version == 4)
    {
        // obscure version in use by Antares, contains bjorken x and y.
        // Those should never have been in I3Particle. Load them and
        // forget them.
        double bjorkenx, bjorkeny;
        ar& make_nvp("bjorkenx", bjorkenx);
        ar& make_nvp("bjorkeny", bjorkeny);
    }
}

// specialize save and load for XML archives in order to display
// enums as strings instead of their numerical value
// (except for the particleType which is ignored when reading,
// as it is derived from the pdgEncoding)
template<>
void I3Particle::save(icecube::archive::xml_oarchive& ar, unsigned version) const
{
    std::string tempString;

    ar& make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar& make_nvp("minorID", ID_.minorID);
    ar& make_nvp("majorID", ID_.majorID);

    ar& make_nvp("pdgEncoding", pdgEncoding_);
    tempString = GetTypeString();
    ar& make_nvp("particleType", tempString);

    tempString = GetShapeString();
    ar& make_nvp("shape", tempString);

    tempString = GetFitStatusString();
    ar& make_nvp("fitStatus", tempString);

    ar& make_nvp("pos", pos_);
    ar& make_nvp("dir", dir_);
    ar& make_nvp("time", time_);
    ar& make_nvp("energy", energy_);
    ar& make_nvp("length", length_);
    ar& make_nvp("speed", speed_);

    tempString = GetLocationTypeString();
    ar& make_nvp("locationType", tempString);
}

template<>
void I3Particle::load(icecube::archive::xml_iarchive& ar, unsigned version)
{
    if (version != i3particle_version_)
        log_fatal("Cannot load XML data for I3Particle from an archive with version %u. Only the current version (%u) "
                  "is supported.",
                  version,
                  i3particle_version_);

    std::string tempString;

    ar& make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar& make_nvp("minorID", ID_.minorID);
    ar& make_nvp("majorID", ID_.majorID);

    ar& make_nvp("pdgEncoding", pdgEncoding_);
    ar& make_nvp("particleType", tempString); // ignored (I3Particle pdgEncoding)

    ar& make_nvp("shape", tempString);
    SetShapeString(tempString);

    ar& make_nvp("fitStatus", tempString);
    SetFitStatusString(tempString);

    ar& make_nvp("pos", pos_);
    ar& make_nvp("dir", dir_);
    ar& make_nvp("time", time_);
    ar& make_nvp("energy", energy_);
    ar& make_nvp("length", length_);
    ar& make_nvp("speed", speed_);

    ar& make_nvp("locationType", tempString);
    SetLocationTypeString(tempString);
}

std::ostream& operator<<(std::ostream& oss, const I3Particle& p)
{
    oss << "[ I3Particle MajorID : " << p.GetMajorID() << std::endl
        << "             MinorID : " << p.GetMinorID() << std::endl
        << "              Zenith : " << p.GetZenith() << std::endl
        << "             Azimuth : " << p.GetAzimuth() << std::endl
        << "                   X : " << p.GetX() << std::endl
        << "                   Y : " << p.GetY() << std::endl
        << "                   Z : " << p.GetZ() << std::endl
        << "                Time : " << p.GetTime() << std::endl
        << "              Energy : " << p.GetEnergy() << std::endl
        << "               Speed : " << p.GetSpeed() << std::endl
        << "              Length : " << p.GetLength() << std::endl
        << "                Type : " << p.GetTypeString() << std::endl
        << "        PDG encoding : " << p.GetPdgEncoding() << std::endl
        << "               Shape : " << p.GetShapeString() << std::endl
        << "              Status : " << p.GetFitStatusString() << std::endl
        << "            Location : " << p.GetLocationTypeString() << std::endl
        << "]";
    return oss;
}

std::ostream& operator<<(std::ostream& oss, const I3ParticleID& pid)
{
    oss << "I3ParticleID(" << pid.majorID << ", " << pid.minorID << ")";
    return oss;
}

I3_SPLIT_SERIALIZABLE(I3Particle);

// The name passed to this macro is set for
// all eternity.  Since all past I3Vector<I3Particle>
// containers were serialized as I3ParticleVect,
// this is the way it has to always be.
typedef I3Vector<I3Particle> I3ParticleVect;
I3_SERIALIZABLE(I3ParticleVect);
I3_SERIALIZABLE(I3ParticleID);
