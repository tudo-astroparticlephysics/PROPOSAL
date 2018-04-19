/**
    copyright  (C) 2004
    the icecube collaboration
    @version $Id: I3Particle.h 151491 2016-11-14 17:44:55Z kjmeagher $
    @file I3Particle.h
    @date $Date: 2016-11-14 18:44:55 +0100 (Mo, 14. Nov 2016) $
*/

#ifndef I3PARTICLE_H_INCLUDED
#define I3PARTICLE_H_INCLUDED

#include "dataclasses/I3Direction.h"
#include "dataclasses/I3Position.h"
#include "dataclasses/I3Vector.h"
#include "dataclasses/external/CompareFloatingPoint.h"
#include "dataclasses/physics/I3ParticleID.h"
#include "icetray/I3Units.h"
#include <map>
#include <string>

#ifndef __CINT__
#include <archive/xml_iarchive.hpp>
#include <archive/xml_oarchive.hpp>
#include <boost/optional.hpp>
#endif

/**
 * @brief
 */
class I3Particle : public I3FrameObject
{
public:
    enum ParticleType
    { // NB: These match the PDG codes. Keep it that way!
        unknown       = 0,
        Gamma         = 22,
        EPlus         = -11,
        EMinus        = 11,
        MuPlus        = -13,
        MuMinus       = 13,
        Pi0           = 111,
        PiPlus        = 211,
        PiMinus       = -211,
        K0_Long       = 130,
        KPlus         = 321,
        KMinus        = -321,
        Neutron       = 2112,
        PPlus         = 2212,
        PMinus        = -2212,
        K0_Short      = 310,
        Eta           = 221,
        Lambda        = 3122,
        SigmaPlus     = 3222,
        Sigma0        = 3212,
        SigmaMinus    = 3112,
        Xi0           = 3322,
        XiMinus       = 3312,
        OmegaMinus    = 3334,
        NeutronBar    = -2112,
        LambdaBar     = -3122,
        SigmaMinusBar = -3222,
        Sigma0Bar     = -3212,
        SigmaPlusBar  = -3112,
        Xi0Bar        = -3322,
        XiPlusBar     = -3312,
        OmegaPlusBar  = -3334,
        DPlus         = 411,
        DMinus        = -411,
        D0            = 421,
        D0Bar         = -421,
        DsPlus        = 431,
        DsMinusBar    = -431,
        LambdacPlus   = 4122,
        WPlus         = 24,
        WMinus        = -24,
        Z0            = 23,
        NuE           = 12,
        NuEBar        = -12,
        NuMu          = 14,
        NuMuBar       = -14,
        TauPlus       = -15,
        TauMinus      = 15,
        NuTau         = 16,
        NuTauBar      = -16,

        /* Nuclei */
        He3Nucleus  = 1000020030,
        He4Nucleus  = 1000020040,
        Li6Nucleus  = 1000030060,
        Li7Nucleus  = 1000030070,
        Be9Nucleus  = 1000040090,
        B10Nucleus  = 1000050100,
        B11Nucleus  = 1000050110,
        C12Nucleus  = 1000060120,
        C13Nucleus  = 1000060130,
        N14Nucleus  = 1000070140,
        N15Nucleus  = 1000070150,
        O16Nucleus  = 1000080160,
        O17Nucleus  = 1000080170,
        O18Nucleus  = 1000080180,
        F19Nucleus  = 1000090190,
        Ne20Nucleus = 1000100200,
        Ne21Nucleus = 1000100210,
        Ne22Nucleus = 1000100220,
        Na23Nucleus = 1000110230,
        Mg24Nucleus = 1000120240,
        Mg25Nucleus = 1000120250,
        Mg26Nucleus = 1000120260,
        Al26Nucleus = 1000130260,
        Al27Nucleus = 1000130270,
        Si28Nucleus = 1000140280,
        Si29Nucleus = 1000140290,
        Si30Nucleus = 1000140300,
        Si31Nucleus = 1000140310,
        Si32Nucleus = 1000140320,
        P31Nucleus  = 1000150310,
        P32Nucleus  = 1000150320,
        P33Nucleus  = 1000150330,
        S32Nucleus  = 1000160320,
        S33Nucleus  = 1000160330,
        S34Nucleus  = 1000160340,
        S35Nucleus  = 1000160350,
        S36Nucleus  = 1000160360,
        Cl35Nucleus = 1000170350,
        Cl36Nucleus = 1000170360,
        Cl37Nucleus = 1000170370,
        Ar36Nucleus = 1000180360,
        Ar37Nucleus = 1000180370,
        Ar38Nucleus = 1000180380,
        Ar39Nucleus = 1000180390,
        Ar40Nucleus = 1000180400,
        Ar41Nucleus = 1000180410,
        Ar42Nucleus = 1000180420,
        K39Nucleus  = 1000190390,
        K40Nucleus  = 1000190400,
        K41Nucleus  = 1000190410,
        Ca40Nucleus = 1000200400,
        Ca41Nucleus = 1000200410,
        Ca42Nucleus = 1000200420,
        Ca43Nucleus = 1000200430,
        Ca44Nucleus = 1000200440,
        Ca45Nucleus = 1000200450,
        Ca46Nucleus = 1000200460,
        Ca47Nucleus = 1000200470,
        Ca48Nucleus = 1000200480,
        Sc44Nucleus = 1000210440,
        Sc45Nucleus = 1000210450,
        Sc46Nucleus = 1000210460,
        Sc47Nucleus = 1000210470,
        Sc48Nucleus = 1000210480,
        Ti44Nucleus = 1000220440,
        Ti45Nucleus = 1000220450,
        Ti46Nucleus = 1000220460,
        Ti47Nucleus = 1000220470,
        Ti48Nucleus = 1000220480,
        Ti49Nucleus = 1000220490,
        Ti50Nucleus = 1000220500,
        V48Nucleus  = 1000230480,
        V49Nucleus  = 1000230490,
        V50Nucleus  = 1000230500,
        V51Nucleus  = 1000230510,
        Cr50Nucleus = 1000240500,
        Cr51Nucleus = 1000240510,
        Cr52Nucleus = 1000240520,
        Cr53Nucleus = 1000240530,
        Cr54Nucleus = 1000240540,
        Mn52Nucleus = 1000250520,
        Mn53Nucleus = 1000250530,
        Mn54Nucleus = 1000250540,
        Mn55Nucleus = 1000250550,
        Fe54Nucleus = 1000260540,
        Fe55Nucleus = 1000260550,
        Fe56Nucleus = 1000260560,
        Fe57Nucleus = 1000260570,
        Fe58Nucleus = 1000260580,

        /* The following are fake particles used in Icetray and have no official codes */
        /* The section abs(code) > 2000000000 is reserved for this kind of use */
        CherenkovPhoton      = 2000009900,
        Nu                   = -2000000004,
        Monopole             = -2000000041,
        Brems                = -2000001001,
        DeltaE               = -2000001002,
        PairProd             = -2000001003,
        NuclInt              = -2000001004,
        MuPair               = -2000001005,
        Hadrons              = -2000001006,
        ContinuousEnergyLoss = -2000001111,
        FiberLaser           = -2000002100,
        N2Laser              = -2000002101,
        YAGLaser             = -2000002201,
        STauPlus             = -2000009131,
        STauMinus            = -2000009132,
        SMP                  = -2000009500,
    };

    enum ParticleShape
    {
        Null           = 0,
        Primary        = 10,
        TopShower      = 20,
        Cascade        = 30,
        CascadeSegment = 31,
        InfiniteTrack  = 40,
        StartingTrack  = 50,
        StoppingTrack  = 60,
        ContainedTrack = 70,
        MCTrack        = 80,
        Dark           = 90
    };

    enum FitStatus
    {
        NotSet              = -1,
        OK                  = 0,
        GeneralFailure      = 10,
        InsufficientHits    = 20,
        FailedToConverge    = 30,
        MissingSeed         = 40,
        InsufficientQuality = 50
    };

    enum LocationType
    {
        Anywhere       = 0,
        IceTop         = 10,
        InIce          = 20,
        InActiveVolume = 30
    };

public:
    /** @brief Constructor for a simple particle with generated ID
     * useful in applications where the shape and type of an particle are in focus
     * @param shape Shape of the track
     * @param type Particle type
     */
    I3Particle(ParticleShape shape = Null, ParticleType type = unknown);

    /** @brief Constructor for I3Particle, focusing on directional properties
     * @param pos Position of the vertex
     * @param dir Direction of the track
     * @param vertextime time that the vertex is happening
     * @param shape Shape of the track
     * @param length
     */
#ifndef __CINT__
    I3Particle(const I3Position pos,
               const I3Direction dir,
               const double vertextime,
               ParticleShape shape = Null,
               double length       = NAN);
#else
    I3Particle(const I3Position pos,
               const I3Direction dir,
               const double vertextime,
               ParticleShape shape,
               double length);
#endif

#ifndef __CINT__
    /** @brief Constructor for particle from boost::optional<I3Particle>
     * @param p A particle
     */
    I3Particle(const boost::optional<I3Particle>& p);
#endif

    ~I3Particle();

    bool IsTrack() const;
    bool IsCascade() const;
    bool IsTopShower() const;
    bool IsNeutrino() const;
    bool IsNucleus() const;

    bool HasPosition() const;
    bool HasDirection() const;
    bool HasEnergy() const;

    operator I3ParticleID() const { return ID_; }
    I3Particle& operator=(const I3Particle& rhs)
    {
        ID_           = rhs.ID_;
        pdgEncoding_  = rhs.pdgEncoding_;
        shape_        = rhs.shape_;
        status_       = rhs.status_;
        pos_          = rhs.pos_;
        dir_          = rhs.dir_;
        time_         = rhs.time_;
        energy_       = rhs.energy_;
        length_       = rhs.length_;
        speed_        = rhs.speed_;
        locationType_ = rhs.locationType_;
        return *this;
    }
    bool operator==(const I3Particle& rhs) const
    {
        return (ID_ == rhs.ID_ && pdgEncoding_ == rhs.pdgEncoding_ && shape_ == rhs.shape_ && status_ == rhs.status_ &&
                pos_ == rhs.pos_ && dir_ == rhs.dir_ && CompareFloatingPoint::Compare_NanEqual(time_, rhs.time_) &&
                CompareFloatingPoint::Compare_NanEqual(energy_, rhs.energy_) &&
                CompareFloatingPoint::Compare_NanEqual(length_, rhs.length_) &&
                CompareFloatingPoint::Compare_NanEqual(speed_, rhs.speed_) && locationType_ == rhs.locationType_);
    }
    bool operator!=(const I3Particle& rhs) const { return !(*this == rhs); }

    /**
     * Returns a particle with a new I3ParticleID, but will all other properties unchanged.
     */
    I3Particle Clone() const;

    I3ParticleID GetID() const { return ID_; }
    int32_t GetMinorID() const { return ID_.minorID; }
    uint64_t GetMajorID() const { return ID_.majorID; }

    int32_t GetPdgEncoding() const { return pdgEncoding_; }
    void SetPdgEncoding(int32_t newid) { pdgEncoding_ = newid; }

    ParticleType GetType() const { return ParticleType(pdgEncoding_); };
    std::string GetTypeString() const;

    void SetType(ParticleType type) { pdgEncoding_ = type; };
    void SetTypeString(const std::string& str);

    ParticleShape GetShape() const { return shape_; }
    std::string GetShapeString() const;
    void SetShape(ParticleShape shape) { shape_ = shape; }
    void SetShapeString(const std::string& str);

    FitStatus GetFitStatus() const { return status_; }
    std::string GetFitStatusString() const;
    void SetFitStatus(FitStatus status) { status_ = status; }
    void SetFitStatusString(const std::string& str);

    LocationType GetLocationType() const { return locationType_; }
    std::string GetLocationTypeString() const;
    void SetLocationType(LocationType type) { locationType_ = type; }
    void SetLocationTypeString(const std::string& str);

    const I3Position& GetPos() const { return pos_; }
    void SetPos(const I3Position& p) { pos_ = p; }
    void SetPos(double p1, double p2, double p3, I3Position::RefFrame frame) { pos_ = I3Position(p1, p2, p3, frame); }
    void SetPos(double x, double y, double z) { pos_ = I3Position(x, y, z); }

    const I3Direction& GetDir() const { return dir_; }
    void SetDir(const I3Direction& d) { dir_ = d; }
    void SetDir(double zen, double azi) { dir_ = I3Direction(zen, azi); }
    void SetDir(double x, double y, double z) { dir_ = I3Direction(x, y, z); }

    void SetThetaPhi(double theta, double phi) { dir_.SetThetaPhi(theta, phi); }

    // handy shortcuts for components
    double GetZenith() const { return dir_.GetZenith(); }
    double GetAzimuth() const { return dir_.GetAzimuth(); }
    double GetX() const { return pos_.GetX(); }
    double GetY() const { return pos_.GetY(); }
    double GetZ() const { return pos_.GetZ(); }

    double GetTime() const { return time_; }
    void SetTime(double t) { time_ = t; }

    double GetLength() const { return length_; }
    void SetLength(double length) { length_ = length; }

    /** @brief Returns the kinetic energy of the particle. */
    double GetEnergy() const { return energy_; }

    /** @brief Returns the kinetic energy of the particle. */
    double GetKineticEnergy() const { return GetEnergy(); }

    /** @brief Returns the total energy of the particle. */
    double GetTotalEnergy() const;

    /** @brief Returns the mass of the particle. */
    double GetMass() const;

    /** @brief Returns the mass of a particle of given type. */
    static double GetMassForType(ParticleType type);

    /** @brief Check if particle has mass. */
    bool HasMass() const;

    /** @brief Sets the kinetic energy of the particle. */
    void SetEnergy(double energy) { energy_ = energy; }

    /** @brief Sets the kinetic energy of the particle. */
    void SetKineticEnergy(double energy) { SetEnergy(energy); }

    /** @brief Sets the total energy of the particle. */
    void SetTotalEnergy(double total_energy);

    double GetSpeed() const { return speed_; }
    void SetSpeed(double s) { speed_ = s; }

    /** @brief get the position of the particle at this time
     *  @note ignores the shape of the track (start/stopping point), so the particle might not be definded at that very
     * position
     *  @param time the time in ns
     */
    I3Position ShiftTimeTrack(const double time) const;

    I3Position GetStartPos() const;

    double GetStartTime() const;

    I3Position GetStopPos() const;

    double GetStopTime() const;

private:
    I3ParticleID ID_;
    int32_t pdgEncoding_;
    ParticleShape shape_;
    FitStatus status_;
    I3Position pos_;
    I3Direction dir_;
    double time_;
    double energy_;
    double length_;
    double speed_;
    LocationType locationType_;

    /** @brief Constructor for particle as a unique identifier by majorID and minorID
     * useful in iteration processes (in I3MCTree) where unique tagging is required
     * @param major MajorID of that particle
     * @param minor MinorID of that particle
     */
    I3Particle(const uint64_t major, const int32_t minor);

    /** @brief generate a new, unique ID combination
     */
    void generateID();

    friend class I3Stochastic;
    // since that ctor is private we need to give this
    // class friend access to test it
    friend class test_particle_id_private_ctor;

    friend class icecube::serialization::access;
    template<class Archive>
    void save(Archive& ar, unsigned version) const;
    template<class Archive>
    void load(Archive& ar, unsigned version);
    I3_SERIALIZATION_SPLIT_MEMBER();
};

static const unsigned i3particle_version_ = 5;

// let other code know that I3Particle stores PDG encodings internally
#define I3PARTICLE_SUPPORTS_PDG_ENCODINGS

/**
 * List the names of enumeration members defined in this file
 * here. These can be used for e.g. pybindings, which require
 * the names of the enumeration members to be known. This list
 * should be updated whenver members or new enums are added to
 * the class.
 */
#define I3PARTICLE_H_I3Particle_ParticleType                                                                           \
    (unknown)(Gamma)(EPlus)(EMinus)(MuPlus)(MuMinus)(Pi0)(PiPlus)(PiMinus)(K0_Long)(KPlus)(KMinus)(Neutron)(PPlus)(    \
        PMinus)(K0_Short)(Eta)(Lambda)(SigmaPlus)(Sigma0)(SigmaMinus)(Xi0)(XiMinus)(OmegaMinus)(NeutronBar)(           \
        LambdaBar)(SigmaMinusBar)(Sigma0Bar)(SigmaPlusBar)(Xi0Bar)(XiPlusBar)(OmegaPlusBar)(DPlus)(DMinus)(D0)(D0Bar)( \
        DsPlus)(DsMinusBar)(LambdacPlus)(WPlus)(WMinus)(Z0)(NuE)(NuEBar)(NuMu)(NuMuBar)(TauPlus)(TauMinus)(NuTau)(     \
        NuTauBar)(He3Nucleus)(He4Nucleus)(Li6Nucleus)(Li7Nucleus)(Be9Nucleus)(B10Nucleus)(B11Nucleus)(C12Nucleus)(     \
        C13Nucleus)(N14Nucleus)(N15Nucleus)(O16Nucleus)(O17Nucleus)(O18Nucleus)(F19Nucleus)(Ne20Nucleus)(Ne21Nucleus)( \
        Ne22Nucleus)(Na23Nucleus)(Mg24Nucleus)(Mg25Nucleus)(Mg26Nucleus)(Al26Nucleus)(Al27Nucleus)(Si28Nucleus)(       \
        Si29Nucleus)(Si30Nucleus)(Si31Nucleus)(Si32Nucleus)(P31Nucleus)(P32Nucleus)(P33Nucleus)(S32Nucleus)(           \
        S33Nucleus)(S34Nucleus)(S35Nucleus)(S36Nucleus)(Cl35Nucleus)(Cl36Nucleus)(Cl37Nucleus)(Ar36Nucleus)(           \
        Ar37Nucleus)(Ar38Nucleus)(Ar39Nucleus)(Ar40Nucleus)(Ar41Nucleus)(Ar42Nucleus)(K39Nucleus)(K40Nucleus)(         \
        K41Nucleus)(Ca40Nucleus)(Ca41Nucleus)(Ca42Nucleus)(Ca43Nucleus)(Ca44Nucleus)(Ca45Nucleus)(Ca46Nucleus)(        \
        Ca47Nucleus)(Ca48Nucleus)(Sc44Nucleus)(Sc45Nucleus)(Sc46Nucleus)(Sc47Nucleus)(Sc48Nucleus)(Ti44Nucleus)(       \
        Ti45Nucleus)(Ti46Nucleus)(Ti47Nucleus)(Ti48Nucleus)(Ti49Nucleus)(Ti50Nucleus)(V48Nucleus)(V49Nucleus)(         \
        V50Nucleus)(V51Nucleus)(Cr50Nucleus)(Cr51Nucleus)(Cr52Nucleus)(Cr53Nucleus)(Cr54Nucleus)(Mn52Nucleus)(         \
        Mn53Nucleus)(Mn54Nucleus)(Mn55Nucleus)(Fe54Nucleus)(Fe55Nucleus)(Fe56Nucleus)(Fe57Nucleus)(Fe58Nucleus)(       \
        CherenkovPhoton)(Nu)(Monopole)(Brems)(DeltaE)(PairProd)(NuclInt)(MuPair)(Hadrons)(ContinuousEnergyLoss)(       \
        FiberLaser)(N2Laser)(YAGLaser)(STauPlus)(STauMinus)(SMP)

#define I3PARTICLE_H_I3Particle_ParticleShape                                                                          \
    (Null)(Primary)(TopShower)(Cascade)(CascadeSegment)(InfiniteTrack)(StartingTrack)(StoppingTrack)(ContainedTrack)(  \
        MCTrack)(Dark)

#define I3PARTICLE_H_I3Particle_FitStatus                                                                              \
    (NotSet)(OK)(GeneralFailure)(InsufficientHits)(FailedToConverge)(MissingSeed)(InsufficientQuality)

#define I3PARTICLE_H_I3Particle_LocationType (Anywhere)(IceTop)(InIce)(InActiveVolume)

#ifndef __CINT__
// template specialization for XML i/o
template<>
void I3Particle::save(icecube::archive::xml_oarchive& ar, unsigned version) const;
template<>
void I3Particle::load(icecube::archive::xml_iarchive& ar, unsigned version);
#endif

std::ostream& operator<<(std::ostream& oss, const I3Particle& d);

I3_POINTER_TYPEDEFS(I3Particle);
I3_CLASS_VERSION(I3Particle, i3particle_version_);

typedef I3Vector<I3Particle> I3VectorI3Particle;
I3_POINTER_TYPEDEFS(I3VectorI3Particle);

#endif
