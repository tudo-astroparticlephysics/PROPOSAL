
#pragma once

#include <string>
#include <vector>

namespace PROPOSAL
{

class Vector3D;
class DynamicData;
class Particle;

class DecayChannel
{

    public:
    typedef std::vector<Particle*> DecayProducts;

    DecayChannel();
    DecayChannel(const DecayChannel&);
    virtual ~DecayChannel();

    virtual DecayChannel* clone() const = 0;

    bool operator==(const DecayChannel&) const;
    bool operator!=(const DecayChannel&) const;

    friend std::ostream& operator<<(std::ostream&, DecayChannel const&);

    // --------------------------------------------------------------------- //
    // Public methods
    // --------------------------------------------------------------------- //

    virtual DecayProducts Decay(Particle&) = 0;

    // ----------------------------------------------------------------------------
    /// @brief Boost the particle along a direction
    ///
    /// Internally the energy / momentum and direction of the particle
    /// will be set according to the boost.
    ///
    /// @param Particle
    /// @param direction
    /// @param beta = v/c
    // ----------------------------------------------------------------------------
    static void Boost(Particle&, const Vector3D& direction, double beta);

    // ----------------------------------------------------------------------------
    /// @brief Boost a set of particles along a direction
    ///
    /// Same as Boost(Particle&, const Vector3D&, double) but with
    /// a std::vector of Particles.
    ///
    /// @param DecayProducts
    /// @param direction
    /// @param beta = v/c
    // ----------------------------------------------------------------------------
    static void Boost(DecayProducts&, const Vector3D& direction, double beta);

    // ----------------------------------------------------------------------------
    /// @brief Calculate the momentum in a two-body-phase-space decay
    ///
    // Returns:
    //
    // $P = \frac{\sqrt{\lambda(m_1, m_2, m_3)}}{2 m_1}$
    //
    /// @param m1
    /// @param m2
    /// @param m3
    ///
    /// @return
    // ----------------------------------------------------------------------------
    static double Momentum(double m1, double m2, double m3);

    // ----------------------------------------------------------------------------
    /// @brief Set common particle properties to the decay products
    ///
    /// Sets:
    ///  - Position
    ///  - Time
    ///  - ParentParticle energy
    ///  - Particle Id
    ///  - ParentParticle Id
    ///
    /// @return
    // ----------------------------------------------------------------------------
    static Vector3D GenerateRandomDirection();

    virtual const std::string& GetName() const = 0;

    protected:
    virtual bool compare(const DecayChannel&) const = 0;
    virtual void print(std::ostream&) const {};

    static void CopyParticleProperties(DecayProducts&, const Particle&);
};

} /* PROPOSAL */
