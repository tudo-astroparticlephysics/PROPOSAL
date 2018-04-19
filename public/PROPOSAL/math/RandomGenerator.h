/**
 * root of pure math functions
 */
#pragma once

#ifndef RANDOMGENERATOR_H_
#define RANDOMGENERATOR_H_

#include <boost/function.hpp>
#include <boost/random.hpp>

namespace PROPOSAL {

// ----------------------------------------------------------------------------
/// @brief Random number generator
// ----------------------------------------------------------------------------
class RandomGenerator
{
public:
    static RandomGenerator& Get()
    {
        static RandomGenerator instance;
        return instance;
    }

    // ----------------------------------------------------------------------------
    /// @brief Execute the given rng to get a random number
    ///
    /// @return random number
    // ----------------------------------------------------------------------------
    double RandomDouble();

    void SetSeed(int seed);

    // ----------------------------------------------------------------------------
    /// @brief Serialize the rng to a stream
    ///
    /// Useful for debuging to save a specific state.
    /// Only supported for internal used rng from the boost libraries.
    ///
    /// @param std::ostream
    // ----------------------------------------------------------------------------
    void Serialize(std::ostream&);

    // ----------------------------------------------------------------------------
    /// @brief Deserialize the rng from a stream
    ///
    /// Useful for debuging to get back a specific state.
    /// Only supported for internal used rng from the boost libraries.
    ///
    /// @param std::ostream
    // ----------------------------------------------------------------------------
    void Deserialize(std::istream&);

    /** @brief Set a custom random number generator
     *
     * Classes that contain other subclasses of MathModel should
     * override this to pass the new RNG on to their members.
     */
    virtual void SetRandomNumberGenerator(boost::function<double()>& f);

private:
    RandomGenerator();
    virtual ~RandomGenerator();

    static double DefaultRandomDouble();

    static boost::random::mt19937 rng_;
    static boost::variate_generator<boost::mt19937&, boost::uniform_real<> > variate_real;
    boost::function<double()> random_function;
};

} // namespace PROPOSAL

#endif // RANDOMGENERATOR_H_
