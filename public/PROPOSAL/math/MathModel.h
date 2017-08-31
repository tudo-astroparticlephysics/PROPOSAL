/**
 * root of pure math functions
 */
#pragma once

#ifndef MATHMODEL_H_
#define MATHMODEL_H_

// #include <fstream>

#include <boost/function.hpp>
#include <boost/random.hpp>


namespace PROPOSAL{

class MathModel
{

public:

//----------------------------------------------------------------------------//
    //Constructor

    MathModel();
    MathModel(const MathModel&);
    MathModel& operator=(const MathModel&);

//----------------------------------------------------------------------------//
    //Destructor

    virtual ~MathModel();

//----------------------------------------------------------------------------//
    //Memberfunctions

    double RandomDouble();

//----------------------------------------------------------------------------//

    static void set_seed(int seed);

    /** @brief Set a custom random number generator
     *
     * Classes that contain other subclasses of MathModel should
     * override this to pass the new RNG on to their members.
     */
    virtual void SetRandomNumberGenerator(boost::function<double ()> &f);

private:
        static double DefaultRandomDouble();

	boost::function<double ()> rng_;
	static boost::mt19937 *default_rng_;
};

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

    double RandomDouble() { return var_real(); }
    void SetSeed(int seed)
    {
        boost::mt19937::result_type s = static_cast<boost::mt19937::result_type>(seed);
        rng_.seed(s);
    }

    private:
    RandomGenerator()
        : rng_()
        , var_real(rng_, boost::uniform_real<>(0.0, 1.0))
    {
    }
    virtual ~RandomGenerator() {}

    boost::random::mt19937 rng_;
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> > var_real;
};

} // namespace PROPOSAL

#endif // MATHMODEL_H_
