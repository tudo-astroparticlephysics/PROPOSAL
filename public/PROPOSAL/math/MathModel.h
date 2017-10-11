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

// class MathModel
// {
//
// public:
//
// //----------------------------------------------------------------------------//
//     //Constructor
//
//     MathModel();
//     MathModel(const MathModel&);
//     MathModel& operator=(const MathModel&);
//
// //----------------------------------------------------------------------------//
//     //Destructor
//
//     virtual ~MathModel();
//
// //----------------------------------------------------------------------------//
//     //Memberfunctions
//
//     double RandomDouble();
//
// //----------------------------------------------------------------------------//
//
//     static void set_seed(int seed);
//
//     #<{(|* @brief Set a custom random number generator
//      *
//      * Classes that contain other subclasses of MathModel should
//      * override this to pass the new RNG on to their members.
//      |)}>#
//     virtual void SetRandomNumberGenerator(boost::function<double ()> &f);
//
// private:
//         static double DefaultRandomDouble();
//
// 	boost::function<double ()> rng_;
// 	static boost::mt19937 *default_rng_;
// };

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

    double RandomDouble();

    void SetSeed(int seed);

    /** @brief Set a custom random number generator
     *
     * Classes that contain other subclasses of MathModel should
     * override this to pass the new RNG on to their members.
     */
    virtual void SetRandomNumberGenerator(boost::function<double ()> &f);

    private:
    RandomGenerator();
    virtual ~RandomGenerator();

    static double DefaultRandomDouble();

    static boost::random::mt19937 rng_;
    static boost::variate_generator<boost::mt19937&, boost::uniform_real<> > variate_real;
	boost::function<double ()> random_function;
};


} // namespace PROPOSAL

#endif // MATHMODEL_H_
