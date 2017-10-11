/*! \file   MathModel.cxx
*   \brief  Source file for the used MathModels in PROPOSAL.
*
*   For more details see the class documentation.

*   \author Jan-Hendrik Koehne
*/


// #include <iostream>

// #include <boost/generator_iterator.hpp>

#include "PROPOSAL/math/MathModel.h"

using namespace PROPOSAL;

boost::random::mt19937 RandomGenerator::rng_;
boost::variate_generator<boost::mt19937&, boost::uniform_real<> > RandomGenerator::variate_real(RandomGenerator::rng_, boost::uniform_real<>(0.0, 1.0));


RandomGenerator::RandomGenerator()
    : random_function(&RandomGenerator::DefaultRandomDouble)
{
}

RandomGenerator::~RandomGenerator()
{
}

double RandomGenerator::RandomDouble()
{
    return random_function();
}

void RandomGenerator::SetSeed(int seed)
{
    boost::mt19937::result_type s = static_cast<boost::mt19937::result_type>(seed);
    rng_.seed(s);
}

void RandomGenerator::SetRandomNumberGenerator(boost::function<double ()> &f)
{
    random_function = f;
}

double RandomGenerator::DefaultRandomDouble()
{
    return variate_real();
}
