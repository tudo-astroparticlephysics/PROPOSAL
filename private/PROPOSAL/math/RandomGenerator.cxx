/*! \file   MathModel.cxx
*   \brief  Source file for the used MathModels in PROPOSAL.
*
*   For more details see the class documentation.

*   \author Jan-Hendrik Koehne
*/

#include "PROPOSAL/math/RandomGenerator.h"

using namespace PROPOSAL;

boost::random::mt19937 RandomGenerator::rng_;
boost::variate_generator<boost::mt19937&, boost::uniform_real<> > RandomGenerator::variate_real(
    RandomGenerator::rng_,
    boost::uniform_real<>(0.0, 1.0));

// ------------------------------------------------------------------------- //
// Constructor & destructor
// ------------------------------------------------------------------------- //

RandomGenerator::RandomGenerator()
    : random_function(&RandomGenerator::DefaultRandomDouble)
{
}

RandomGenerator::~RandomGenerator() {}

// ------------------------------------------------------------------------- //
// Methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double RandomGenerator::RandomDouble()
{
    return random_function();
}

// ------------------------------------------------------------------------- //
void RandomGenerator::SetSeed(int seed)
{
    boost::mt19937::result_type s = static_cast<boost::mt19937::result_type>(seed);
    rng_.seed(s);
}

// ------------------------------------------------------------------------- //
void RandomGenerator::Serialize(std::ostream& os)
{
    os << rng_;
}

// ------------------------------------------------------------------------- //
void RandomGenerator::Deserialize(std::istream& is)
{
    is >> rng_;
}

// ------------------------------------------------------------------------- //
void RandomGenerator::SetRandomNumberGenerator(boost::function<double()>& f)
{
    random_function = f;
}

// ------------------------------------------------------------------------- //
double RandomGenerator::DefaultRandomDouble()
{
    return variate_real();
}
