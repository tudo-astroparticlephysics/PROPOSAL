/*! \file   MathModel.cxx
*   \brief  Source file for the used MathModels in PROPOSAL.
*
*   For more details see the class documentation.

*   \author Jan-Hendrik Koehne
*/

#include "PROPOSAL/math/RandomGenerator.h"

using namespace PROPOSAL;

std::mt19937 RandomGenerator::rng_;
std::uniform_real_distribution<double> RandomGenerator::uniform_distribution(0.0, 1.0);

// ------------------------------------------------------------------------- //
// Constructor & destructor
// ------------------------------------------------------------------------- //

RandomGenerator::RandomGenerator()
    : random_function(&RandomGenerator::DefaultRandomDouble)
#ifdef ICECUBE_PROJECT
    , i3random_gen_(NULL)
#endif
{
}

RandomGenerator::~RandomGenerator() {}

// ------------------------------------------------------------------------- //
// Methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double RandomGenerator::RandomDouble()
{
#ifdef ICECUBE_PROJECT
    if (i3random_gen_)
    {
        return i3random_gen_->Uniform(0.0, 1.0);
    }
    else
    {
        return random_function();
    }
#else
    return random_function();
#endif
}

// ------------------------------------------------------------------------- //
void RandomGenerator::SetSeed(int seed)
{
    rng_.seed(seed);
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
void RandomGenerator::SetRandomNumberGenerator(std::function<double()>& f)
{
    random_function = f;
}

#ifdef ICECUBE_PROJECT
void RandomGenerator::SetI3RandomNumberGenerator(I3RandomServicePtr random)
{
    i3random_gen_ = random.get();
}
#endif

void RandomGenerator::SetDefaultRandomNumberGenerator()
{
    random_function = &RandomGenerator::DefaultRandomDouble;
}

// ------------------------------------------------------------------------- //
double RandomGenerator::DefaultRandomDouble()
{
    // return variate_real();
    return uniform_distribution(rng_);
}
