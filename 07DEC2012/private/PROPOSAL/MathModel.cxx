/*! \file   MathModel.cxx
*   \brief  Source file for the used MathModels in PROPOSAL.
*
*   For more details see the class documentation.

*   \author Jan-Hendrik Koehne
*/


#include "PROPOSAL/MathModel.h"
#include "iostream"

RNGType* MathModel::rng = new RNGType();

//----------------------------------------------------------------------------//

MathModel::MathModel()
{

}

//----------------------------------------------------------------------------//

void MathModel::set_seed(int seed)
{
    boost::mt19937::result_type s = static_cast<boost::mt19937::result_type>(seed);
    rng->seed(s);
    //rng->reset()

}

//----------------------------------------------------------------------------//

double MathModel::RandomDouble(){


    boost::uniform_real<> dist(0, 1);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> > die(*rng, dist);
    double x = die();
    return x;
}

//----------------------------------------------------------------------------//

//double MathModel::RandomDouble()

//{

// // boost::mt19937 rng(43);
//  static boost::uniform_01<boost::mt19937> zeroone(*rng);
//  return zeroone();
//}
