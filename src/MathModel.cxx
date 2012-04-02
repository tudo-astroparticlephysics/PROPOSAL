#include "MathModel.h"
#include "iostream"
#include "global.h"
RNGType* MathModel::rng = new RNGType();
//boost::mt19937* MathModel::rng = new boost::mt19937(4);



MathModel::MathModel(){


}

//double MathModel::RandomDouble(){

// double x;
// randomIn>>x;
////    boost::uniform_real<> dist(0, 1);
////    boost::variate_generator<boost::mt19937&, boost::uniform_real<> > die(*rng, dist);
////    double x = die();
// //std::cout<<"x="<<x<<std::endl;
//    return x;
//  //  return 0.7;
//}

void MathModel::set_seed(int seed){

    boost::mt19937::result_type s = static_cast<boost::mt19937::result_type>(seed);
    rng->seed(s);
    //rng->reset()

}

double MathModel::RandomDouble(){


    boost::uniform_real<> dist(0, 1);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> > die(*rng, dist);
    double x = die();
    return x;
}

//double MathModel::RandomDouble()

//{

// // boost::mt19937 rng(43);
//  static boost::uniform_01<boost::mt19937> zeroone(*rng);
//  return zeroone();
//}
