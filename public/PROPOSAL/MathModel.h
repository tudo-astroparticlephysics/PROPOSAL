/**
 * root of pure math functions
 */
#ifndef MATHMODEL_H_
#define MATHMODEL_H_

#include "iostream"
#include "fstream"
#include "boost/function.hpp"


#include "boost/random.hpp"

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

#endif // MATHMODEL_H_ 
