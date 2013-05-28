/*! \file   ContinuousRandomization.h
*   \brief  Headerfile for the randomization of continous energy losses.
*
*   For more details see the class documentation.
*
*   \date   2013.05.28
*   \author Jan-Hednrik Köhne
*/

#ifndef CONTINUOUSRANDOMIZATION_H_
#define CONTINUOUSRANDOMIZATION_H_


/**
 * \brief Class containing the functions to randomize the continuous energy losses
 *
 */


class ContinuousRandomization
{

private:

public:

 //----------------------------------------------------------------------------//
    // Constructors

    ContinuousRandomization();
    ContinuousRandomization(const ContinuousRandomization&);
    ContinuousRandomization& operator=(const ContinuousRandomization& continuous_randomization);
    bool operator==(const ContinuousRandomization &continuous_randomization) const;
    bool operator!=(const ContinuousRandomization &continuous_randomization) const;

//----------------------------------------------------------------------------//

    //Memberfunction


//----------------------------------------------------------------------------//

    void swap(ContinuousRandomization &continuous_randomization);

//----------------------------------------------------------------------------//
    //Getter

//----------------------------------------------------------------------------//
    //Setter
};

#endif /* CONTINUOUSRANDOMIZATION_H_ */
