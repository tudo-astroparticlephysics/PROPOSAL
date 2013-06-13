/*! \file   Scattering.h
*   \brief  Header file for the Scattering routines.
*
*   For more details see the class documentation.
*
*   \date   2013.06.13
*   \author Tomasz Fuchs
*/


#ifndef SCATTERING_H
#define SCATTERING_H
#include "vector"
#include <string>

/**
  * \brief This class provides the scattering routine provided by moliere.
  *
  * More precise scattering angles will be added soon.
  */


class Scattering
{


private:


//----------------------------------------------------------------------------//

public:

    /**
     * \brief Default Constructor
     *
     * Constructor which sets "default" settings.
     */
    Scattering();

//----------------------------------------------------------------------------//

    Scattering(const Scattering&);
    Scattering& operator=(const Scattering&);
    bool operator==(const Scattering &scattering) const;
    bool operator!=(const Scattering &scattering) const;
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
    // Memberfunctions


//----------------------------------------------------------------------------//

    void swap(Scattering &scattering);

//----------------------------------------------------------------------------//

    //Setter



//----------------------------------------------------------------------------//
    // Getter


//----------------------------------------------------------------------------//
    // destructors
    ~Scattering() {}


};



#endif //SCATTERING_H
