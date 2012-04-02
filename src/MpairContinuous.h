#ifndef MPAIRCONTINUOUS_H
#define MPAIRCONTINUOUS_H
#include"Mpairproduction.h"

class MpairContinuous : public Mpairproduction{


protected:

    Interpolate *interpolateJ_;

    // funtion to integrate/interpolate

    double function(double v);

    double functionInt(double e);

public:

   // usual constructor

    MpairContinuous(Mpairproduction *MPAIR);

    // This calculates the contribution to dEdx() up to v_max.
    double dEdx();



    void activate_interpolation();



};

#endif // MPAIRCONTINUOUS_H
