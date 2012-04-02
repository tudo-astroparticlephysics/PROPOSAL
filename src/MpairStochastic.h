#ifndef MPAIRSTOCHASTIC_H
#define MPAIRSTOCHASTIC_H

#include "Mpairproduction.h"

class MpairStochastic : public Mpairproduction {

protected:

    std::vector<Integral*> integral_;
    std::vector<double> H_;
    std::vector<Interpolate*>  interpolateJ_;
    std::vector<Interpolate*>  interpolateJo_;
    double sum;
    double rnd;

    double functionInt(double v);
    double functionInt(double e, double v);
    double function(double v);

public:

    MpairStochastic(Mpairproduction *MPAIR);
    double dNdx();
    double dNdx(double rnd);
    double e(double rnd);
    void activate_interpolation();


};

#endif // MPAIRSTOCHASTIC_H
