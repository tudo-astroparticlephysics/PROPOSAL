#include "PROPOSAL/density_distr/density_distr.h"
#include "PROPOSAL/math/Vector3D.h"
#include <iostream>

Density_distr::Density_distr()
{
    Vector3D fAxis_(0,0,0);
    Vector3D fp0_(0,0,0);
}

Density_distr::Density_distr(const Density_distr& density_distr):
    fAxis_(density_distr.fAxis_),
    fp0_(density_distr.fp0_)
{}
                                                    
Density_distr::Density_distr(Vector3D fAxis, Vector3D fp0):
    fAxis_(fAxis),                               
    fp0_(fp0)
{}

double Density_distr::GetDepth(Vector3D xi) const
{
    return GetAxis()* (xi - GetFp0());
}
                                                                                              
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%      Getter      %%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

std::function<double(double)> Density_distr::GetDensityDistribution()
{                                                                
    return density_distribution_;                                
}                                
 

