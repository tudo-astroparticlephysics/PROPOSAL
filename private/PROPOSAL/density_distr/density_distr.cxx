#include "PROPOSAL/density_distr/density_distr.h"
#include "PROPOSAL/math/Vector3D.h"

Density_distr::Density_distr()
{
    Vector3D fAxis_(0,0,0);
    Vector3D fp0_(0,0,0);
    std::function<double(double)> density_distribution_ = nullptr;
}

Density_distr::Density_distr(const Density_distr& density_distr)
{
    Vector3D fAxis_(density_distr.fAxis_);
    Vector3D fp0_(density_distr.fp0_);
    std::function<double(double)> density_distribution_ = density_distr.density_distribution_;
}
                                                    
Density_distr::Density_distr(Vector3D fAxis, Vector3D fp0, std::function<double(double)> density_distribution):
    fAxis_(fAxis),                               
    fp0_(fp0),
    density_distribution_(density_distribution)
{
}

double Density_distr::GetDepth(Vector3D xi) const
{
    return fAxis_ * (xi - fp0_);
}
                                                                                              
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%      Getter      %%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

std::function<double(double)> Density_distr::GetDensityDistribution()
{                                                                
    return density_distribution_;                                
}                                
                                                                                                                                                    
Vector3D Density_distr::GetAxis() 
{ 
    return fAxis_;
}
                                                                
Vector3D Density_distr::GetFp0()                                
{                                                                
    return fp0_;                                
}
