#include "PROPOSAL/density_distr/density_distr.h"
#include "PROPOSAL/math/Vector3D.h"
                                                    
Density_distr::Density_distr(Vector3D fAxis, Vector3D fp0, std::function<double(double)> density_distribution):
    fAxis_(fAxis),                               
    fp0_(fp0),
    density_distribution_(density_distribution)
{
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
