#include "PROPOSAL/density_distr/density_distr.h"
                                                    
Density_distr::Density_distr(double rho, std::function<double(double)> density_distribution):                                
    rho_(rho),                               
    density_distribution_(density_distribution)                                
{
}
                                                                                              
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                             
// %%%%%%%%%%%%%%%%%%%      Setter      %%%%%%%%%%%%%%%%%%%%%%                                
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                
                                                                                                 
void Density_distr::SetXi(double xi)                                
{                                                                        
    xi_ = xi;                                                                                     
}
                                 
                                                                                              
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                               
// %%%%%%%%%%%%%%%%%%%      Getter      %%%%%%%%%%%%%%%%%%%%%%                                                           
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                            
                                                                                                                                              
std::function<double(double)> Density_distr::GetDensityDistribution()
{                                                                
    return density_distribution_;                                
}                                
                                                                                                                                                    
double Density_distr::GetXi()                                                                                                          
{                                                        
    return xi_;                                                                            
}                                
                                                                
double Density_distr::GetRho()                                
{                                                                
    return rho_;                                
}
