
#include <cmath>

#include "PROPOSAL/crossection/parametrization/PhotoPairProduction.h"

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"


using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

PhotoPairProduction::PhotoPairProduction(const ParticleDef& particle_def,
                                 const Medium& medium,
                                 double multiplier)
        : Parametrization(particle_def, medium, EnergyCutSettings(), multiplier)
{
}

PhotoPairProduction::PhotoPairProduction(const PhotoPairProduction& param)
        : Parametrization(param)
{
}

PhotoPairProduction::~PhotoPairProduction() {}

bool PhotoPairProduction::compare(const Parametrization& parametrization) const
{
    return Parametrization::compare(parametrization);
}

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

Parametrization::IntegralLimits PhotoPairProduction::GetIntegralLimits(double energy)
{
    IntegralLimits limits;

    // x is the integration variable here
    if(energy <= 2. * ME){
        limits.vMin = 0.5;
        limits.vUp = limits.vMin;
        limits.vMax = 0.5;
    }
    else {
        limits.vMin = ME / energy;
        limits.vUp = limits.vMin; //treat interaction as fully stochastic
        limits.vMax = 1 - ME / energy;
    }

    return limits;
}

size_t PhotoPairProduction::GetHash() const
{
    size_t seed = Parametrization::GetHash();
    hash_combine(seed);

    return seed;
}

// ------------------------------------------------------------------------- //
// Specific implementations
// ------------------------------------------------------------------------- //

PhotoPairTsai::PhotoPairTsai(const ParticleDef& particle_def,
                                                 const Medium& medium,
                                                 double multiplier)
        : PhotoPairProduction(particle_def, medium, multiplier)
{

}

PhotoPairTsai::PhotoPairTsai(const PhotoPairTsai& param)
        : PhotoPairProduction(param)
{

}

PhotoPairTsai::~PhotoPairTsai()
{

}

bool PhotoPairTsai::compare(const Parametrization& parametrization) const
{
    return PhotoPairProduction::compare(parametrization);
}


double PhotoPairTsai::DifferentialCrossSection(double energy, double x)
{
    // Pair production and bremsstrahlung of chraged leptons, Yung-Su Tsai,
    // Review of Modern Physics, Vol. 46, No. 4, October 1974
    // see formula (3.9)

    double Phi1, Phi2, Psi1, Psi2;
    double Z = components_[component_index_]->GetNucCharge();
    double logZ = std::log(components_[component_index_]->GetNucCharge());
    double eta;

    double k = energy;

    double delta = std::pow(ME, 2.) * k / (2. * energy );


    if(Z<=2.5){
        // Hydrogen or Helium
        if(Z<=1.5){
            eta = 1;
        }
        else{
            eta = 1.6875;
        }

        double C = delta / ( 2. * ALPHA * ME * eta );
        delta = std::pow(ME, 2.) / ( 2. * k * x * (1. - x)); // (3.20);

        Phi1 = 4./3. * logZ + 4. * std::log(1. / (2. * eta * ALPHA)) + 13./3. - 2. * std::log( 1. + std::pow(C, 2.))
                - (13./2.) * C * std::atan(1./C) + (1./6.) / (1 + std::pow(C, -2.)); // (3.25), with erratum
        Phi2 = 4./3. * logZ + 4. * std::log(1. / (2. * eta * ALPHA)) + 11./3. - 2. * std::log(1 + std::pow(C, 2.))
                + 25. * std::pow(C, 2.) * (1. - C * std::atan(1./C)) - 14. * std::pow(C, 2.) * std::log(1 + std::pow(C, -2.)); // (3.26)
        Psi1 = 8./3. * logZ + 4. * std::log(1. / (2. * eta * ALPHA)) + 23./3. - 2. * std::log(1 + std::pow(C, 2.))
                - 17.5 * C * std::atan(1./C) + 8 * std::pow(C, 2.) * std::log(1 + std::pow(C, -2.)) - 1./6. / (1 + std::pow(C, -2.)); // (3.27)
        Psi2 = 8./3. * logZ + 4. * std::log(1. / (2. * eta * ALPHA)) + 21./3. - 2. * std::log(1 + std::pow(C, 2.))
                - 105. * std::pow(C, 2.) * (1. - C * std::atan(1./C)) + 50 * std::pow(C, 2.) * std::log(1 + std::pow(C, -2.))
                - 24. * std::pow(C, 2.) * ( -std::log(std::pow(C, 2.)) * std::log(1 + std::pow(C, -2.)) + dilog(1 + std::pow(C, -2.)) - dilog(1) ); // (3.28)
    }else if (Z<=4.5){
        // Lithium or Berylium
        double a, b, ap, bp;
        if(Z<=3.5){
            a = 100. * std::pow(Z, -1./3.) / ME;
            ap = 418.6 * std::pow(Z, -2./3.) / ME;
        }
        else{
            a = 106. * std::pow(Z, -1./3.) / ME;
            ap = 571.4 * std::pow(Z, -2./3.) / ME;
        }
        b = delta * a;
        bp = delta * ap;

        Phi1 = 2. * (1. + std::log(std::pow(a, 2.) * std::pow(Z, 2./3.) * std::pow(ME, 2.))) - 2. * std::log(1. + std::pow(b, 2.))
                - 4. * b * std::atan(1./b); // (3.46)
        Phi2 = 2. * (2./3. + std::log(std::pow(a, 2.) * std::pow(Z, 2./3.) * std::pow(ME, 2.))) - 2. * std::log(1. + std::pow(b, 2.))
                + 8. * std::pow(b, 2.) * (1 - b * std::atan(1./b) - 0.75 * std::log(1 + std::pow(b, -2.)) ); // (3.47)
        Psi1 = 2. * (1. + std::log(std::pow(ap, 2.) * std::pow(Z, 4./3.) * std::pow(ME, 2.))) - 2. * std::log(1 + std::pow(bp, 2.))
                - 4. * bp * std::atan(1./b); // (3.48)
        Psi2 = 2. * (2./3. + std::log(std::pow(ap, 2.) * std::pow(Z, 4./3.) * std::pow(ME, 2.))) - 2. * std::log(1. + std::pow(bp, 2.))
                + 8. * std::pow(bp, 2.) * (1. - b * std::atan(1./bp) - 0.75 * std::log(1 + std::pow(bp, -2.))); // (3.49)
    }else{
        //Heavier elements
        double gamma = 200. * delta / (ME * std::pow(Z, 1./3.)); // (3.30)
        double epsilon = 200. * delta / (ME * std::pow(Z, 2./3.)); // (3.31)

        Phi1 = 20.863 - 2. * std::log( 1. + std::pow(0.55846 * gamma, 2.) )
                - 4. * ( 1. - 0.6 * std::exp(-0.9 * gamma) - 0.4 * std::exp(-1.5 * gamma) ); // (3.38)
        Phi2 = Phi1 - 2./3. * 1. / (1. + 6.5 * gamma + 6. * std::pow(gamma, 2.)); // (3.39)
        Psi1 = 28.340 - 2. * std::log(1. + std::pow(3.621 * epsilon, 2.))
                - 4. * (1. - 0.7 * std::exp(-8. * epsilon) - 0.3 * std::exp(-29.2 * epsilon)); // (3.40)
        Psi2 = Psi1 - 2./3. * 1. / (1. + 40. * epsilon + 400 * std::pow(epsilon, 2.)); // (3.41)
    }

    double z = std::pow(Z / 137., 2.);
    double f = 1.202 * z - 1.0369 * std::pow(z, 2.) + 1.008 * std::pow(z, 3.)/(1. + z); // (3.3)

    double aux = ( 4./3. * std::pow(x, 2.) - 4./3. * x + 1. )
            * ( std::pow(Z, 2.) * (Phi1 - 4./3. * logZ - 4. * f) + Z * (Psi1 - 8./3. * logZ) )
            - 2./3. * x * (1. - x) * ( std::pow(Z, 2.) * (Phi1 - Phi2) + Z * (Psi1 - Psi2) ); // (3.9)

    aux *= ALPHA * std::pow(RE, 2.) / k; // (3.9)

    double p = std::sqrt( std::pow(x * k, 2.) - std::pow(ME, 2.) ); // electron momentum
    aux *= x * std::pow(k, 2.) / p; // conversion from differential cross section in electron momentum to x

    return medium_->GetMolDensity() * components_[component_index_]->GetAtomInMolecule() * aux; //TODO what are the real factors here, those are just guesses
}


const std::string PhotoPairTsai::name_ = "PhotoPairTsai";

