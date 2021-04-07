
/******************************************************************************
 *                                                                            *
 * This file is part of the simulation tool PROPOSAL.                         *
 *                                                                            *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,          *
 *                    Chair Experimental Physics 5b                           *
 *                                                                            *
 * This software may be modified and distributed under the terms of a         *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE".                                     *
 *                                                                            *
 * Modifcations to the LGPL License:                                          *
 *                                                                            *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the       *
 *         following reference:                                               *
 *                                                                            *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001                                          *
 *                                                                            *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *         GitHub webpage                                                     *
 *                                                                            *
 *         "https://github.com/tudo-astroparticlephysics/PROPOSAL"            *
 *                                                                            *
 ******************************************************************************/


#pragma once

#include <functional>
#include <random>
#include <iostream>


#ifdef ICECUBE_PROJECT
#include <phys-services/I3RandomService.h>
#endif

namespace PROPOSAL {

// ----------------------------------------------------------------------------
/// @brief Random number generator
// ----------------------------------------------------------------------------
class RandomGenerator
{
public:
    static RandomGenerator& Get()
    {
        static RandomGenerator instance;
        return instance;
    }

    // ----------------------------------------------------------------------------
    /// @brief Execute the given rng to get a random number
    ///
    /// @return random number
    // ----------------------------------------------------------------------------
    double RandomDouble();

    void SetSeed(int seed);

    // ----------------------------------------------------------------------------
    /// @brief Serialize the rng to a stream
    ///
    /// Useful for debuging to save a specific state.
    /// Only supported for internal used rng from the standard libraries.
    ///
    /// @param std::ostream
    // ----------------------------------------------------------------------------
    void Serialize(std::ostream&);

    // ----------------------------------------------------------------------------
    /// @brief Deserialize the rng from a stream
    ///
    /// Useful for debuging to get back a specific state.
    /// Only supported for internal used rng from the standard libraries.
    ///
    /// @param std::ostream
    // ----------------------------------------------------------------------------
    void Deserialize(std::istream&);

    /** @brief Set a custom random number generator
     *
     * Classes that contain other subclasses of MathModel should
     * override this to pass the new RNG on to their members.
     */
    virtual void SetRandomNumberGenerator(std::function<double()>& f);

#ifdef ICECUBE_PROJECT
    virtual void SetI3RandomNumberGenerator(I3RandomServicePtr random);
#endif

    /** @brief Set a the default random number generator
     *
     * The default random number generator is the
     * std mersenne twister std::mt19937.
     */
    virtual void SetDefaultRandomNumberGenerator();

private:
    RandomGenerator();
    virtual ~RandomGenerator();

    static double DefaultRandomDouble();

    static std::mt19937 rng_;
    static std::uniform_real_distribution<double> uniform_distribution;
    std::function<double()> random_function;
#ifdef ICECUBE_PROJECT
    I3RandomService* i3random_gen_;
#endif
};

} // namespace PROPOSAL
