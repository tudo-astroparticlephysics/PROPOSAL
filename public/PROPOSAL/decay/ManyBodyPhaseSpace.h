
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

#include <unordered_map>
#include <functional>

#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/particle/ParticleDef.h"

namespace PROPOSAL {

class ManyBodyPhaseSpace : public DecayChannel
{
public:
    class Builder;

    struct PhaseSpaceParameters
    {
        double normalization;
        double weight_max;
        double weight_min;
    };

    struct PhaseSpaceKinematics
    {
        std::vector<double> virtual_masses;
        std::vector<double> momenta;
        double weight;
    };

    typedef std::unordered_map<ParticleDef, PhaseSpaceParameters> ParameterMap;
    typedef std::function<double(const Particle&, const DecayProducts&)> MatrixElementFunction;
    typedef std::function<void(PhaseSpaceParameters&, const ParticleDef&)> EstimateFunction;

public:
    ManyBodyPhaseSpace(std::vector<const ParticleDef*> daughters, MatrixElementFunction ME = nullptr);
    ManyBodyPhaseSpace(const ManyBodyPhaseSpace& mode);
    virtual ~ManyBodyPhaseSpace();

    // No copy and assignemnt -> done by clone
    DecayChannel* clone() const { return new ManyBodyPhaseSpace(*this); }

    // ----------------------------------------------------------------------------
    /// @brief Many body phase space decay
    ///
    /// Calculate decay products with the help of the Raubold Lynch algorithm.
    ///
    /// @param Particle
    ///
    /// @return Vector of particles, the decay products
    // ----------------------------------------------------------------------------
    DecayProducts Decay(const Particle&);

    // ----------------------------------------------------------------------------
    /// @brief Evalutate the matrix element of this channel
    ///
    /// This is the default matrix element for the many body decay.
    /// Its just constant one, which results to a uniform sampling
    /// of phase space events.
    ///
    /// @return matrix element
    // ----------------------------------------------------------------------------
    static double DefaultEvaluate(const Particle&, const DecayProducts&);

    // ----------------------------------------------------------------------------
    /// @brief Evalutate the matrix element of this channel
    ///
    /// @return matrix element
    // ----------------------------------------------------------------------------
    double Evaluate(const Particle&, const DecayProducts&);

    // ----------------------------------------------------------------------------
    /// @brief Sets the uniform flag
    ///
    /// If uniform is true, the momenta will be sampled uniform in the phase space.
    /// This is done by rejection, since the pure raubold lynch algorithm does not
    /// create a uniform phase space. So enabling uniform sampling comes in with
    /// a cost of performance.
    ///
    /// @param uniform
    // ----------------------------------------------------------------------------
    void SetUniformSampling(bool uniform) {uniform_ = uniform;}

    const std::string& GetName() const { return name_; }

private:
    ManyBodyPhaseSpace& operator=(const ManyBodyPhaseSpace&); // Undefined & not allowed


    // ----------------------------------------------------------------------------
    /// @brief Many body phase space decay
    ///
    /// Calculate decay products with the help of the Raubold Lynch algorithm.
    ///
    /// @param Particle
    ///
    /// @return Vector of particles, the decay products
    // ----------------------------------------------------------------------------
    double GenerateEvent(DecayProducts& products, const PhaseSpaceParameters&, const Particle&);

    // ----------------------------------------------------------------------------
    /// @brief Calculate the normalization of the phase space density
    ///
    /// @param parent_mass
    ///
    /// @return \f$ \rho(\Phi) = \frac{{(M - \mu_n)}^{n-2}}{(n-2)!} \cdot \prod_{i=2}^{n} (2\pi P_i)~. \f$
    // ----------------------------------------------------------------------------
    void CalculateNormalization(PhaseSpaceParameters&, double parent_mass);

    // ----------------------------------------------------------------------------
    /// @brief Calculate the maximum weight for the phase space
    ///
    /// @param normalization
    /// @param parent
    ///
    /// This value is need for the rejection method to create a uniform distribution
    /// of the phase space.
    ///
    /// @return maximum weight
    // ----------------------------------------------------------------------------
    void EstimateMaxWeight(PhaseSpaceParameters&, const ParticleDef& parent);

    // ----------------------------------------------------------------------------
    /// @brief Calculate the maximum weight for the phase space
    ///
    /// @param normalization
    /// @param parent_mass
    ///
    /// This value is need for the rejection method to create a uniform distribution
    /// of the phase space.
    ///
    /// @return maximum weight
    // ----------------------------------------------------------------------------
    void SampleEstimateMaxWeight(PhaseSpaceParameters&, const ParticleDef& parent);

    // ----------------------------------------------------------------------------
    /// @brief Calculate the normalization and maximum weight
    ///
    /// @param parent
    ///
    /// For every particle definition the normalization and maximum weight is unique.
    /// Both values will be created and stored in an hash table.
    ///
    /// @return struct containing the normalization and maximum weight
    // ----------------------------------------------------------------------------
    PhaseSpaceParameters GetPhaseSpaceParams(const ParticleDef& parent);


    // ----------------------------------------------------------------------------
    /// @brief Calculate the kinematics for the use in the raubold lynch algorithm
    ///
    /// @param normalization
    /// @param parent_mass
    ///
    /// @return struct containing the weight of the phase space point,
    ///         intermediate momenta and virtual masses for the algorithm.
    // ----------------------------------------------------------------------------
    PhaseSpaceKinematics CalculateKinematics(double normalization, double parent_mass);

    bool compare(const DecayChannel&) const;
    void print(std::ostream&) const;

    std::vector<const ParticleDef*> daughters_;
    std::vector<double> daughter_masses_;

    int number_of_daughters_;
    double sum_daughter_masses_;
    bool uniform_;
    int broad_phase_statistic_;

    MatrixElementFunction matrix_element_;
    bool use_default_matrix_element_;
    EstimateFunction estimate_;

    static const std::string name_;

    ParameterMap parameter_map_;
};

class ManyBodyPhaseSpace::Builder
{
public:
    Builder();
    Builder(const Builder&);
    ~Builder();

    // --------------------------------------------------------------------- //
    // Setter
    // --------------------------------------------------------------------- //

    Builder& addDaughter(const ParticleDef& daughter);
    Builder& setMatrixElement(MatrixElementFunction);
    ManyBodyPhaseSpace build();

private:
    std::vector<const ParticleDef*> daughters_;
    MatrixElementFunction matrix_element_;
};

} // namespace PROPOSAL
