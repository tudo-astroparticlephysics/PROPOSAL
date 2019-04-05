
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

#include <string>

#define COMPONENT_DEC(cls, ATOMS)                                                                                      \
    class cls : public Component                                                                                       \
    {                                                                                                                  \
    public:                                                                                                            \
        cls(double atomInMolecule = ATOMS);                                                                            \
        cls(const cls&);                                                                                               \
        virtual ~cls() {}                                                                                              \
                                                                                                                       \
        virtual Component* clone() const { return new cls(*this); };                                                   \
    };

namespace PROPOSAL {

namespace Components {

class Component
{
public:
    Component(std::string name, double charge, double atomicNum, double atomInMolecule);
    Component(const Component&) = default;
    virtual ~Component() = default;
    Component& operator=(const Component&) = default;

    virtual Component* clone() const { return new Component(*this); };

    bool operator==(const Component&) const;
    bool operator!=(const Component&) const;
    friend std::ostream& operator<<(std::ostream&, Component const&);

    // Getter
    std::string GetName() const { return name_; }
    double GetNucCharge() const { return nucCharge_; }
    double GetAtomicNum() const { return atomicNum_; }
    double GetAtomInMolecule() const { return atomInMolecule_; }
    double GetLogConstant() const { return logConstant_; }
    double GetBPrime() const { return bPrime_; }
    double GetAverageNucleonWeight() const { return M_; }
    double GetMN() const { return mN_; }
    double GetR0() const { return r0_; }

protected:
    /*!
     * set the value of radiation logarithm constant B
     *
     * \param   i   nucleon number
     * \return  value of radiation logarithm constant B
     */
    void SetLogConstant();

    /*!
     * set the value of radiation logarithm constant bPrime
     *
     * \param   i   nucleon number
     * \return  value of radiation logarithm constant bPrime
     */
    void SetBPrime();

    /*!
     * Woods-Saxon potential calculation - function to integrate
     *
     * \param   r
     * \return  value of the Woods-Saxon potential
     */
    double FunctionToIntegral(double r);

    // Passed to constructor
    std::string name_;
    double nucCharge_;      ///< nucleus charge
    double atomicNum_;      ///< molar mass [g/mol]
    double atomInMolecule_; ///< number of atoms in a molecule

    // Calculated in constructor
    double logConstant_; ///< radiation logarithm constant B
    double bPrime_;      ///< radiation logarithm constant bPrime
    double M_;           ///< average nucleon weight in a nucleus [MeV]
    double mN_;          ///< Woods-Saxon potential factor
    double r0_;          // //TODO(mario): Must really be stored? Thu 2017/08/03
};

COMPONENT_DEC(Hydrogen, 2.0)
COMPONENT_DEC(Carbon, 1.0)
COMPONENT_DEC(Nitrogen, 1.0)
COMPONENT_DEC(Oxygen, 1.0)
COMPONENT_DEC(Sodium, 1.0)
COMPONENT_DEC(Magnesium, 1.0)
COMPONENT_DEC(Sulfur, 1.0)
COMPONENT_DEC(Chlorine, 1.0)
COMPONENT_DEC(Argon, 1.0)
COMPONENT_DEC(Potassium, 1.0)
COMPONENT_DEC(Calcium, 1.0)
COMPONENT_DEC(Iron, 1.0)
COMPONENT_DEC(Copper, 1.0)
COMPONENT_DEC(Lead, 1.0)
COMPONENT_DEC(Uranium, 1.0)
COMPONENT_DEC(StandardRock, 1.0)
COMPONENT_DEC(FrejusRock, 1.0)

} // namespace Components

} // namespace PROPOSAL

#undef COMPONENT_DEC
