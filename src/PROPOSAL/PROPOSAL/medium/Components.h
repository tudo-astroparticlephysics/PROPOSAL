
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

#include <map>
#include <memory>
#include <string>
#include <vector>

#define COMPONENT_DEC(cls, ATOMS)                                              \
    class cls : public PROPOSAL::Component {                                   \
    public:                                                                    \
        cls(double atomInMolecule = ATOMS);                                    \
    };

namespace PROPOSAL {

class Component {
public:
    Component() = default;
    Component(std::string name, double charge, double atomicNum,
        double atomInMolecule, double ionizationEnergy);
    virtual ~Component() = default;

    friend bool operator==(Component const&, Component const&) noexcept;
    bool operator!=(const Component&) const;
    friend std::ostream& operator<<(std::ostream&, Component const&) noexcept;

    std::string GetName() const { return name_; }
    double GetNucCharge() const { return nucCharge_; }
    double GetAtomicNum() const { return atomicNum_; }
    double GetAtomInMolecule() const { return atomInMolecule_; }
    double GetLogConstant() const { return logConstant_; }
    double GetBPrime() const { return bPrime_; }
    double GetAverageNucleonWeight() const { return averageNucleonWeight_; }
    double GetWoodSaxon() const { return wood_saxon_; }
    double GetIonizationEnergy() const { return ionizationEnergy_; }

    size_t GetHash() const noexcept;
    static Component GetComponentForHash(size_t);

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

    /*!
     * Woods-Saxon potential calculation
     *
     * \param   r
     * \return  value of the Woods-Saxon potential
     *
     * Calculation of the integral defined in
     * Butkevich Mikheyev JETP 95 (2002), 11 eq. 46
     * This can be integrated analytically.
     */
    double WoodSaxonPotential(double r0);

    // Passed to constructor
    std::string name_;
    double nucCharge_;        ///< nucleus charge
    double atomicNum_;        ///< molar mass [g/mol]
    double atomInMolecule_;   ///< number of atoms in a molecule
    double ionizationEnergy_; ///< ionization energy in eV

    // Calculated in constructor
    double logConstant_ = 0;          ///< radiation logarithm constant B
    double bPrime_ = 0;               ///< radiation logarithm constant bPrime
    double averageNucleonWeight_ = 0; ///< average nucleon weight in a nucleus
                                      ///< [MeV]
    double wood_saxon_ = 0;           ///< Woods-Saxon potential factor
    size_t hash;

private:
    static std::unique_ptr<std::map<size_t, Component>> component_map;
};

bool operator==(Component const&, Component const&) noexcept;
inline bool operator==(std::shared_ptr<Component> const& lhs,
    std::shared_ptr<Component> const& rhs) noexcept
{
    if (lhs != nullptr && rhs != nullptr)
        return *lhs == *rhs;
    if (lhs == nullptr && rhs == nullptr)
        return true;
    return false;
}

namespace Components {

    COMPONENT_DEC(Hydrogen, 2.0)  // 1
    COMPONENT_DEC(Carbon, 1.0)    // 6
    COMPONENT_DEC(Nitrogen, 1.0)  // 7
    COMPONENT_DEC(Oxygen, 1.0)    // 8
    COMPONENT_DEC(Sodium, 1.0)    // 11
    COMPONENT_DEC(Magnesium, 1.0) // 12
    COMPONENT_DEC(Aluminium, 1.0) // 13
    COMPONENT_DEC(Silicon, 1.0)   // 14
    COMPONENT_DEC(Sulfur, 1.0)    // 16
    COMPONENT_DEC(Chlorine, 1.0)  // 17
    COMPONENT_DEC(Argon, 1.0)     // 18
    COMPONENT_DEC(Potassium, 1.0) // 19
    COMPONENT_DEC(Calcium, 1.0)   // 20
    COMPONENT_DEC(Iron, 1.0)      // 26
    COMPONENT_DEC(Cobalt, 1.0)    // 27
    COMPONENT_DEC(Nickel, 1.0)    // 28
    COMPONENT_DEC(Copper, 1.0)    // 29
    COMPONENT_DEC(Arsenic, 1.0)   // 33
    COMPONENT_DEC(Lead, 1.0)      // 82
    COMPONENT_DEC(Uranium, 1.0)   // 92
    COMPONENT_DEC(StandardRock, 1.0)
    COMPONENT_DEC(FrejusRock, 1.0)

} // namespace Components

using component_list = std::vector<Component>;

double calculate_proton_massnumber_fraction(
    const component_list& comp_list) noexcept;

} // namespace PROPOSAL

namespace std {
template <> struct hash<PROPOSAL::Component> {
    std::size_t operator()(PROPOSAL::Component const& comp) const noexcept
    {
        return comp.GetHash();
    }
};

template <> struct hash<std::shared_ptr<PROPOSAL::Component>> {
    std::size_t operator()(
        std::shared_ptr<PROPOSAL::Component> const& comp) const noexcept
    {
        if (comp)
            return comp->GetHash();
        return 0;
    }
};
}

#undef COMPONENT_DEC
