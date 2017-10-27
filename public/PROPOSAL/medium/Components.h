/*! \file   Components.h
*   \brief  Header file for definition of the components class object.
*
*   For more details see the class documentation.
*
*   \date   Sat Aug  5 14:47:16 CEST 2017
*   \author Mario Dunsch
*/

#pragma once

#include <string>

#define COMPONENT_DEC(cls, ATOMS)                                                                                      \
    class cls : public Component                                                                                       \
    {                                                                                                                  \
        public:                                                                                                        \
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
    Component(const Component&);
    virtual ~Component(){};

    void swap(Component&);
    virtual Component* clone() const { return new Component(*this); }; // Prototyping/Virtual constructor idiom (used for deep copies)

    Component& operator=(const Component&);
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
    double atomicNum_;      ///< atomic number
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

} // namesapce Components

} // namesapce PROPOSAL

#undef COMPONENT_DEC
