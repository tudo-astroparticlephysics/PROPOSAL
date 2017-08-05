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

namespace PROPOSAL {

namespace Components {

class Component
{
    public:
    Component(std::string name, double charge, double atomicNum, double atomInMolecule);
    Component(const Component&);
    virtual ~Component(){};

    void swap(Component&);
    virtual Component* clone() const = 0;

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

class Oxygen : public Component
{
    public:
    Oxygen(double atomInMolecule = 1.0);
    virtual ~Oxygen() {}

    Component* clone() const { return new Oxygen(*this); };
};

class Hydrogen : public Component
{
    public:
    Hydrogen(double atomInMolecule = 2.0);
    virtual ~Hydrogen() {}

    Component* clone() const { return new Hydrogen(*this); };
};

class Natrium : public Component
{
    public:
    Natrium(double atomInMolecule = 1.0);
    virtual ~Natrium() {}

    Component* clone() const { return new Natrium(*this); };
};

class Chloride : public Component
{
    public:
    Chloride(double atomInMolecule = 1.0);
    virtual ~Chloride() {}

    Component* clone() const { return new Chloride(*this); };
};

class StandardRock : public Component
{
    public:
    StandardRock(double atomInMolecule = 1.0);
    virtual ~StandardRock() {}

    Component* clone() const { return new StandardRock(*this); };
};

class FrejusRock : public Component
{
    public:
    FrejusRock(double atomInMolecule = 1.0);
    virtual ~FrejusRock() {}

    Component* clone() const { return new FrejusRock(*this); };
};

class Iron : public Component
{
    public:
    Iron(double atomInMolecule = 1.0);
    virtual ~Iron() {}

    Component* clone() const { return new Iron(*this); };
};

class Lead : public Component
{
    public:
    Lead(double atomInMolecule = 1.0);
    virtual ~Lead() {}

    Component* clone() const { return new Lead(*this); };
};

class Copper : public Component
{
    public:
    Copper(double atomInMolecule = 1.0);
    virtual ~Copper() {}

    Component* clone() const { return new Copper(*this); };
};

class Uranium : public Component
{
    public:
    Uranium(double atomInMolecule = 1.0);
    virtual ~Uranium() {}

    Component* clone() const { return new Uranium(*this); };
};

class Nitrogen : public Component
{
    public:
    Nitrogen(double atomInMolecule = 1.0);
    virtual ~Nitrogen() {}

    Component* clone() const { return new Nitrogen(*this); };
};

class Arsenic : public Component
{
    public:
    Arsenic(double atomInMolecule = 1.0);
    virtual ~Arsenic() {}

    Component* clone() const { return new Arsenic(*this); };
};

class Carbon : public Component
{
    public:
    Carbon(double atomInMolecule = 1.0);
    virtual ~Carbon() {}

    Component* clone() const { return new Carbon(*this); };
};

class Potassium : public Component
{
    public:
    Potassium(double atomInMolecule = 1.0);
    virtual ~Potassium() {}

    Component* clone() const { return new Potassium(*this); };
};

class Magnesium : public Component
{
    public:
    Magnesium(double atomInMolecule = 1.0);
    virtual ~Magnesium() {}

    Component* clone() const { return new Magnesium(*this); };
};

class Calcium : public Component
{
    public:
    Calcium(double atomInMolecule = 1.0);
    virtual ~Calcium() {}

    Component* clone() const { return new Calcium(*this); };
};

class Sulfur : public Component
{
    public:
    Sulfur(double atomInMolecule = 1.0);
    virtual ~Sulfur() {}

    Component* clone() const { return new Sulfur(*this); };
};

} // namesapce Components

} // namesapce PROPOSAL
