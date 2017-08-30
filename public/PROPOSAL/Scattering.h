/*! \file   Scattering.h
*   \brief  Header file for the Scattering bug routines.
*
*   This version has a major bug and produces too small scattering angles.
*
*   \date   2013.08.19
*   \author Tomasz Fuchs
*/
#pragma once

#include <boost/function.hpp>

#include <vector>
#include <map>
#include <string>

namespace PROPOSAL {

class PROPOSALParticle;
class CrossSections;

class Scattering
{
    public:
    Scattering();
    virtual ~Scattering();

    virtual Scattering* clone() const = 0; // virtual constructor idiom (used for deep copies)

    void Scatter(PROPOSALParticle&, const std::vector<CrossSections*>&, double dr, double ei, double ef);

    virtual void EnableInterpolation(const PROPOSALParticle&, const std::vector<CrossSections*>&, std::string path = "") = 0;
    virtual void DisableInterpolation() = 0;

    protected:
    struct RandomAngles
    {
        double sx, sy, tx, ty;
    };

    virtual RandomAngles CalculateRandomAngle(const PROPOSALParticle&, const std::vector<CrossSections*>&, double dr, double ei, double ef) = 0;
};

class ScatteringFactory
{
    public:

    struct ScatteringModel
    {
        enum Enum
        {
            Default = 0,
            Moliere,
            MoliereFirstOrder
        };
    };

    typedef boost::function<Scattering* (void)> RegisterFunction;
    typedef std::map<std::string, boost::function<Scattering* (void)> > ScatteringMapString;
    typedef std::map<ScatteringModel::Enum, boost::function<Scattering* (void)> > ScatteringMapEnum;
    typedef std::map<std::string, ScatteringModel::Enum> MapStringToEnum;

    Scattering* CreateScattering(const std::string&);
    Scattering* CreateScattering(const ScatteringModel::Enum);

    ScatteringModel::Enum GetEnumFromString(const std::string&);

    static ScatteringFactory& Get()
    {
        static ScatteringFactory instance;
        return instance;
    }

    private:
    ScatteringFactory();
    ~ScatteringFactory();

    void Register(const std::string& name, RegisterFunction);
    void Register(ScatteringModel::Enum, RegisterFunction);
    void Register(const std::string&, ScatteringModel::Enum);

    // std::map<std::string, Scattering* (*)(void)> scattering_map;
    ScatteringMapString scattering_map_str_;
    ScatteringMapEnum scattering_map_enum_;
    MapStringToEnum map_string_to_enum;
};

}
