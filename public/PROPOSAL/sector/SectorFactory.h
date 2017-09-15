
#pragma once

#include <vector>
#include <map>
#include <string>

#include "PROPOSAL/sector/Sector.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/geometry/GeometryFactory.h"

namespace  PROPOSAL
{

class SectorFactory
{
    public:

    struct Definition : Sector::Definition
    {
        double e_cut;
        double v_cut;
        bool do_interpolation;
        ScatteringFactory::Enum scattering_model; //!< if true moliere scattering is enabled
        MediumFactory::Enum medium;
        double density_correction;
        GeometryFactory::Enum geometry;
        Vector3D position;
        double inner_radius;
        double radius;
        double width;
        double height;
        double depth;

        Definition();
        ~Definition();
    };

    Sector* CreateSector(PROPOSALParticle&, const Definition&);

    static SectorFactory& Get()
    {
        static SectorFactory instance;
        return instance;
    }

    private:
    SectorFactory(){};
    ~SectorFactory(){};
};

} /*  PROPOSAL */

