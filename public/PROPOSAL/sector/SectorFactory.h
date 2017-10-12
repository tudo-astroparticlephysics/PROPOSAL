
#pragma once


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

        MediumFactory::Definition medium_def;
        GeometryFactory::Definition geometry_def;

        Definition();
        ~Definition();
    };

    Sector* CreateSector(Particle&, const Definition&);
    Sector* CreateSector(Particle&, const Definition&, const InterpolationDef&);

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

