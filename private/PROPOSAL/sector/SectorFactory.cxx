
#include "PROPOSAL/sector/SectorFactory.h"
#include "PROPOSAL/sector/Sector.h"
// #include "PROPOSAL/geometry/Sphere.h"
// #include "PROPOSAL/geometry/Box.h"
// #include "PROPOSAL/geometry/Cylinder.h"
// #include "PROPOSAL/Output.h"
#include "PROPOSAL/medium/Medium.h"

using namespace PROPOSAL;

SectorFactory::Definition::Definition()
    : e_cut(500)
    , v_cut(0.05)
    , medium_def()
    , geometry_def()
{
}

SectorFactory::Definition::~Definition()
{
}

Sector* SectorFactory::CreateSector(Particle& particle, const Definition& def)
{
    Medium* med = MediumFactory::Get().CreateMedium(def.medium_def);
    Geometry* geometry = GeometryFactory::Get().CreateGeometry(def.geometry_def);
    EnergyCutSettings cuts(def.e_cut, def.v_cut);

    Sector* sec = new Sector(particle, *med, cuts, *geometry, def);

    delete med;
    delete geometry;

    return sec;
}

Sector* SectorFactory::CreateSector(Particle& particle, const Definition& def, const InterpolationDef& interpolation_def)
{
    Medium* med = MediumFactory::Get().CreateMedium(def.medium_def);
    Geometry* geometry = GeometryFactory::Get().CreateGeometry(def.geometry_def);
    EnergyCutSettings cuts(def.e_cut, def.v_cut);

    Sector* sec = new Sector(particle, *med, cuts, *geometry, def, interpolation_def);

    delete med;
    delete geometry;

    return sec;
}
