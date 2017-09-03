
#include "PROPOSAL/sector/SectorFactory.h"
#include "PROPOSAL/sector/Sector.h"
#include "PROPOSAL/sector/SectorIntegral.h"
#include "PROPOSAL/sector/SectorInterpolant.h"
#include "PROPOSAL/geometry/Sphere.h"
#include "PROPOSAL/geometry/Box.h"
#include "PROPOSAL/geometry/Cylinder.h"
#include "PROPOSAL/Output.h"

using namespace PROPOSAL;

SectorFactory::Definition::Definition()
    : e_cut(500)
    , v_cut(0.05)
    , do_interpolation(true)
    , medium(MediumFactory::Water)
    , density_correction(1.0)
    , geometry(GeometryFactory::Sphere)
    , position()
    , inner_radius(0.0)
    , radius(0.0)
    , width(0.0)
    , height(0.0)
    , depth(0.0)
{
}

SectorFactory::Definition::~Definition()
{
}

Sector* SectorFactory::CreateSector(PROPOSALParticle& particle, const Definition& def)
{
    Medium* med = MediumFactory::Get().CreateMedium(def.medium, def.density_correction);
    Geometry* geometry = GeometryFactory::Get().CreateGeometry(def.geometry);
    EnergyCutSettings cuts(def.e_cut, def.v_cut);

    if (Sphere* sphere = dynamic_cast<PROPOSAL::Sphere*>(geometry))
    {
        sphere->SetPosition(def.position);
        sphere->SetRadius(def.radius);
        sphere->SetInnerRadius(def.inner_radius);
    }
    else if (Box* box = dynamic_cast<PROPOSAL::Box*>(geometry))
    {
        box->SetPosition(def.position);
        box->SetX(def.width);
        box->SetY(def.height);
        box->SetZ(def.depth);
    }

    else if (Cylinder* cylinder = dynamic_cast<PROPOSAL::Cylinder*>(geometry))
    {
        cylinder->SetPosition(def.position);
        cylinder->SetRadius(def.radius);
        cylinder->SetInnerRadius(def.inner_radius);
    }
    else
    {
        log_fatal("Geometry %s not registerd!", typeid(geometry).name());
    }

    Sector* sec;

    if (def.do_interpolation)
    {
        sec = new SectorInterpolant(particle, *med, *geometry, cuts, def);
    }
    else
    {
        sec = new SectorIntegral(particle, *med, *geometry, cuts, def);
    }

    delete med;
    delete geometry;

    return sec;
}
