
#pragma once

#include "PROPOSAL/crossection/parametrization/Photonuclear.h"

namespace PROPOSAL {

/******************************************************************************
*                         PhotoRealPhotonAssumption                           *
******************************************************************************/

class PhotoRealPhotonAssumption : public Photonuclear
{
    public:
    PhotoRealPhotonAssumption(const ParticleDef&,
                              const Medium&,
                              const EnergyCutSettings&,
                              const RealPhoton&,
                              Definition = Definition());
    PhotoRealPhotonAssumption(const PhotoRealPhotonAssumption&);
    virtual ~PhotoRealPhotonAssumption();

    virtual Parametrization* clone() const = 0;

    // ----------------------------------------------------------------- //
    // Public methods
    // ----------------------------------------------------------------- //

    virtual double DifferentialCrossSection(double energy, double v);

    virtual double CalculateParametrization(double nu) = 0;
    double NucleusCrossSectionCaldwell(double nu);

    protected:
    RealPhoton* hardBB_;
};

/******************************************************************************
*                            Zeus Parametrization                            *
******************************************************************************/

class PhotoZeus: public PhotoRealPhotonAssumption
{
    public:
    PhotoZeus(const ParticleDef&,
              const Medium&,
              const EnergyCutSettings&,
              const RealPhoton& hardBB,
              Definition = Definition());
    PhotoZeus(const PhotoZeus&);
    virtual ~PhotoZeus();

    virtual Parametrization* clone() const { return new PhotoZeus(*this); }

    virtual double CalculateParametrization(double nu);
};

/******************************************************************************
*                      Bezrukov Bugaev Parametrization                       *
******************************************************************************/

class PhotoBezrukovBugaev: public PhotoRealPhotonAssumption
{
    public:
    PhotoBezrukovBugaev(const ParticleDef&,
                        const Medium&,
                        const EnergyCutSettings&,
                        const RealPhoton& hardBB,
                        Definition = Definition());
    PhotoBezrukovBugaev(const PhotoBezrukovBugaev&);
    virtual ~PhotoBezrukovBugaev();

    virtual Parametrization* clone() const { return new PhotoBezrukovBugaev(*this); }

    virtual double CalculateParametrization(double nu);
};

/******************************************************************************
*                           Rhode Parametrization                            *
******************************************************************************/

class PhotoRhode : public PhotoRealPhotonAssumption
{
    public:
    PhotoRhode(const ParticleDef&,
               const Medium&,
               const EnergyCutSettings&,
               const RealPhoton& hardBB,
               Definition = Definition());
    PhotoRhode(const PhotoRhode&);
    virtual ~PhotoRhode();

    virtual Parametrization* clone() const { return new PhotoRhode(*this); }

    virtual double CalculateParametrization(double nu);

    private:
    double MeasuredSgN(double e);

    Interpolant* interpolant_;
};

/******************************************************************************
*                          Kokoulin Parametrization                           *
******************************************************************************/

class PhotoKokoulin : public PhotoBezrukovBugaev
{
    public:
    PhotoKokoulin(const ParticleDef&,
                  const Medium&,
                  const EnergyCutSettings&,
                  const RealPhoton& hardBB,
                  Definition = Definition());
    PhotoKokoulin(const PhotoKokoulin&);
    virtual ~PhotoKokoulin();

    virtual Parametrization* clone() const { return new PhotoKokoulin(*this); }

    virtual double CalculateParametrization(double nu);
};

} /* PROPOSAL */
