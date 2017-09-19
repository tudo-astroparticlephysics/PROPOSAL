
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
                              double multiplier);
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
              double multiplier);
    PhotoZeus(const PhotoZeus&);
    virtual ~PhotoZeus();

    Parametrization* clone() const { return new PhotoZeus(*this); }
    static Parametrization* create(const ParticleDef&,
                            const Medium&,
                            const EnergyCutSettings&,
                            const RealPhoton&,
                            double multiplier);

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
                        double multiplier);
    PhotoBezrukovBugaev(const PhotoBezrukovBugaev&);
    virtual ~PhotoBezrukovBugaev();

    Parametrization* clone() const { return new PhotoBezrukovBugaev(*this); }
    static Parametrization* create(const ParticleDef&,
                            const Medium&,
                            const EnergyCutSettings&,
                            const RealPhoton&,
                            double multiplier);

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
               double multiplier);
    PhotoRhode(const PhotoRhode&);
    virtual ~PhotoRhode();

    Parametrization* clone() const { return new PhotoRhode(*this); }
    static Parametrization* create(const ParticleDef&,
                            const Medium&,
                            const EnergyCutSettings&,
                            const RealPhoton&,
                            double multiplier);

    double CalculateParametrization(double nu);

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
                  double multiplier);
    PhotoKokoulin(const PhotoKokoulin&);
    virtual ~PhotoKokoulin();

    Parametrization* clone() const { return new PhotoKokoulin(*this); }
    static Parametrization* create(const ParticleDef&,
                            const Medium&,
                            const EnergyCutSettings&,
                            const RealPhoton&,
                            double multiplier);

    double CalculateParametrization(double nu);
};

} /* PROPOSAL */
