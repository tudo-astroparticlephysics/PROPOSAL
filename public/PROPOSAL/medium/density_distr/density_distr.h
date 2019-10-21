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
#include <exception>
#include <functional>
#include <string>
#include "PROPOSAL/math/Vector3D.h"

namespace PROPOSAL {
class Axis {
   public:
    Axis();
    Axis(const Vector3D& fp0, const Vector3D& fAxis);
    Axis(const Axis&);

    virtual ~Axis() {};

    bool operator==(const Axis& axis) const;
    bool operator!=(const Axis& axis) const;

    virtual Axis* clone() const = 0;

    virtual double GetDepth(const Vector3D& xi) const = 0;
    virtual double GetEffectiveDistance(const Vector3D& xi,
                                        const Vector3D& direction) const = 0;

    Vector3D GetAxis() const { return fAxis_; };
    Vector3D GetFp0() const { return fp0_; };

   protected:
    Vector3D fAxis_;
    Vector3D fp0_;
};
}  // namespace PROPOSAL

namespace PROPOSAL {
class RadialAxis : public Axis {
   public:
    RadialAxis();
    RadialAxis(const Vector3D& fAxis, const Vector3D& fp0);

    Axis* clone() const override { return new RadialAxis(*this); };
    ~RadialAxis() {};

    double GetDepth(const Vector3D& xi) const override;
    double GetEffectiveDistance(const Vector3D& xi, const Vector3D& direction) const override;
};
}  // namespace PROPOSAL

namespace PROPOSAL {
class CartesianAxis : public Axis {
   public:
    CartesianAxis();
    CartesianAxis(const Vector3D& fAxis, const Vector3D& fp0);
    ~CartesianAxis() {};

    Axis* clone() const override { return new CartesianAxis(*this); };

    double GetDepth(const Vector3D& xi) const override;
    double GetEffectiveDistance(const Vector3D& xi, const Vector3D& direction) const override;
};
}  // namespace PROPOSAL

namespace PROPOSAL {
class DensityException : public std::exception {
   public:
    DensityException(const char* m) : message_(m){};
    const char* what() const throw() { return message_.c_str(); };

   private:
    std::string message_;
};
}  // namespace PROPOSAL

namespace PROPOSAL {
class Density_distr {
   public:
    Density_distr();
    Density_distr(const Axis& axis);
    Density_distr(const Density_distr&);

    virtual ~Density_distr() { delete axis_; };

    virtual bool operator==(const Density_distr& dens_distr) const;
    virtual bool operator!=(const Density_distr& dens_distr) const;
    virtual bool compare(const Density_distr& dens_distr) const = 0;


    virtual Density_distr* clone() const = 0;

    virtual double Correct(const Vector3D& xi,
                           const Vector3D& direction,
                           double res,
                           double distance_to_border) const = 0;
    virtual double Integrate(const Vector3D& xi,
                             const Vector3D& direction,
                             double l) const = 0;
    virtual double Calculate(const Vector3D& xi,
                             const Vector3D& direction,
                             double distance) const = 0;
    virtual double Evaluate(const Vector3D& xi) const = 0;

    const Axis& GetAxis() const { return *axis_; }

   protected:
    Axis* axis_;
};
}  // namespace PROPOSAL
