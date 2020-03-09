
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

#include <functional>
#include <cmath>
#include <fstream>

#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"


namespace PROPOSAL {

    class Interpolant;

    class PhotoPairProduction : public Parametrization
    {
    public:
        PhotoPairProduction(const ParticleDef&, std::shared_ptr<const Medium>, double multiplier);
        PhotoPairProduction(const PhotoPairProduction&);
        virtual ~PhotoPairProduction();

        virtual Parametrization* clone() const = 0;

        // ----------------------------------------------------------------- //
        // Public methods
        // ----------------------------------------------------------------- //

        virtual double DifferentialCrossSection(double energy, double x) = 0;

        virtual IntegralLimits GetIntegralLimits(double energy);

        virtual size_t GetHash() const;

    protected:
        bool compare(const Parametrization&) const;

    };


    class PhotoPairTsai : public PhotoPairProduction
    {
    public:
        typedef std::vector<Interpolant*> InterpolantVec;

        PhotoPairTsai(const ParticleDef&, std::shared_ptr<const Medium>, double multiplier);
        PhotoPairTsai(const PhotoPairTsai&);
        virtual ~PhotoPairTsai();

        virtual Parametrization* clone() const { return new PhotoPairTsai(*this); }
        static PhotoPairProduction* create(const ParticleDef& particle_def,
                                       std::shared_ptr<const Medium> medium,
                                       double multiplier)
        {
            return new PhotoPairTsai(particle_def, medium, multiplier);
        }

        virtual double DifferentialCrossSection(double energy, double x);

        const std::string& GetName() const { return name_; }

    protected:
        static const std::string name_;
        virtual bool compare(const Parametrization&) const;

    };

    class PhotoAngleDistribution
    {
    public:
        PhotoAngleDistribution(const ParticleDef&, std::shared_ptr<const Medium>);
        PhotoAngleDistribution(const PhotoAngleDistribution&);
        virtual ~PhotoAngleDistribution();

        virtual PhotoAngleDistribution* clone() const = 0;

        bool operator==(const PhotoAngleDistribution&) const;
        bool operator!=(const PhotoAngleDistribution&) const;

        struct DeflectionAngles
        {
            double cosphi0, theta0, cosphi1, theta1;
        };

        virtual DeflectionAngles SampleAngles(double energy, double rho, int component_index) = 0;

        // Getter
        const ParticleDef& GetParticleDef() const { return particle_def_; }
        std::shared_ptr<const Medium> GetMedium() const { return medium_; }

        virtual const std::string& GetName() const = 0;
        virtual size_t GetHash() const = 0;

        // Setter
        void SetCurrentComponent(int index) { component_index_ = index; }

    protected:
        virtual bool compare(const PhotoAngleDistribution&) const;

        const ParticleDef particle_def_;
        std::shared_ptr<const Medium> medium_;
        const std::vector<Components::Component>& components_;
        int component_index_;
    };

    class PhotoAngleTsaiIntegral : public PhotoAngleDistribution
    {
    public:
        PhotoAngleTsaiIntegral(const ParticleDef& particle_def, std::shared_ptr<const Medium> medium)
        : PhotoAngleDistribution(particle_def, medium)
        , integral_(IROMB, IMAXS, IPREC)

        {
        }
        PhotoAngleTsaiIntegral(const PhotoAngleTsaiIntegral& photoangle)
        : PhotoAngleDistribution(photoangle)
        , integral_(photoangle.integral_)
        {
        }
        virtual ~PhotoAngleTsaiIntegral() {}

        PhotoAngleDistribution* clone() const {return new PhotoAngleTsaiIntegral(*this); }
        static PhotoAngleDistribution* create(const ParticleDef& particle_def, std::shared_ptr<const Medium> medium) {
            return new PhotoAngleTsaiIntegral(particle_def, medium);
        }

        virtual DeflectionAngles SampleAngles(double energy, double rho, int component_index) ;
        double FunctionToIntegral(double energy, double x, double theta);

        // Getter

        virtual const std::string& GetName() const { return name_; }
        virtual size_t GetHash() const;

    protected:
        Integral integral_;

    private:
        static const std::string name_;

    };

    class PhotoAngleNoDeflection : public PhotoAngleDistribution
    {
    public:
        PhotoAngleNoDeflection(const ParticleDef& particle_def, std::shared_ptr<const Medium> medium)
        : PhotoAngleDistribution(particle_def, medium)
        {
        }
        PhotoAngleNoDeflection(const PhotoAngleNoDeflection& photoangle)
        : PhotoAngleDistribution(photoangle)
        {
        }
        virtual ~PhotoAngleNoDeflection() {}

        PhotoAngleDistribution* clone() const {return new PhotoAngleNoDeflection(*this); }

        static PhotoAngleDistribution* create(const ParticleDef& particle_def, std::shared_ptr<const Medium> medium) {
            return new PhotoAngleNoDeflection(particle_def, medium);
        }

        virtual DeflectionAngles SampleAngles(double energy, double rho, int component_index);

        // Getter

        virtual const std::string& GetName() const { return name_; }
        virtual size_t GetHash() const;

    private:
        static const std::string name_;

    };

    class PhotoAngleEGS : public PhotoAngleDistribution
    {
    public:
        PhotoAngleEGS(const ParticleDef& particle_def, std::shared_ptr<const Medium> medium)
                : PhotoAngleDistribution(particle_def, medium)
        {
        }
        PhotoAngleEGS(const PhotoAngleNoDeflection& photoangle)
                : PhotoAngleDistribution(photoangle)
        {
        }
        virtual ~PhotoAngleEGS() {}

        PhotoAngleDistribution* clone() const {return new PhotoAngleEGS(*this); }

        static PhotoAngleDistribution* create(const ParticleDef& particle_def, std::shared_ptr<const Medium> medium) {
            return new PhotoAngleEGS(particle_def, medium);
        }

        virtual DeflectionAngles SampleAngles(double energy, double rho, int component_index);

        // Getter

        virtual const std::string& GetName() const { return name_; }
        virtual size_t GetHash() const;

    private:
        static const std::string name_;

    };

} // namespace PROPOSAL
