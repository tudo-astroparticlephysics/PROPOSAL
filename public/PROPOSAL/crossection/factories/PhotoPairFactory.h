
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
#include <map>
#include <string>
#include <memory>
#include <sstream>

#include "PROPOSAL/methods.h"

namespace PROPOSAL {

    class CrossSection;
    class PhotoPairProduction;
    struct ParticleDef;
    class EnergyCutSettings;
    class Medium;
    class PhotoAngleDistribution;

    class PhotoPairFactory
    {
    public:
        // --------------------------------------------------------------------- //
        // Declare usable enums
        // --------------------------------------------------------------------- //

        enum Enum
        {
            Fail = 0,
            None,
            Tsai
        };

        enum PhotoAngle
        {
            PhotoAngleTsaiIntegral,
            PhotoAngleNoDeflection,
            PhotoAngleEGS
        };

        struct Definition
        {
            Definition()
                    : parametrization(None)
                    , photoangle(PhotoAngleNoDeflection)
                    , multiplier(1.0)
            {
            }

            bool operator==(const PhotoPairFactory::Definition& def) const
            {
                if (parametrization != def.parametrization)
                    return false;
                else if (multiplier != def.multiplier)
                    return false;
                else if (photoangle != def.photoangle)
                    return false;

                return true;
            }

            bool operator!=(const PhotoPairFactory::Definition& def) const
            {
                return !(*this == def);
            }

            friend std::ostream& operator<<(std::ostream& os, PROPOSAL::PhotoPairFactory::Definition const& definition)
            {
                std::stringstream ss;
                ss << " Photo Pair Production Definition (" << &definition << ") ";
                os << Helper::Centered(60, ss.str()) << '\n';

                os << "Parametrization: " << definition.parametrization << std::endl;
                os << "Multiplier: " << definition.multiplier << std::endl;
                os << "PhotoAngle: " << definition.photoangle << std::endl;

                os << Helper::Centered(60, "");
                return os;
            }

            Enum parametrization;
            PhotoAngle photoangle;
            double multiplier;
        };

        // --------------------------------------------------------------------- //
        // Typedefs for readablitiy
        // --------------------------------------------------------------------- //

        typedef std::function<
                PhotoPairProduction*(const ParticleDef&, std::shared_ptr<const Medium>, double multiplier)>
                RegisterFunction;

        typedef std::map<std::string, RegisterFunction> PhotoPairMapString;
        typedef std::map<Enum, RegisterFunction> PhotoPairMapEnum;
        typedef Helper::Bimap<std::string, Enum> BimapStringEnum;

        typedef std::function<PhotoAngleDistribution*(const ParticleDef&, std::shared_ptr<const Medium>)> RegisterPhotoAngleFunction;

        typedef std::map<std::string, RegisterPhotoAngleFunction> PhotoAngleMapString;
        typedef std::map<PhotoAngle, RegisterPhotoAngleFunction> PhotoAngleMapEnum;
        typedef Helper::Bimap<std::string, PhotoAngle> BimapStringPhotoAngleEnum;

        // --------------------------------------------------------------------- //
        // Most general creation
        // --------------------------------------------------------------------- //

        CrossSection* CreatePhotoPair(const ParticleDef&,
                                    std::shared_ptr<const Medium>,
                                    const Definition&) const;

        CrossSection* CreatePhotoPair(const ParticleDef&,
                                    std::shared_ptr<const Medium>,
                                    const Definition&,
                                    InterpolationDef) const;

        // ----------------------------------------------------------------------------
        /// @brief string to enum conversation for photopair parametrizations
        // ----------------------------------------------------------------------------
        Enum GetEnumFromString(const std::string&);

        // ----------------------------------------------------------------------------
        /// @brief enum to string conversation for photopair parametrizations
        // ----------------------------------------------------------------------------
        std::string GetStringFromEnum(const Enum&);

        // --------------------------------------------------------------------- //
        // Singleton pattern
        // --------------------------------------------------------------------- //

        static PhotoPairFactory& Get()
        {
            static PhotoPairFactory instance;
            return instance;
        }

        // --------------------------------------------------------------------- //
        // PhotoAngle
        // --------------------------------------------------------------------- //

        PhotoAngleDistribution* CreatePhotoAngleDistribution(const std::string&, const ParticleDef&, std::shared_ptr<const Medium> medium);
        PhotoAngleDistribution* CreatePhotoAngleDistribution(const PhotoAngle&, const ParticleDef&, std::shared_ptr<const Medium> medium);

        PhotoAngle GetPhotoAngleEnumFromString(const std::string&);
        std::string GetStringFromPhotoAngleEnum(const PhotoAngle&);

    private:
        PhotoPairFactory();
        ~PhotoPairFactory();

        // ----------------------------------------------------------------------------
        /// @brief Register PhotoPair parametrizations
        ///
        /// @param name
        /// @param Enum
        /// @param RegisterFunction
        // ----------------------------------------------------------------------------
        void Register(const std::string& name, Enum, RegisterFunction);

        PhotoPairMapString photopair_map_str_;
        PhotoPairMapEnum photopair_map_enum_;
        BimapStringEnum string_enum_;

        // PhotoAngle

        void RegisterPhotoAngle(const std::string& name, const PhotoAngle&, RegisterPhotoAngleFunction);

        PhotoAngleMapString photo_angle_map_str_;
        PhotoAngleMapEnum photo_angle_map_enum_;
        BimapStringPhotoAngleEnum string_photo_angle_enum_;
    };

} // namespace PROPOSAL
