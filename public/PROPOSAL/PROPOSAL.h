
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

#include "PROPOSAL/version.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"

#include "PROPOSAL/geometry/Box.h"
#include "PROPOSAL/geometry/Cylinder.h"
#include "PROPOSAL/geometry/GeometryFactory.h"
#include "PROPOSAL/geometry/Sphere.h"

#include "PROPOSAL/crossection/factories/AnnihilationFactory.h"
#include "PROPOSAL/crossection/factories/BremsstrahlungFactory.h"
#include "PROPOSAL/crossection/factories/ComptonFactory.h"
#include "PROPOSAL/crossection/factories/EpairProductionFactory.h"
#include "PROPOSAL/crossection/factories/IonizationFactory.h"
#include "PROPOSAL/crossection/factories/MupairProductionFactory.h"
#include "PROPOSAL/crossection/factories/PhotonuclearFactory.h"
#include "PROPOSAL/crossection/factories/PhotoPairFactory.h"
#include "PROPOSAL/crossection/factories/WeakInteractionFactory.h"

#include "PROPOSAL/crossection/parametrization/Annihilation.h"
#include "PROPOSAL/crossection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/crossection/parametrization/Compton.h"
#include "PROPOSAL/crossection/parametrization/EpairProduction.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"
#include "PROPOSAL/crossection/parametrization/MupairProduction.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/crossection/parametrization/ParamTables.h"
#include "PROPOSAL/crossection/parametrization/Photonuclear.h"
#include "PROPOSAL/crossection/parametrization/PhotoPairProduction.h"
#include "PROPOSAL/crossection/parametrization/PhotoQ2Integration.h"
#include "PROPOSAL/crossection/parametrization/PhotoRealPhotonAssumption.h"
#include "PROPOSAL/crossection/parametrization/WeakInteraction.h"

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/CrossSectionIntegral.h"
#include "PROPOSAL/crossection/CrossSectionInterpolant.h"

#include "PROPOSAL/crossection/AnnihilationIntegral.h"
#include "PROPOSAL/crossection/AnnihilationInterpolant.h"
#include "PROPOSAL/crossection/BremsIntegral.h"
#include "PROPOSAL/crossection/BremsInterpolant.h"
#include "PROPOSAL/crossection/ComptonIntegral.h"
#include "PROPOSAL/crossection/ComptonInterpolant.h"
#include "PROPOSAL/crossection/EpairIntegral.h"
#include "PROPOSAL/crossection/EpairInterpolant.h"
#include "PROPOSAL/crossection/IonizIntegral.h"
#include "PROPOSAL/crossection/IonizInterpolant.h"
#include "PROPOSAL/crossection/MupairIntegral.h"
#include "PROPOSAL/crossection/MupairInterpolant.h"
#include "PROPOSAL/crossection/PhotoIntegral.h"
#include "PROPOSAL/crossection/PhotoInterpolant.h"
#include "PROPOSAL/crossection/PhotoPairIntegral.h"
#include "PROPOSAL/crossection/PhotoPairInterpolant.h"
#include "PROPOSAL/crossection/WeakIntegral.h"
#include "PROPOSAL/crossection/WeakInterpolant.h"

#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/decay/DecayTable.h"
#include "PROPOSAL/decay/LeptonicDecayChannel.h"
#include "PROPOSAL/decay/ManyBodyPhaseSpace.h"
#include "PROPOSAL/decay/StableChannel.h"
#include "PROPOSAL/decay/TwoBodyPhaseSpace.h"

#include "PROPOSAL/medium/density_distr/density_distr.h"
#include "PROPOSAL/medium/density_distr/density_exponential.h"
#include "PROPOSAL/medium/density_distr/density_homogeneous.h"
#include "PROPOSAL/medium/density_distr/density_polynomial.h"
#include "PROPOSAL/medium/density_distr/density_splines.h"

#include "PROPOSAL/math/Function.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/math/Spline.h"
#include "PROPOSAL/math/TableWriter.h"
#include "PROPOSAL/math/Vector3D.h"

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/particle/ParticleDef.h"

#include "PROPOSAL/propagation_utility/ContinuousRandomizer.h"
#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

#include "PROPOSAL/scattering/Coefficients.h"
#include "PROPOSAL/scattering/Scattering.h"
#include "PROPOSAL/scattering/ScatteringFactory.h"
#include "PROPOSAL/scattering/ScatteringHighland.h"
#include "PROPOSAL/scattering/ScatteringHighlandIntegral.h"
#include "PROPOSAL/scattering/ScatteringMoliere.h"
#include "PROPOSAL/scattering/ScatteringNoScattering.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/PropagatorService.h"
#include "PROPOSAL/Sector.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Secondaries.h"

#if ROOT_SUPPORT
    #include "PROPOSAL/interfaces/root.h"
#endif
