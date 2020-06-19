
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

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"

#include "PROPOSAL/geometry/Box.h"
#include "PROPOSAL/geometry/Cylinder.h"
#include "PROPOSAL/geometry/GeometryFactory.h"
#include "PROPOSAL/geometry/Sphere.h"

#include "PROPOSAL/crossection/parametrization/Annihilation.h"
#include "PROPOSAL/crossection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/crossection/parametrization/Compton.h"
#include "PROPOSAL/crossection/parametrization/EpairProduction.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"
#include "PROPOSAL/crossection/parametrization/MupairProduction.h"
#include "PROPOSAL/crossection/parametrization/ParamTables.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/crossection/parametrization/PhotoPairProduction.h"
#include "PROPOSAL/crossection/parametrization/PhotoQ2Integration.h"
#include "PROPOSAL/crossection/parametrization/PhotoRealPhotonAssumption.h"
#include "PROPOSAL/crossection/parametrization/Photonuclear.h"
#include "PROPOSAL/crossection/parametrization/WeakInteraction.h"

#include "PROPOSAL/crossection/CrossSectionDNDX/CrossSectionDNDX.h"
#include "PROPOSAL/crossection/CrossSectionDNDX/CrossSectionDNDXBuilder.h"
#include "PROPOSAL/crossection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"
#include "PROPOSAL/crossection/CrossSectionDNDX/CrossSectionDNDXInterpolant.h"

#include "PROPOSAL/secondaries/Parametrization.h"
#include "PROPOSAL/secondaries/SecondariesBuilder.h"
#include "PROPOSAL/secondaries/annihilation/Annihilation.h"
#include "PROPOSAL/secondaries/bremsstrahlung/Bremsstrahlung.h"
#include "PROPOSAL/secondaries/bremsstrahlung/NaivBremsstrahlung.h"
#include "PROPOSAL/secondaries/compton/Compton.h"
#include "PROPOSAL/secondaries/compton/NaivCompton.h"
#include "PROPOSAL/secondaries/ionization/Ionization.h"
#include "PROPOSAL/secondaries/ionization/NaivIonization.h"
#include "PROPOSAL/secondaries/mupairproduction/KelnerKokoulinPetrukhinMupairProduction.h"
#include "PROPOSAL/secondaries/mupairproduction/MupairProduction.h"

#include "PROPOSAL/secondaries/SecondariesCalculator.h"
#include "PROPOSAL/secondaries/photopairproduction/PhotoTsai.h"
#include "PROPOSAL/secondaries/photopairproduction/PhotopairProduction.h"

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/CrossSectionBuilder.h"
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

#include "PROPOSAL/crossection/ParticleDefaultCrossSectionList.h"

#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/decay/DecayTable.h"
#include "PROPOSAL/decay/LeptonicDecayChannel.h"
#include "PROPOSAL/decay/ManyBodyPhaseSpace.h"
#include "PROPOSAL/decay/StableChannel.h"
#include "PROPOSAL/decay/TwoBodyPhaseSpace.h"

#include "PROPOSAL/density_distr/density_distr.h"
#include "PROPOSAL/density_distr/density_exponential.h"
#include "PROPOSAL/density_distr/density_homogeneous.h"
#include "PROPOSAL/density_distr/density_polynomial.h"
#include "PROPOSAL/density_distr/density_splines.h"

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

#include "PROPOSAL/Propagator.h"

#include "PROPOSAL/propagation_utility/ContRand.h"
#include "PROPOSAL/propagation_utility/Decay.h"
#include "PROPOSAL/propagation_utility/DecayBuilder.h"
#include "PROPOSAL/propagation_utility/Displacement.h"
#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include "PROPOSAL/propagation_utility/Interaction.h"
#include "PROPOSAL/propagation_utility/InteractionBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include "PROPOSAL/propagation_utility/Time.h"
#include "PROPOSAL/propagation_utility/TimeBuilder.h"

#include "PROPOSAL/scattering/Coefficients.h"
#include "PROPOSAL/scattering/Scattering.h"
// #include "PROPOSAL/scattering/ScatteringFactory.h"
#include "PROPOSAL/scattering/ScatteringHighland.h"
#include "PROPOSAL/scattering/ScatteringHighlandIntegral.h"
#include "PROPOSAL/scattering/ScatteringMoliere.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/Secondaries.h"
#include "PROPOSAL/methods.h"

#if ROOT_SUPPORT
#include "PROPOSAL/interfaces/root.h"
#endif
