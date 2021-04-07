
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
#include "PROPOSAL/version.h"

#include "PROPOSAL/geometry/Box.h"
#include "PROPOSAL/geometry/Cylinder.h"
#include "PROPOSAL/geometry/GeometryFactory.h"
#include "PROPOSAL/geometry/Sphere.h"

#include "PROPOSAL/crosssection/parametrization/Annihilation.h"
#include "PROPOSAL/crosssection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/crosssection/parametrization/Compton.h"
#include "PROPOSAL/crosssection/parametrization/EpairProduction.h"
#include "PROPOSAL/crosssection/parametrization/Ionization.h"
#include "PROPOSAL/crosssection/parametrization/MupairProduction.h"
#include "PROPOSAL/crosssection/parametrization/ParamTables.h"
#include "PROPOSAL/crosssection/parametrization/PhotoPairProduction.h"
#include "PROPOSAL/crosssection/parametrization/PhotoQ2Integration.h"
#include "PROPOSAL/crosssection/parametrization/PhotoRealPhotonAssumption.h"
#include "PROPOSAL/crosssection/parametrization/Photonuclear.h"
#include "PROPOSAL/crosssection/parametrization/WeakInteraction.h"

#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXBuilder.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXInterpolant.h"

#include "PROPOSAL/crosssection/Factories/AnnihilationFactory.h"
#include "PROPOSAL/crosssection/Factories/BremsstrahlungFactory.h"
#include "PROPOSAL/crosssection/Factories/ComptonFactory.h"
#include "PROPOSAL/crosssection/Factories/EpairProductionFactory.h"
#include "PROPOSAL/crosssection/Factories/IonizationFactory.h"
#include "PROPOSAL/crosssection/Factories/MupairProductionFactory.h"
#include "PROPOSAL/crosssection/Factories/PhotoPairProductionFactory.h"
#include "PROPOSAL/crosssection/Factories/PhotonuclearFactory.h"
#include "PROPOSAL/crosssection/Factories/WeakInteractionFactory.h"

#include "PROPOSAL/secondaries/parametrization/Parametrization.h"
#include "PROPOSAL/secondaries/parametrization/annihilation/Annihilation.h"
#include "PROPOSAL/secondaries/parametrization/annihilation/HeitlerAnnihilation.h"
#include "PROPOSAL/secondaries/parametrization/annihilation/SingleDifferentialAnnihilation.h"
#include "PROPOSAL/secondaries/parametrization/bremsstrahlung/Bremsstrahlung.h"
#include "PROPOSAL/secondaries/parametrization/bremsstrahlung/NaivBremsstrahlung.h"
#include "PROPOSAL/secondaries/parametrization/compton/Compton.h"
#include "PROPOSAL/secondaries/parametrization/compton/NaivCompton.h"
#include "PROPOSAL/secondaries/parametrization/epairproduction/EpairProduction.h"
#include "PROPOSAL/secondaries/parametrization/epairproduction/KelnerKokoulinPetrukhinEpairProduction.h"
#include "PROPOSAL/secondaries/parametrization/epairproduction/NaivEpairProduction.h"
#include "PROPOSAL/secondaries/parametrization/ionization/Ionization.h"
#include "PROPOSAL/secondaries/parametrization/ionization/NaivIonization.h"
#include "PROPOSAL/secondaries/parametrization/mupairproduction/KelnerKokoulinPetrukhinMupairProduction.h"
#include "PROPOSAL/secondaries/parametrization/mupairproduction/MupairProduction.h"
#include "PROPOSAL/secondaries/parametrization/photonuclear/Photonuclear.h"
#include "PROPOSAL/secondaries/parametrization/photopairproduction/PhotoTsai.h"
#include "PROPOSAL/secondaries/parametrization/photopairproduction/PhotoTsaiForwardPeaked.h"
#include "PROPOSAL/secondaries/parametrization/photopairproduction/PhotopairProduction.h"

#include "PROPOSAL/secondaries/SecondariesCalculator.h"

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"
#include "PROPOSAL/crosssection/CrossSectionIntegral.h"
#include "PROPOSAL/crosssection/CrossSectionInterpolant.h"
#include "PROPOSAL/crosssection/CrossSectionMultiplier.h"

#include "PROPOSAL/crosssection/ParticleDefaultCrossSectionList.h"

#include "PROPOSAL/decay/DecayTable.h"
#include "PROPOSAL/decay/LeptonicDecayChannel.h"
#include "PROPOSAL/decay/ManyBodyPhaseSpace.h"
#include "PROPOSAL/decay/StableChannel.h"
#include "PROPOSAL/decay/TwoBodyPhaseSpace.h"

#include "PROPOSAL/density_distr/density_exponential.h"
#include "PROPOSAL/density_distr/density_homogeneous.h"
#include "PROPOSAL/density_distr/density_polynomial.h"
#include "PROPOSAL/density_distr/density_splines.h"

#include "PROPOSAL/math/Cartesian3D.h"
#include "PROPOSAL/math/Function.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/math/Spherical3D.h"
#include "PROPOSAL/math/Spline.h"
#include "PROPOSAL/math/TableWriter.h"

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/particle/ParticleDef.h"

#include "PROPOSAL/Propagator.h"

#include "PROPOSAL/propagation_utility/ContRand.h"
#include "PROPOSAL/propagation_utility/ContRandBuilder.h"
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

#include "PROPOSAL/scattering/stochastic_deflection/ScatteringFactory.h"
#include "PROPOSAL/scattering/stochastic_deflection/bremsstrahlung/Bremsstrahlung.h"
#include "PROPOSAL/scattering/stochastic_deflection/bremsstrahlung/TsaiApproximationBremsstrahlung.h"
#include "PROPOSAL/scattering/stochastic_deflection/nuclearInteraction/NuclearInteraction.h"
#include "PROPOSAL/scattering/stochastic_deflection/nuclearInteraction/BorogPetrukhinNuclearInteraction.h"
#include "PROPOSAL/scattering/stochastic_deflection/pairProd/PairProduction.h"
#include "PROPOSAL/scattering/stochastic_deflection/pairProd/KelnerPairProduction.h"
#include "PROPOSAL/scattering/stochastic_deflection/ionization/Ionization.h"
#include "PROPOSAL/scattering/stochastic_deflection/ionization/NaivIonization.h"

#include "PROPOSAL/scattering/multiple_scattering/Coefficients.h"
#include "PROPOSAL/scattering/multiple_scattering/Highland.h"
#include "PROPOSAL/scattering/multiple_scattering/HighlandIntegral.h"
#include "PROPOSAL/scattering/multiple_scattering/Moliere.h"
#include "PROPOSAL/scattering/multiple_scattering/ScatteringFactory.h"

#include "PROPOSAL/scattering/Scattering.h"
#include "PROPOSAL/scattering/ScatteringFactory.h"
#include "PROPOSAL/scattering/ScatteringMultiplier.h"
#include "PROPOSAL/scattering/ScatteringFactory.h"
#include "PROPOSAL/scattering/DirectionChange.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/Secondaries.h"
#include "PROPOSAL/methods.h"

#if ROOT_SUPPORT
#include "PROPOSAL/interfaces/root.h"
#endif
