#include <string>

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/Secondaries.h"
#include "PROPOSAL/geometry/Geometry.h"
#include "pyPROPOSAL/pyBindings.h"

#define PARTICLE_DEF(module, cls)                                                    \
    py::class_<cls##Def, ParticleDef, std::shared_ptr<cls##Def>>(module, #cls "Def") \
        .def(py::init<>());                                                          \

namespace py = pybind11;
using namespace PROPOSAL;

void init_particle(py::module& m) {
    py::module m_sub = m.def_submodule("particle");

    m_sub.def("get_ParticleDef_for_type", &ParticleDef::GetParticleDefForType,
              py::arg("particle_type"), R"pbdoc(
                Get ParticleDef for given particle_type
    )pbdoc");

    m_sub.doc() = R"pbdoc(
        For each propagation a defined particle is needed.
        You have the possibility to define one by your own or select one of
        the listed below. Particle data are based on source ??.

        +----------+----------+----------+----------+----------+----------+
        |MuMinus   |MuPlus    |EMinus    |EPlus     |TauMinus  |TauPlus   |
        +----------+----------+----------+----------+----------+----------+
        |NuTau     |NuTauBar  |NuE       |NuEBar    |NuMu      |NuMuBar   |
        +----------+----------+----------+----------+----------+----------+
        |StauMinus |StauPlus  |Pi0       |PiMinus   |PiPlus    |K0        |
        +----------+----------+----------+----------+----------+----------+
        |KMinus    |KPlus     |Monopole  |Gamma     |SMPMinus  |SMPPlus   |
        +----------+----------+----------+----------+----------+----------+

        Particle are static objects, so that the properties can not be
        changed once it is initalized.  Own particle can be created if you
        initalize a complete new with the :meth:`ParticleDef` class or
        modifie a copy of a existing one with the :meth:`ParticleDefBuilder`.

        A predefined particle can be initialize for example with

        >>> muon = proposal.particle.MuMinusDef()

        The :meth:`Particle` class is a container for partilce related data.
        There can be for example initial energy set and losses read out.
        The :meth:`particle.Data` class is a container for interaction a
        particle does.
    )pbdoc";

    py::class_<ParticleDef, std::shared_ptr<ParticleDef>>(m_sub, "ParticleDef")
        // .def(py::init<>())
        .def(py::init<std::string, double, double, double, double,
                      const HardComponentTables::VecType&, const DecayTable&,
                      const int, const int>(),
             py::arg("name"), py::arg("mass"), py::arg("low"),
             py::arg("lifetime"), py::arg("charge"), py::arg("hard_component"),
             py::arg("decay_table"), py::arg("particle_type"), py::arg("weak_partner"))
        .def(py::init<const ParticleDef&>())

        .def("__str__", &py_print<ParticleDef>)
        .def("__eq__", &ParticleDef::operator==)
        .def("__ne__", &ParticleDef::operator!=)
        .def_readonly("name", &ParticleDef::name,
                      R"pbdoc(
                name ot the particle
            )pbdoc")
        .def_readonly("mass", &ParticleDef::mass,
                      R"pbdoc(
                mass of the particle in MeV
            )pbdoc")
        .def_readonly("low", &ParticleDef::low,
                      R"pbdoc(
                lower stochastic sampling energy limit.
            )pbdoc")
        .def_readonly("lifetime", &ParticleDef::lifetime,
                      R"pbdoc(
                average lifetime of the particle in seconds.
            )pbdoc")
        .def_readonly("charge", &ParticleDef::charge,
                      R"pbdoc(
                charge of the particle
            )pbdoc")
        .def_readonly("decay_table", &ParticleDef::decay_table,
                      R"pbdoc(
                generate a static particle with in ParticleDefBuilder properties
            )pbdoc")
        .def_readonly("particle_type", &ParticleDef::particle_type,
                      R"pbdoc(
                particle type of the particle
            )pbdoc")
        .def_readonly("weak_partner", &ParticleDef::weak_partner,
                      R"pbdoc(
                particle type of the weak partner particle
            )pbdoc")
        // .def_readonly("hard_component_table",
        // &ParticleDef::hard_component_table)
        // .add_property("hard_component_table",
        // make_function(&get_hard_component, return_internal_reference<>()))
        // //TODO(mario): shit Fri 2017/10/13
        ;

    py::class_<ParticleDef::Builder, std::shared_ptr<ParticleDef::Builder>>(
        m_sub, "ParticleDefBuilder",
        R"pbdoc(
                With the ParticleDefBuilder it is possible to define new particles
                or change porperties of existing one. First you have to initalize
                a ParticleDefBuilder. The new properties of the particle can be
                collected before a new static particle is build.

                Example:
                    A usecase for the ParticleDefBuilder might be if your
                    particle will be propagated stochastically until a
                    certain lower limit which is higher than the particle
                    mass is reached.

                    >>> builder = pp.particle.ParticleDefBuilder()
                    >>> builder.SetParticleDef(pp.particle.MuMinusDef.get())
                    >>> builder.SetLow(1e6)  # 1 Tev
                    >>> mu_def = mu_def_builder.build()

                    Therby muons which energy is lower than one TeV will be
                    handled continously.
            )pbdoc")
        .def(py::init<>(),
             R"pbdoc(
                Before you can create or modify a particle a definition builder has
                to be initalized. It collect the properties of the new or changed
                particle.
            )pbdoc")
        .def("SetName", &ParticleDef::Builder::SetName,
             R"pbdoc(
                Args:
                    arg1 (str): name of the particle
            )pbdoc")
        .def("SetMass", &ParticleDef::Builder::SetMass,
             R"pbdoc(
                Args:
                    arg1 (float): mass of the particle
            )pbdoc")
        .def("SetLow", &ParticleDef::Builder::SetLow,
             R"pbdoc(
                Args:
                    arg1 (float): lower stochastic sampling energy limit.
            )pbdoc")
        .def("SetLifetime", &ParticleDef::Builder::SetLifetime,
             R"pbdoc(
                Args:
                    arg1 (float): lifetime of the particle
            )pbdoc")
        .def("SetCharge", &ParticleDef::Builder::SetCharge,
             R"pbdoc(
                Args:
                    arg1 (float): charge of the particle in units of coulomb
            )pbdoc")
        .def("SetDecayTable", &ParticleDef::Builder::SetDecayTable,
             R"pbdoc(
                Args:
                    arg1 (???): ???
            )pbdoc")
        .def("SetParticleType", &ParticleDef::Builder::SetParticleType,
             R"pbdoc(
                Args:
                    arg1 (int): particle type of the particle
            )pbdoc")
        .def("SetWeakPartner", &ParticleDef::Builder::SetWeakPartner,
             R"pbdoc(
                Args:
                    arg1 (ParticleType): particle type of the weak partner particle
            )pbdoc")
        .def("SetParticleDef", &ParticleDef::Builder::SetParticleDef,
             R"pbdoc(
                Args:
                    arg1 (ParticleDef): a pre defined particle which values should
                        be take over.
            )pbdoc")
        .def("build", &ParticleDef::Builder::build,
             R"pbdoc(
                Return:
                    ParticleDef: generate a static particle with in ParticleDefBuilder properties
            )pbdoc");

    PARTICLE_DEF(m_sub, MuMinus)
    PARTICLE_DEF(m_sub, MuPlus)
    PARTICLE_DEF(m_sub, EMinus)
    PARTICLE_DEF(m_sub, EPlus)
    PARTICLE_DEF(m_sub, TauMinus)
    PARTICLE_DEF(m_sub, TauPlus)
    PARTICLE_DEF(m_sub, StauMinus)
    PARTICLE_DEF(m_sub, StauPlus)
    PARTICLE_DEF(m_sub, Pi0)
    PARTICLE_DEF(m_sub, PiMinus)
    PARTICLE_DEF(m_sub, PiPlus)
    PARTICLE_DEF(m_sub, K0)
    PARTICLE_DEF(m_sub, KMinus)
    PARTICLE_DEF(m_sub, KPlus)
    PARTICLE_DEF(m_sub, NuE)
    PARTICLE_DEF(m_sub, NuEBar)
    PARTICLE_DEF(m_sub, NuMu)
    PARTICLE_DEF(m_sub, NuMuBar)
    PARTICLE_DEF(m_sub, NuTau)
    PARTICLE_DEF(m_sub, NuTauBar)
    PARTICLE_DEF(m_sub, Monopole)
    PARTICLE_DEF(m_sub, Gamma)
    PARTICLE_DEF(m_sub, SMPMinus)
    PARTICLE_DEF(m_sub, SMPPlus)

    py::class_<ParticleState, std::shared_ptr<ParticleState>>(m_sub, "ParticleState",
                                                              R"pbdoc(
                Dynamic data objects store information about particle states,
                for example the initial particle state for the propagator,
                specific states during propagation or secondary particles.
            )pbdoc")
        .def(py::init<const Vector3D&, const Vector3D&, const double&, const double&, const double&>(),
                py::arg("position"), py::arg("direction"), py::arg("energy"),
                py::arg("time"), py::arg("propagated_distance"))
        .def(py::init<>())
        .def(py::init<const ParticleType&, const Vector3D&, const Vector3D&, const double&, const double&, const double&>(),
             py::arg("particle_type"), py::arg("position"), py::arg("direction"),
             py::arg("energy"), py::arg("time"), py::arg("propagated_distance"))
        .def(py::init<const ParticleState&>())
        .def("__str__", &py_print<ParticleState>)
        .def_readwrite("type", &ParticleState::type,
                               R"pbdoc(
                Type of Dynamic Data. Describes the ParticleType of instance.
            )pbdoc")
        .def_property_readonly("particle_def", &ParticleState::GetParticleDef,
                               R"pbdoc(
                Get corresponding particle_def to ParticleState.
            )pbdoc")
        .def_readwrite("position", &ParticleState::position,
                      R"pbdoc(
                Position of particle (in cm).
            )pbdoc")
        .def_readwrite("direction", &ParticleState::direction,
                      R"pbdoc(
                Direction of particle.
            )pbdoc")
        .def_readwrite("energy", &ParticleState::energy,
                      R"pbdoc(
                Energy of particle in MeV.
            )pbdoc")
        .def_property("momentum", &ParticleState::GetMomentum,
                      &ParticleState::SetMomentum,
                      R"pbdoc(
                Momentum of particle in MeV
            )pbdoc")
        .def_readwrite("time", &ParticleState::time,
                      R"pbdoc(
                Time (in sec) since beginning of propagation the particle.
            )pbdoc")
        .def_readwrite("propagated_distance", &ParticleState::propagated_distance,
                      R"pbdoc(
                Propagated distance of primary particle.
            )pbdoc");

    py::class_<Loss, std::shared_ptr<Loss>>(m_sub, "Loss", R"pbdoc(
                Loss objects are the base class for energy losses (stochastic
                and continuous losses) in PROPOSAL. They are defined by their
                type and energy. More specific information are stored in the
                derived classes.
            )pbdoc")
            .def(py::init<const int&, const double&, const double&>(),
                    py::arg("type"), py::arg("energy"), py::arg("parent_particle_energy"))
            .def_readwrite("type", &Loss::type, R"pbdoc(Type of energy loss.)pbdoc")
            .def_readwrite("energy", &Loss::energy, R"pbdoc(Total energy loss in MeV.)pbdoc")
            .def_readwrite("parent_particle_energy", &Loss::parent_particle_energy, R"pbdoc(Particle energy in MeV at the time of the stochastic loss or at the beginning of the continuous loss.)pbdoc");

    py::class_<StochasticLoss, Loss, std::shared_ptr<StochasticLoss>>(m_sub, "StochasticLoss")
            .def(py::init<const int&, const double&, const Vector3D&, const Vector3D&, const double&, const double&, const double &, const size_t&>(),
                 py::arg("type"), py::arg("loss_energy"), py::arg("position"),
                 py::arg("direction"), py::arg("time"),
                 py::arg("propagated_distance"),
                 py::arg("parent_particle_energy"),
                 py::arg("target_hash") = 0)
            .def_readwrite("position", &StochasticLoss::position, R"pbdoc(Position of stochastic interaction.)pbdoc")
            .def_readwrite("direction", &StochasticLoss::direction, R"pbdoc(Direction of stochastic loss.)pbdoc")
            .def_readwrite("time", &StochasticLoss::time, R"pbdoc(Time when stochastic loss occured.)pbdoc")
            .def_readwrite("target_hash", &StochasticLoss::target_hash, R"pbdoc(Hash of target with which the stochastic interaction happened. Can be the hash of a Component or the hash of a Medium`.)pbdoc")
            .def_readwrite("propagated_distance", &StochasticLoss::propagated_distance, R"pbdoc(Distance (in cm) the parent particle has propagated when the stochastic loss occured.)pbdoc");

    py::class_<ContinuousLoss, Loss, std::shared_ptr<ContinuousLoss>>(m_sub, "ContinuousLoss")
            .def(py::init<const double&, const double&, const Vector3D&, const Vector3D&, const Vector3D&, const Vector3D&, const double&, const double&>(),
                 py::arg("energy"), py::arg("parent_particle_energy"), py::arg("start_position"), py::arg("end_position"),
                 py::arg("direction_initial"), py::arg("direction_final"), py::arg("time_initial"), py::arg("time_final"))
            .def_readwrite("start_position", &ContinuousLoss::start_position, R"pbdoc(Position where the continuous energy loss started.)pbdoc")
            .def_readwrite("end_position", &ContinuousLoss::end_position, R"pbdoc(Position where the continuous energy loss ended.)pbdoc")
            .def_readwrite("direction_initial", &ContinuousLoss::direction_initial, R"pbdoc(Direction of the particle at the beginning of the continuous energy loss.)pbdoc")
            .def_readwrite("direction_final", &ContinuousLoss::direction_final, R"pbdoc(Direction of the particle at the end of the continuous energy loss.)pbdoc")
            .def_readwrite("time_initial", &ContinuousLoss::time_initial, R"pbdoc(Time when the continuous energy loss started.)pbdoc")
            .def_readwrite("time_final", &ContinuousLoss::time_final, R"pbdoc(Time when the continuous energy loss ended.)pbdoc");

    py::class_<Secondaries, std::shared_ptr<Secondaries>>(m_sub, "Secondaries",
            R"pbdoc(Output of the propagator. The secondaries class stores all information about interactions during the propagation. For this, it stores all
                    intermediate particle states during propagation in a track object, which is a list of
                    :meth:`particle.ParticleState` objects.
            )pbdoc")
            .def("initial_state",
                 &Secondaries::GetInitialState,
                 R"pbdoc(Get information on the initial state of the propagated particle.)pbdoc")
            .def("final_state",
                 &Secondaries::GetFinalState,
                 R"pbdoc(Get information on the final state of the propagated particle.)pbdoc")
            .def("ELost",
                 &Secondaries::GetELost, py::arg("geometry"),
                 R"pbdoc(

                Calculate energy that the particle has lost inside a specific geometry.

                Args:
                    geometry (Geometry): Lost energy in MeV
                Returns:
                    Lost energy in MeV

                )pbdoc")
            .def("entry_point",
                 &Secondaries::GetEntryPoint, py::arg("geometry"),
                 R"pbdoc(
                Get the particle state when the particle has entered a specific geometry. If the particle already
                started in the specific geometry or has never entered the specific geometry, a None is returned.
                Otherwise, a continuous step based on a known particle state will be calculated so that the returned
                particle state describes the state when the particle entered the geometry. Note that this only works
                reliably when continuous randomization is disabled. Otherwise, a warning will be prompted and results
                might be inconsistent.

                Args:
                    geometry (Geometry): Geometry for which we want to check when the propagated particle has entered it

                Returns:
                    ParticleState, describing the state of the particle when it has entered the specific geometry, or a
                    None.
                )pbdoc")
            .def("exit_point",
                 &Secondaries::GetExitPoint, py::arg("geometry"),
                 R"pbdoc(
                Get the particle state when the particle has left a specific geometry. If the particle has never left
                or never crossed the geometry, a None is returned. Otherwise, a continuous step based on a known
                particle state will be calculated so that the returned particle state describes the state when the
                particle exited the geometry. Note that this only works reliably when continuous randomization is
                disabled. Otherwise, a warning will be prompted and results might be inconsistent.

                Args:
                    geometry (Geometry): Geometry for which we want to check when the propagated particle has left it

                Returns:
                    ParticleState, describing the state of the particle when it has left the specific geometry, or a
                    None.
                 )pbdoc")
            .def("closest_approach_point",
                 &Secondaries::GetClosestApproachPoint, py::arg("geometry"),
                 R"pbdoc(
                Get the particle state when the particle had the smallest distance to the center of a specific geometry.
                This will not extrapolate the particle before the initial or after the final state. This will calculate
                a continuous step, based on a known particle state, so that the returned particle state describes the
                particle at is has the smallest distance to the center of the specified geometry.

                Args:
                    geometry (Geometry): Geometry for which we want to calculate the closest approach point
                Returns:
                    A ParticleState, describing the closest approach point of the particle.
                )pbdoc")
            .def("hit_geometry",
                 &Secondaries::HitGeometry, py::arg("geometry"),
                 R"pbdoc(
                Check if particle has hit a specific geometry. Particle tracks ending at the border of a geometry
                count as a hit.
                Args:
                    geometry (Geometry): Geometry for which we want to check whether our particle has hit it
                Returns:
                    Bool - True if geometry has been hit, false otherwise.
                )pbdoc")
            .def("track",
                 overload_cast_<>()(&Secondaries::GetTrack, py::const_),
                 R"pbdoc(
                Returns particle track. This is a list of particle states that describe the particle during propagation.
                The first element of this list describes the initial particle state, the last element of this list
                describes the particle at the end of propagation. All intermediate elements of the list describe the
                particle during propagation. This includes particle states immediately before and after stochastic
                losses and particle states when there has been a sector transition.

                Returns:
                    List of ParticleState objects describing the states of the particle during propagation.
                )pbdoc")
            .def("track",
                 overload_cast_<const Geometry&>()(&Secondaries::GetTrack, py::const_), py::arg("geometry"), R"pbdoc()pbdoc")
            .def("get_state_for_energy",
                 &Secondaries::GetStateForEnergy, py::arg("energy"),
                 R"pbdoc(
                Get information on the state of the propagated particle when it has reached a specific energy.

                If the specified energy is within a continuous loss, a continuous propagation step will be simulated so
                that the energy of the returned particle state is equal to the specified energy. Note that this only
                works consistently if continuous randomization is disabled. Otherwise, a warning will be prompted and
                results might be inconsistent.

                If the specified energy lies within a stochastic loss (i.e. the particle energy has been bigger than
                the specified energy before a stochastic loss but smaller than the specified energy after a stochastic
                loss), the particle state before the stochastic loss is returned.

                If the specified energy is below the final energy of the particle, the final particle state will be returned.

                If the specified energy is above the initial energy of the particle, the initial particle state will be returned.

                Args:
                    energy: Particle energy (in MeV and in total energies) for which we want to get the particle state

                Returns:
                    ParticleState object, describing the particle state at energy
                )pbdoc")
            .def("get_state_for_distance",
                 &Secondaries::GetStateForDistance, py::arg("distance"),
                 R"pbdoc(
                Get information on the state of the propagated particle when it has reaches a specific propagation
                distance.

                This will calculate a continuous propagation step based on a known particle state so that the propagated
                distance of the returned particle state is equal to the specified propagated distance.

                If the specified distance is below the propagated distance of the initial particle state, the initial
                particle state is returned.

                If the specified distance is above the propagated distance of the final
                particle state, the final particle state is returned.

                Note that by default, the particle starts with a propagation_distance of 0 unless specified otherwise in
                the initial_particle state.

                Args:
                    propagated_distance: Propagated distance (in cm) for which we want to get the particle state

                Returns:
                    ParticleState object, describing the particle state when it has propagated the propagated_distance
                )pbdoc")
            .def("track_positions",
                 &Secondaries::GetTrackPositions,
                 R"pbdoc(
                Returns the list of all positions of the propagated particle. The first element corresponds to the
                position of the initial particle, the last element to the final position of the particle after propagation.
                The elements in between describe the intermediate positions of the particle during propagation.

                 Returns:
                    List of Cartesian3D objects, describing the particle positions during propagation. Note that the
                    units of the Cartesian3D object are in cm.
                )pbdoc")
            .def("track_directions",
                 &Secondaries::GetTrackDirections,
                 R"pbdoc(
                Returns the list of all directions of the propagated particle. The first element corresponds to the
                direction of the initial particle, the last element to the final direction of the particle after
                propagation. The elements in between describe the intermediate directions of the particle during
                propagation.

                Returns:
                    List of Cartesian3D objects, describing the particle directions during propagation.
                )pbdoc")
            .def("track_energies",
                 &Secondaries::GetTrackEnergies,
                 R"pbdoc(
                Returns the list of all particle energies (total energies in MeV) of the propagated particle. The first
                element corresponds to the energy of the initial particle, the last element to the final energy of the
                particle after propagation. The elements in between describe the energies of the particle during
                propagation.

                Returns:
                    List of doubles, describing the particle energies during propagation (in total energies, in MeV)
                )pbdoc")
            .def("track_times",
                 &Secondaries::GetTrackTimes,
                 R"pbdoc(
                Returns the list of all particle times (times in s) of the propagated particle. The first element
                corresponds to the time of the initial particle, the last element to the time of the particle after
                propagation. The elements in between describe the times of the particle during propagation.

                Returns:
                    List of doubles, describing the particle times during propagation (in s)

                )pbdoc")
            .def("track_propagated_distances",
                 &Secondaries::GetTrackPropagatedDistances,
                 R"pbdoc(
                Returns the list of the propagated distances (in cm) of the propagated particle. The first element
                orresponds to the propagated distance of the initial particle (0cm if not specified otherwise in the
                initial state), the last element to the propagated distance of the particle after propagation. The
                elements in between describe the propagated distances during propagation

                Returns:
                    List of doubles, describing the particle propagated distances during propgation (in cm)
                 )pbdoc")
            .def("track_types",
                 &Secondaries::GetTrackTypes,
                 R"pbdoc(
                Returns a list of interaction types describing the interactions of the particle during propagation.

                Returns:
                    List of InteractionType objects
                )pbdoc")
            .def("track_length",
                 &Secondaries::GetTrackLength,
                 R"pbdoc(
                Length of the track objects, i.e. the number of particle states that are stored in this Secondaries class.

                Returns:
                    unsigned int which describes the length of the track object.
                 )pbdoc")
            .def("target_hashes",
                 &Secondaries::GetTargetHashes,
                 R"pbdoc(
                Returns a list of hashes, describing the media and components that the particle interacted with during
                propagation. For continuous losses, 0 is returned (since the particle did not interact with a specific
                medium).

                For ionization losses, the hash of the medium that the particle interacted with is returned. We can get
                the medium to the corresponding hash using :meth:`~proposal.medium.get_medium_for_hash`.

                For other stochastic losses, the hash of the component that the particle interacted with is returned.
                We can get the component to the corresponding hash using :meth:`~proposal.component.get_component_for_hash`.

                Returns:
                    List of size_t objects, which are hashed to Medium or Component objects. These are the targets that our particle interacted with during propagation.
                )pbdoc")
            .def("stochastic_losses",
                 overload_cast_<>()(&Secondaries::GetStochasticLosses, py::const_),
                 R"pbdoc(
                Get all stochastic losses (i.e. losses bigger than the energy cut) of the particle that occurred during
                propagation

                Returns:
                    List of StochasticLoss objects
                )pbdoc")
            .def("stochastic_losses",
                 overload_cast_<const Geometry&>()(&Secondaries::GetStochasticLosses, py::const_), py::arg("geometry"),
                 R"pbdoc(
                Get stochastic losses (i.e. losses bigger than the energy cut) of the particle that are within a
                specified geometry.

                Args:
                    geometry: Geometry object, find stochastic losses within this geometry

                Returns:
                    List of StochasticLoss objects
                )pbdoc")
            .def("stochastic_losses",
                 overload_cast_<const InteractionType&>()(&Secondaries::GetStochasticLosses, py::const_), py::arg("interaction_type"),
                 R"pbdoc(
                Get stochastic losses (i.e. losses bigger than the energy cut) of the particle that correspond to a
                specific interaction type

                Args:
                    interaction_type: Interaction type of stochastic loss as an InteractionType object

                Returns:
                    List of StochasticLoss objects
                )pbdoc")
            .def("stochastic_losses",
                 overload_cast_<const std::string&>()(&Secondaries::GetStochasticLosses, py::const_), py::arg("interaction_type"),
                 R"pbdoc(
                Get stochastic losses (i.e. losses bigger than the energy cut) of the particle that correspond to a
                specific interaction type

                Args:
                    interaction_type: Name of interaction type of stochastic loss

                Returns:
                    List of StochasticLoss objects
                )pbdoc")
            .def("continuous_losses",
                 overload_cast_<>()(&Secondaries::GetContinuousLosses, py::const_),
                 R"pbdoc(
                Get all continuous losses (i.e. losses smaller than the energy cut) of the particle that occurred
                during propagation

                Returns:
                    List of ContinousLoss objects
                )pbdoc")
            .def("continuous_losses",
                 overload_cast_<const Geometry&>()(&Secondaries::GetContinuousLosses, py::const_), py::arg("geometry"),
                 R"pbdoc(
                Get all continuous losses (i.e. losses smaller than the energy cut) of the particle that occurred
                during propagation

                Args:
                    geometry: Geometry object, find continuous losses within this geometry
                )pbdoc")
            .def("decay_products",
                 &Secondaries::GetDecayProducts,
                 R"pbdoc(
                If the particle has decayed at the end of propagation, this function calculated the decay products as a
                list of particle states. If the particle did not decay during propagation, the returned list will be
                empty.

                Returns:
                    List of ParticleStates, describing the decay products.
                )pbdoc");

    py::enum_<InteractionType>(m_sub, "Interaction_Type")
        .value("particle", InteractionType::Particle)
        .value("brems", InteractionType::Brems)
        .value("ioniz", InteractionType::Ioniz)
        .value("epair", InteractionType::Epair)
        .value("photonuclear", InteractionType::Photonuclear)
        .value("mupair", InteractionType::MuPair)
        .value("hadrons", InteractionType::Hadrons)
        .value("continuousenergyloss", InteractionType::ContinuousEnergyLoss)
        .value("weakint", InteractionType::WeakInt)
        .value("compton", InteractionType::Compton)
        .value("decay", InteractionType::Decay)
        .value("annihilation", InteractionType::Annihilation)
        .value("photopair", InteractionType::Photopair)
        .value("photoproduction", InteractionType::Photoproduction)
        .value("photomupair", InteractionType::PhotoMuPair)
        .value("photoeffect", InteractionType::Photoeffect);

   py::enum_<ParticleType>(m_sub, "Particle_Type")
        .value("None", ParticleType::None)
        .value("EMinus", ParticleType::EMinus)
        .value("EPlus", ParticleType::EPlus)
        .value("NuE", ParticleType::NuE)
        .value("NuEBar", ParticleType::NuEBar)
        .value("MuMinus", ParticleType::MuMinus)
        .value("NuMu", ParticleType::NuMu)
        .value("NuMuBar", ParticleType::NuMuBar)
        .value("MuPlus", ParticleType::MuPlus)
        .value("TauMinus", ParticleType::TauMinus)
        .value("TauPlus", ParticleType::TauPlus)
        .value("NuTau", ParticleType::NuTau)
        .value("NuTauBar", ParticleType::NuTauBar)
        .value("Gamma", ParticleType::Gamma)
        .value("Pi0", ParticleType::Pi0)
        .value("PiPlus", ParticleType::PiPlus)
        .value("PiMinus", ParticleType::PiMinus)
        .value("K0", ParticleType::K0)
        .value("KPlus", ParticleType::KPlus)
        .value("KMinus", ParticleType::KMinus)
        .value("STauMinus", ParticleType::STauMinus)
        .value("STauPlus", ParticleType::STauPlus)
        .value("Monopole", ParticleType::Monopole)
        .value("Hadron", ParticleType::Hadron);
}

#undef PARTICLE_DEF
