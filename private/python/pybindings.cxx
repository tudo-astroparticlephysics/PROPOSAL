#define BOOST_PYTHON_MAX_ARITY 17

// #include <string>

#include <boost/python.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "PROPOSAL/Propagator.h"
// #include "PROPOSAL/PROPOSALParticle.h"
// #include "PROPOSAL/Medium.h"
// #include "PROPOSAL/EnergyCutSettings.h"
// #include "PROPOSAL/ProcessCollection.h"
// #include "PROPOSAL/CrossSections.h"
// #include "PROPOSAL/Photonuclear.h"
// #include "PROPOSAL/Bremsstrahlung.h"
// #include "PROPOSAL/Epairproduction.h"
// #include "PROPOSAL/Ionization.h"
// #include "PROPOSAL/Geometry.h"

using namespace PROPOSAL;

// struct CrossSectionsWrap : CrossSections, boost::python::wrapper<CrossSections>
// {
//     CrossSectionsWrap(): CrossSections(), wrapper<CrossSections>(){}
//
//     // virtual double CalculatedEdx()
//     // {
//     //     return this->get_override("calculate_dEdx")();
//     // }
// };


/******************************************************************************
*                             Register functions                              *
******************************************************************************/

// ------------------------------------------------------------------------- //
// Vectors
// ------------------------------------------------------------------------- //


struct CrossSectionToPython
{
    static PyObject* convert(std::vector<CrossSections*> const& vec)
    {
        boost::python::list python_list;

        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            switch (vec[i]->GetType())
            {
                case ParticleType::Brems:
                    python_list.append(boost::python::ptr((Bremsstrahlung*)vec[i]));
                    break;
                case ParticleType::NuclInt:
                    python_list.append(boost::python::ptr((Photonuclear*)vec[i]));
                    break;
                case ParticleType::DeltaE:
                    python_list.append(boost::python::ptr((Ionization*)vec[i]));
                    break;
                case ParticleType::EPair:
                    python_list.append(boost::python::ptr((Epairproduction*)vec[i]));
                    break;
                default:
                    boost::python::object obj(NULL);
                    break;
            }
        }

        PyObject* py = boost::python::incref(python_list.ptr());
        return py;
    }
};

template<typename T>
struct VectorToPythonList
{
    static PyObject* convert(std::vector<T> const& vec)
    {
        boost::python::list python_list;
        typename std::vector<T>::const_iterator iter;

        for(iter = vec.begin(); iter != vec.end(); ++iter)
        {
            python_list.append(boost::python::object(*iter));
        }

        return boost::python::incref(python_list.ptr());
    }
};

template<typename T>
struct PVectorToPythonList
{
    static PyObject* convert(std::vector<T> const& vec)
    {
        boost::python::list python_list;
        typename std::vector<T>::const_iterator iter;

        for(iter = vec.begin(); iter != vec.end(); ++iter)
        {
            python_list.append(boost::python::ptr(*iter));
        }

        return boost::python::incref(python_list.ptr());
    }
};


template<typename T>
struct VectorFromPythonList
{

    VectorFromPythonList()
    {
        boost::python::converter::registry::push_back(&VectorFromPythonList<T>::convertible,
                            &VectorFromPythonList<T>::construct,
                            boost::python::type_id<std::vector<T> >());
    }

    // Determine if obj_ptr can be converted in a std::vector<T>
    static void* convertible(PyObject* obj_ptr)
    {
        if (!PyList_Check(obj_ptr))
        {
            return 0;
        }

        return obj_ptr;
    }

    // Convert obj_ptr into a std::vector<T>
    static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        // Use borrowed to construct the object so that a reference
        // count will be properly handled.
        boost::python::list python_list(boost::python::handle<>(boost::python::borrowed(obj_ptr)));

        // Grab pointer to memory into which to construct the new std::vector<T>
        void* storage = reinterpret_cast<boost::python::converter::rvalue_from_python_storage<std::vector<T> >*>(data)->storage.bytes;

        // in-place construct the new std::vector<T> using the character data
        // extraced from the python object
        std::vector<T>& vec = *(new (storage) std::vector<T>());

        // populate the vector from list contains !!!
        int lenght = boost::python::len(python_list);
        vec.resize(lenght);

        for(int i = 0; i != lenght; ++i)
        {
            vec[i] = boost::python::extract<T>(python_list[i]);
        }

        // Stash the memory chunk pointer for later use by boost.python
        data->convertible = storage;
    }
};

// ------------------------------------------------------------------------- //
// Pair
// ------------------------------------------------------------------------- //


template<typename T1, typename T2>
struct PairToPythonList
{
    static PyObject* convert(std::pair<T1, T2> const& p)
    {
        boost::python::list python_list;

        python_list.append(boost::python::object(p.first));
        python_list.append(boost::python::object(p.second));
        // typename std::pair<T1, T2>::const_iterator iter;

        // for(iter = vec.begin(); iter != vec.end(); ++iter)
        // {
        //     python_list.append(boost::python::object(*iter));
        // }

        return boost::python::incref(python_list.ptr());
    }
};


template<typename T1, typename T2>
struct PairFromPythonList
{

    PairFromPythonList()
    {
        boost::python::converter::registry::push_back(&PairFromPythonList<T1, T2>::convertible,
                            &PairFromPythonList<T1, T2>::construct,
                            boost::python::type_id<std::pair<T1, T2> >());
    }

    // Determine if obj_ptr can be converted in a std::pair
    static void* convertible(PyObject* obj_ptr)
    {
        if (!PyList_Check(obj_ptr))
        {
            return 0;
        }

        return obj_ptr;
    }

    // Convert obj_ptr into a std::pair<T>
    static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        // Use borrowed to construct the object so that a reference
        // count will be properly handled.
        boost::python::list python_list(boost::python::handle<>(boost::python::borrowed(obj_ptr)));

        // Grab pointer to memory into which to construct the new std::pair<T1, T2>
        void* storage = reinterpret_cast<boost::python::converter::rvalue_from_python_storage<std::pair<T1, T2> >*>(data)->storage.bytes;

        // in-place construct the new std::pair<T1, T2> using the character data
        // extraced from the python object
        std::pair<T1, T2>& p = *(new (storage) std::pair<T1, T2>());

        assert(boost::python::len(python_list == 2));

        // Populate pari from python list
        p.first = boost::python::extract<T1>(python_list[0]);
        p.second = boost::python::extract<T2>(python_list[1]);

        // Stash the memory chunk pointer for later use by boost.python
        data->convertible = storage;
    }
};

/******************************************************************************
*                               Python Module                                 *
******************************************************************************/

BOOST_PYTHON_MODULE(pyPROPOSAL)
{

    using namespace boost::python;

    docstring_options doc_options(true, true, false);

    // --------------------------------------------------------------------- //
    // Vector classes
    // --------------------------------------------------------------------- //

    // register the to-python converter
    to_python_converter< std::vector<double>, VectorToPythonList<double> >();
    to_python_converter< std::vector<std::string>, VectorToPythonList<std::string> >();
    to_python_converter< std::vector<PROPOSALParticle*>, PVectorToPythonList<PROPOSALParticle*> >();
    to_python_converter< std::vector<ProcessCollection*>, PVectorToPythonList<ProcessCollection*> >();

    to_python_converter<std::vector<CrossSections*>, CrossSectionToPython>();

    to_python_converter<std::pair<double, double>, PairToPythonList<double, double> >();

    // register the from-python converter
    VectorFromPythonList<double>();
    VectorFromPythonList<std::string>();
    VectorFromPythonList<PROPOSALParticle*>();
    VectorFromPythonList<CrossSections*>();
    VectorFromPythonList<ProcessCollection*>();

    PairFromPythonList<double, double>();

    // class_<std::vector<CrossSections*> >("CrossSections")
    //         .def(vector_indexing_suite<std::vector<CrossSections*> >())
    //     ;

    // --------------------------------------------------------------------- //
    // ParticleType::Enum
    // --------------------------------------------------------------------- //

    enum_<ParticleType::Enum>("ParticleType")
        .value("EPlus",                 ParticleType::EPlus)
        .value("EMinus",                ParticleType::EPlus)
        .value("MuPlus",                ParticleType::MuPlus)
        .value("MuMinus",               ParticleType::MuMinus)
        .value("TauPlus",               ParticleType::TauPlus)
        .value("TauMinus",              ParticleType::TauMinus)
        .value("NuE",                   ParticleType::NuE)
        .value("NuEBar",                ParticleType::NuEBar)
        .value("NuMu",                  ParticleType::NuMu)
        .value("NuMuBar",               ParticleType::NuMuBar)
        .value("NuTau",                 ParticleType::NuTau)
        .value("NuTauBar",              ParticleType::NuTauBar)
        .value("Brems",                 ParticleType::Brems)
        .value("DeltaE",                ParticleType::DeltaE)
        .value("EPair",                 ParticleType::EPair)
        .value("NuclInt",               ParticleType::NuclInt)
        .value("MuPair",                ParticleType::MuPair)
        .value("Hadrons",               ParticleType::Hadrons)
        .value("ContinuousEnergyLoss",  ParticleType::ContinuousEnergyLoss)
        .value("Monopole",              ParticleType::Monopole)
        .value("STauPlus",              ParticleType::STauPlus)
        .value("STauMinus",             ParticleType::STauMinus)
        .value("StableMassiveParticle", ParticleType::StableMassiveParticle)
    ;

    // --------------------------------------------------------------------- //
    // ParametrizationType
    // --------------------------------------------------------------------- //

    enum_<ParametrizationType::Enum>("ParametrizationType")
        .value("BremsKelnerKokoulinPetrukhin",                  ParametrizationType::BremsKelnerKokoulinPetrukhin)
        .value("BremsAndreevBezrukovBugaev",                    ParametrizationType::BremsAndreevBezrukovBugaev)
        .value("BremsPetrukhinShestakov",                       ParametrizationType::BremsPetrukhinShestakov)
        .value("BremsCompleteScreeningCase",                    ParametrizationType::BremsCompleteScreeningCase)
        .value("PhotoKokoulinShadowBezrukovSoft",               ParametrizationType::PhotoKokoulinShadowBezrukovSoft)
        .value("PhotoKokoulinShadowBezrukovHard",               ParametrizationType::PhotoKokoulinShadowBezrukovHard)
        .value("PhotoRhodeShadowBezrukovSoft",                  ParametrizationType::PhotoRhodeShadowBezrukovSoft)
        .value("PhotoRhodeShadowBezrukovHard",                  ParametrizationType::PhotoRhodeShadowBezrukovHard)
        .value("PhotoBezrukovBugaevShadowBezrukovSoft",         ParametrizationType::PhotoBezrukovBugaevShadowBezrukovSoft)
        .value("PhotoBezrukovBugaevShadowBezrukovHard",         ParametrizationType::PhotoBezrukovBugaevShadowBezrukovHard)
        .value("PhotoZeusShadowBezrukovSoft",                   ParametrizationType::PhotoZeusShadowBezrukovSoft)
        .value("PhotoZeusShadowBezrukovHard",                   ParametrizationType::PhotoZeusShadowBezrukovHard)
        .value("PhotoAbramowiczLevinLevyMaor91ShadowDutta",     ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowDutta)
        .value("PhotoAbramowiczLevinLevyMaor91ShadowButkevich", ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowButkevich)
        .value("PhotoAbramowiczLevinLevyMaor97ShadowDutta",     ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowDutta)
        .value("PhotoAbramowiczLevinLevyMaor97ShadowButkevich", ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich)
        .value("PhotoRenoSarcevicSuShadowDutta",                ParametrizationType::PhotoRenoSarcevicSuShadowDutta)
        .value("PhotoRenoSarcevicSuShadowButkevich",            ParametrizationType::PhotoRenoSarcevicSuShadowButkevich)
        .value("EPairKelnerKokoulinPetrukhin",                  ParametrizationType::EPairKelnerKokoulinPetrukhin)
        .value("IonizBetheBloch",                               ParametrizationType::IonizBetheBloch)
    ;

    // --------------------------------------------------------------------- //
    // MediumType
    // --------------------------------------------------------------------- //

    enum_<MediumType::Enum>("MediumType")
        .value("Water"         , MediumType::Water)
        .value("Ice"           , MediumType::Ice)
        .value("Hydrogen"      , MediumType::Hydrogen)
        .value("Iron"          , MediumType::Iron)
        .value("Copper"        , MediumType::Copper)
        .value("Lead"          , MediumType::Lead)
        .value("Uranium"       , MediumType::Uranium)
        .value("Air"           , MediumType::Air)
        .value("AntaresWater"  , MediumType::AntaresWater)
        .value("StandardRock"  , MediumType::StandardRock)
        .value("FrejusRock"    , MediumType::FrejusRock)
        .value("Salt"          , MediumType::Salt)
        .value("MineralOil"    , MediumType::MineralOil)
    ;

    // --------------------------------------------------------------------- //
    // Particle
    // --------------------------------------------------------------------- //

    std::string (PROPOSALParticle::*getNameParticle)() const = &PROPOSALParticle::GetName;

    class_<PROPOSALParticle, boost::shared_ptr<PROPOSALParticle> >("Particle",
                                                                  init<ParticleType::Enum>(
                                                                  (arg("particle_type") = ParticleType::MuMinus)))

        .def(init<const PROPOSALParticle&>())

        .def(self_ns::str(self_ns::self))
        .def(self_ns::repr(self_ns::self))

        .add_property("energy", &PROPOSALParticle::GetEnergy, &PROPOSALParticle::SetEnergy)
        .add_property("propagated_distance", &PROPOSALParticle::GetPropagatedDistance, &PROPOSALParticle::SetPropagatedDistance)
        .add_property("X", &PROPOSALParticle::GetX, &PROPOSALParticle::SetX)
        .add_property("Y", &PROPOSALParticle::GetY, &PROPOSALParticle::SetY)
        .add_property("Z", &PROPOSALParticle::GetZ, &PROPOSALParticle::SetZ)
        .add_property("T", &PROPOSALParticle::GetT, &PROPOSALParticle::SetT)
        .add_property("theta", &PROPOSALParticle::GetTheta, &PROPOSALParticle::SetTheta)
        .add_property("phi", &PROPOSALParticle::GetPhi, &PROPOSALParticle::SetPhi)
        .add_property("momentum", &PROPOSALParticle::GetMomentum, &PROPOSALParticle::SetMomentum)
        .add_property("mass", &PROPOSALParticle::GetMass, &PROPOSALParticle::SetMass)
        .add_property("lifetime", &PROPOSALParticle::GetLifetime, &PROPOSALParticle::SetLifetime)
        .add_property("charge", &PROPOSALParticle::GetCharge, &PROPOSALParticle::SetCharge)
        .add_property("name", getNameParticle)
        .add_property("low", &PROPOSALParticle::GetLow, &PROPOSALParticle::SetLow)
        .add_property("type", &PROPOSALParticle::GetType, &PROPOSALParticle::SetType)
        .add_property("parent_particle_id", &PROPOSALParticle::GetParentParticleId, &PROPOSALParticle::SetParentParticleId)
        .add_property("parent_particle_energy", &PROPOSALParticle::GetParentParticleEnergy, &PROPOSALParticle::SetParentParticleEnergy)
        .add_property("particle_id", &PROPOSALParticle::GetParticleId, &PROPOSALParticle::SetParticleId)

        .add_property("Xi", &PROPOSALParticle::GetXi, &PROPOSALParticle::SetXi)
        .add_property("Yi", &PROPOSALParticle::GetYi, &PROPOSALParticle::SetYi)
        .add_property("Zi", &PROPOSALParticle::GetZi, &PROPOSALParticle::SetZi)
        .add_property("Ti", &PROPOSALParticle::GetTi, &PROPOSALParticle::SetTi)
        .add_property("Ei", &PROPOSALParticle::GetEi, &PROPOSALParticle::SetEi)

        .add_property("Xf", &PROPOSALParticle::GetXf, &PROPOSALParticle::SetXf)
        .add_property("Yf", &PROPOSALParticle::GetYf, &PROPOSALParticle::SetYf)
        .add_property("Zf", &PROPOSALParticle::GetZf, &PROPOSALParticle::SetZf)
        .add_property("Tf", &PROPOSALParticle::GetTf, &PROPOSALParticle::SetTf)
        .add_property("Ef", &PROPOSALParticle::GetEf, &PROPOSALParticle::SetEf)

        .add_property("Xc", &PROPOSALParticle::GetXc, &PROPOSALParticle::SetXc)
        .add_property("Yc", &PROPOSALParticle::GetYc, &PROPOSALParticle::SetYc)
        .add_property("Zc", &PROPOSALParticle::GetZc, &PROPOSALParticle::SetZc)
        .add_property("Tc", &PROPOSALParticle::GetTc, &PROPOSALParticle::SetTc)
        .add_property("Ec", &PROPOSALParticle::GetEc, &PROPOSALParticle::SetEc)

        .add_property("energy_lost", &PROPOSALParticle::GetElost, &PROPOSALParticle::SetElost)
    ;

    // --------------------------------------------------------------------- //
    // EnergyCutSettings
    // --------------------------------------------------------------------- //

    class_<EnergyCutSettings, boost::shared_ptr<EnergyCutSettings> >("EnergyCutSettings", init<>())

        .def(init<double, double>((arg("ecut"), arg("vcut"))))
        .def(init<const EnergyCutSettings&>())

        .def(self_ns::str(self_ns::self))
        .def(self_ns::repr(self_ns::self))

        .add_property("ecut", &EnergyCutSettings::GetEcut, &EnergyCutSettings::SetEcut)
        .add_property("vcut", &EnergyCutSettings::GetVcut, &EnergyCutSettings::SetVcut)

        .def("get_cut", &EnergyCutSettings::GetCut, "Return the lower from E*v = e")
    ;

    // --------------------------------------------------------------------- //
    // Medium
    // --------------------------------------------------------------------- //

    std::vector<std::string> (Medium::*getNameMed)() const = &Medium::GetElementName;

    class_<Medium, boost::shared_ptr<Medium> >("Medium", init<>())

        .def(init<MediumType::Enum, double>((arg("medium_type"), arg("rho") = 1.0)))
        .def(init<const Medium&>())

        .def(self_ns::str(self_ns::self))
        .def(self_ns::repr(self_ns::self))

        .add_property("num_components", &Medium::GetNumComponents, &Medium::SetNumComponents)
        .add_property("nuc_charge", &Medium::GetNucCharge, &Medium::SetNucCharge)
        .add_property("atomic_num", &Medium::GetAtomicNum, &Medium::SetAtomicNum)
        .add_property("atom_in_molecule", &Medium::GetAtomInMolecule, &Medium::SetAtomInMolecule)
        .add_property("sum_charge", &Medium::GetSumCharge, &Medium::SetSumCharge)
        .add_property("ZA", &Medium::GetZA, &Medium::SetZA)
        .add_property("I", &Medium::GetI, &Medium::SetI)
        .add_property("C", &Medium::GetC, &Medium::SetC)
        .add_property("A", &Medium::GetA, &Medium::SetA)
        .add_property("M", &Medium::GetM, &Medium::SetM)
        .add_property("X0", &Medium::GetX0, &Medium::SetX0)
        .add_property("X1", &Medium::GetX1, &Medium::SetX1)
        .add_property("D0", &Medium::GetD0, &Medium::SetD0)
        .add_property("R", &Medium::GetR, &Medium::SetR)
        .add_property("log_constant", &Medium::GetLogConstant)
        .add_property("b_prime", &Medium::GetBPrime)
        .add_property("rho", &Medium::GetRho, &Medium::SetRho)
        .add_property("mass_density", &Medium::GetMassDensity, &Medium::SetMassDensity)
        .add_property("average_nucleon_weight", &Medium::GetAverageNucleonWeight, &Medium::SetAverageNucleonWeight)
        .add_property("element_name", getNameMed, &Medium::SetElementName)
        .add_property("mol_density", &Medium::GetMolDensity, &Medium::SetMolDensity)
        .add_property("name", &Medium::GetName, &Medium::SetName)
        .add_property("type", &Medium::GetType, &Medium::SetType)
        .add_property("MN", &Medium::GetMN, &Medium::SetMN)
        .add_property("MM", &Medium::GetMM, &Medium::SetMM)
        .add_property("sum_nucleons", &Medium::GetSumNucleons, &Medium::SetSumNucleons)
        .add_property("R0", &Medium::GetR0, &Medium::SetR0)
    ;

    // --------------------------------------------------------------------- //
    // Propagator
    // --------------------------------------------------------------------- //

    class_<Propagator, boost::shared_ptr<Propagator> >("Propagator",
                                                      init<std::string, PROPOSALParticle*, bool>(
                                                      (arg("config"),
                                                       arg("particle"),
                                                       arg("applyoptions") = true)))

        .def(init<
             Medium*,
             EnergyCutSettings*,
             ParticleType::Enum,
             std::string,
             bool,
             bool,
             bool,
             bool,
             ParametrizationType::Enum,
             ParametrizationType::Enum,
             double,
             double,
             double,
             double,
             bool,
             int>(
            (arg("medium"),
             arg("energy_cuts"),
             arg("particle_type"),
             arg("path_to_tables") = "",
             arg("moliere") = true,
             arg("continuous_rand") = true,
             arg("exact_time") = true,
             arg("lpm") = true,
             arg("brems") = ParametrizationType::BremsKelnerKokoulinPetrukhin,
             arg("photo") = ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich,
             arg("brems_multiplier") = 1,
             arg("photo_multiplier") = 1,
             arg("ioniz_multiplier") = 1,
             arg("epair_multiplier") = 1,
             arg("integrate") = false,
             arg("scattering_model") = 0
            )))
        .def(init<
             ParticleType::Enum,
             std::string,
             bool,
             bool,
             bool,
             int>(
            (arg("particle_type"),
             arg("path_to_tables") = "",
             arg("exact_time") = true,
             arg("lpm") = true,
             arg("integrate") = false,
             arg("scattering_model") = 0
            )))
        .def(init<const Propagator&>())

        .def("propagate", &Propagator::propagate, (arg("max_distance_cm") = 1e20))
        .def("apply_options", &Propagator::ApplyOptions)
        .def("reset_particle", &Propagator::ResetParticle)

        // .add_property("particle", &Propagator::SetParticle)
        .add_property("particle", make_function(&Propagator::GetParticle, return_internal_reference<>()), &Propagator::SetParticle)
        .add_property("seed",&Propagator::GetSeed ,&Propagator::SetSeed)
        .add_property("brems",&Propagator::GetBrems ,&Propagator::SetBrems)
        .add_property("photo",&Propagator::GetPhoto ,&Propagator::SetPhoto)
        .add_property("path_to_tables",&Propagator::GetPath_to_tables ,&Propagator::SetPath_to_tables)
        .add_property("stopping_decay",&Propagator::GetStopping_decay ,&Propagator::SetStopping_decay)
        .add_property("current_collection",make_function(&Propagator::GetCurrentCollection, return_internal_reference<>()))
        .add_property("collections",&Propagator::GetCollections, &Propagator::SetCollections)
        .add_property("detector", make_function(&Propagator::GetDetector, return_internal_reference<>()), &Propagator::SetDetector)
    ;

    // --------------------------------------------------------------------- //
    // Cross sections
    // --------------------------------------------------------------------- //

    double (CrossSections::*CalculatedNdx)() = &CrossSections::CalculatedNdx;
    double (CrossSections::*CalculatedNdxRnd)(double) = &CrossSections::CalculatedNdx;

    class_<CrossSections, boost::shared_ptr<CrossSections>, boost::noncopyable>("CrossSections", no_init)

        .def(self_ns::str(self_ns::self))
        .def(self_ns::repr(self_ns::self))

        .add_property("cut_setting", make_function(&CrossSections::GetEnergyCutSettings, return_internal_reference<>()), &CrossSections::SetEnergyCutSettings)
        .add_property("medium", make_function(&CrossSections::GetMedium, return_internal_reference<>()), &CrossSections::SetMedium)
        .add_property("parametrization", &CrossSections::GetParametrization, &CrossSections::SetParametrization)
        .add_property("name", &CrossSections::GetName)
        .add_property("particle", make_function(&CrossSections::GetParticle, return_internal_reference<>()))

        .def("calculate_dEdx", &CrossSections::CalculatedEdx, "Calculates dE/dx")
        .def("calculate_dNdx", CalculatedNdx, "Calculates dN/dx")
        .def("calculate_dNdx", CalculatedNdxRnd, "Calculates dN/dx with random number rnd")
        .def("enable_dEdx_interpolation", &CrossSections::EnableDEdxInterpolation, (arg("path") = "", arg("raw") = false))
        .def("enable_dNdx_interpolation", &CrossSections::EnableDNdxInterpolation, (arg("path") = "", arg("raw") = false))
        .def("disable_dEdx_interpolation", &CrossSections::DisableDEdxInterpolation)
        .def("disable_dNdx_interpolation", &CrossSections::DisableDNdxInterpolation)
        .def("enable_lpm_effect", &CrossSections::EnableLpmEffect)
    ;

    class_<Photonuclear, boost::shared_ptr<Photonuclear>, bases<CrossSections> >("Photonuclear", init<PROPOSALParticle*, Medium*, EnergyCutSettings*>())
        .def(self_ns::str(self_ns::self))
        .def(self_ns::repr(self_ns::self))
    ;
    class_<Epairproduction, boost::shared_ptr<Epairproduction>, bases<CrossSections> >("Epairproduction", init<PROPOSALParticle*, Medium*, EnergyCutSettings*>())
        .def(self_ns::str(self_ns::self))
        .def(self_ns::repr(self_ns::self))
    ;
    class_<Bremsstrahlung, boost::shared_ptr<Bremsstrahlung>, bases<CrossSections> >("Bremsstrahlung", init<PROPOSALParticle*, Medium*, EnergyCutSettings*>())
        .def(self_ns::str(self_ns::self))
        .def(self_ns::repr(self_ns::self))
    ;
    class_<Ionization, boost::shared_ptr<Ionization>, bases<CrossSections> >("Ionization", init<PROPOSALParticle*, Medium*, EnergyCutSettings*>())
        .def(self_ns::str(self_ns::self))
        .def(self_ns::repr(self_ns::self))
    ;

    // --------------------------------------------------------------------- //
    // Decay
    // --------------------------------------------------------------------- //

    class_<Decay, boost::shared_ptr<Decay> >("Decay", init<>())

        .def(init<const Decay&>())
        .def(init<PROPOSALParticle*>((arg("particle"))))

        // .def(self_ns::str(self_ns::self))
        // .def(self_ns::repr(self_ns::self))

        .add_property("out", &Decay::GetOut, &Decay::SetOut)
        .add_property("particle", make_function(&Decay::GetParticle, return_internal_reference<>()), &Decay::SetParticle)
        .add_property("store_neutrinos", &Decay::GetStoreNeutrinos, &Decay::SetStoreNeutrinos)
        .add_property("multiplier", &Decay::GetMultiplier, &Decay::SetMultiplier)

        .def("make_decay", &Decay::MakeDecay, "Cross section describing the decay.")
        .def("calculate_product_energy", &Decay::CalculateProductEnergy, "Calculates the product energy.")
    ;

    // ------------------------------------------------------------------------- //
    // ProcessCollection
    // ------------------------------------------------------------------------- //

    class_<ProcessCollection, boost::shared_ptr<ProcessCollection> >("ProcessCollection", init<>())

        .def(init<const ProcessCollection&>())
        .def(init<PROPOSALParticle*, Medium*, EnergyCutSettings*>((arg("particle"), arg("medium"), arg("energy_cut"))))

        .def(self_ns::str(self_ns::self))
        .def(self_ns::repr(self_ns::self))

        .add_property("cross_sections", &ProcessCollection::GetCrosssections, &ProcessCollection::SetCrosssections)
        .add_property("cut_setting", make_function(&ProcessCollection::GetCutSettings, return_internal_reference<>()), &ProcessCollection::SetCutSettings)
        .add_property("lpm_effect", &ProcessCollection::GetLpmEffectEnabled, &ProcessCollection::SetLpmEffectEnabled)
        .add_property("medium", make_function(&ProcessCollection::GetMedium, return_internal_reference<>()), &ProcessCollection::SetMedium)
        .add_property("geometry", make_function(&ProcessCollection::GetGeometry, return_internal_reference<>()), &ProcessCollection::SetGeometry)
        .add_property("particle", make_function(&ProcessCollection::GetParticle, return_internal_reference<>()), &ProcessCollection::SetParticle)
        .add_property("location", &ProcessCollection::GetLocation, &ProcessCollection::SetLocation)
        .add_property("density_correction", &ProcessCollection::GetDensityCorrection, &ProcessCollection::SetDensityCorrection)
        .add_property("enable_randomization", &ProcessCollection::GetEnableRandomization, &ProcessCollection::SetEnableRandomization)
    ;

    // ------------------------------------------------------------------------- //
    // Geometry
    // ------------------------------------------------------------------------- //

    class_<Geometry, boost::shared_ptr<Geometry> >("Geometry", init<>())

        .def(init<const Geometry&>())

        .def(self_ns::str(self_ns::self))
        .def(self_ns::repr(self_ns::self))

        .def("init_box", &Geometry::InitBox, return_internal_reference<>(), "All units in [m]", (
            arg("x0"),
            arg("y0"),
            arg("z0"),
            arg("x"),
            arg("y"),
            arg("z")
        ))
        .def("init_sphere", &Geometry::InitSphere, return_internal_reference<>(), "All units in [m]", (
            arg("x0"),
            arg("y0"),
            arg("z0"),
            arg("radius"),
            arg("inner_radius")
        ))
        .def("init_cylinder", &Geometry::InitCylinder, return_internal_reference<>(), "All units in [m]", (
            arg("x0"),
            arg("y0"),
            arg("z0"),
            arg("radius"),
            arg("inner_radius"),
            arg("z")
        ))
        .def("is_particle_inside", &Geometry::IsParticleInside)
        .def("is_particle_infront", &Geometry::IsParticleInfront)
        .def("is_particle_behind", &Geometry::IsParticleBehind)
        .def("distance_to_border", &Geometry::DistanceToBorder)

        .add_property("x", &Geometry::GetX, &Geometry::SetX)
        .add_property("y", &Geometry::GetY, &Geometry::SetY)
        .add_property("z", &Geometry::GetZ, &Geometry::SetZ)
        .add_property("x0", &Geometry::GetX0, &Geometry::SetX0)
        .add_property("y0", &Geometry::GetY0, &Geometry::SetY0)
        .add_property("z0", &Geometry::GetZ0, &Geometry::SetZ0)
        .add_property("inner_radius", &Geometry::GetInnerRadius, &Geometry::SetInnerRadius)
        .add_property("radius", &Geometry::GetRadius, &Geometry::SetRadius)
        .add_property("object", &Geometry::GetObject, &Geometry::SetObject)
        .add_property("hirachy", &Geometry::GetHirarchy, &Geometry::SetHirarchy)
    ;
}
