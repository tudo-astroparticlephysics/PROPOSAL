#define BOOST_PYTHON_MAX_ARITY 17

// #include <string>

#include <boost/python.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <map>

// #include "PROPOSAL/Propagator.h"
#include "PROPOSAL/PROPOSAL.h"

#define PARTICLE_DEF(cls)                                                                                              \
    class_<cls##Def, boost::shared_ptr<cls##Def>, bases<ParticleDef>, boost::noncopyable>(#cls "Def", no_init)         \
                                                                                                                       \
        .def("get", make_function(&cls##Def::Get, return_value_policy<reference_existing_object>()))                   \
        .staticmethod("get");

#define COMPONENT_DEF(cls)\
    class_<Components::cls, boost::shared_ptr<Components::cls>, bases<Components::Component> >( #cls, init<double>());

#define MEDIUM_DEF(cls)\
    class_<cls, boost::shared_ptr<cls>, bases<Medium> >( #cls, init<double>());

#define BREMS_DEF(cls)                                                                                                 \
    class_<Brems##cls, boost::shared_ptr<Brems##cls>, bases<Bremsstrahlung> >(                                         \
        #cls,                                                                                                          \
        init<const ParticleDef&, const Medium&, const EnergyCutSettings&, double, bool>(                               \
            (arg("particle_def"), arg("medium"), arg("energy_cuts"), arg("multiplier"), arg("lpm"))));

#define PHOTO_REAL_DEF(cls, parent)                                                                                    \
    class_<Photo##cls, boost::shared_ptr<Photo##cls>, bases<Photo##parent> >(                                          \
        #cls,                                                                                                          \
        init<const ParticleDef&, const Medium&, const EnergyCutSettings&, const RealPhoton&, double>(                  \
            (arg("particle_def"), arg("medium"), arg("energy_cuts"), arg("real_photon"), arg("multiplier"))));

#define PHOTO_Q2_DEF(cls)                                                                                              \
    class_<Photo##cls, boost::shared_ptr<Photo##cls>, bases<PhotoQ2Integral> >(                                        \
        #cls,                                                                                                          \
        init<const ParticleDef&, const Medium&, const EnergyCutSettings&, const ShadowEffect&, double>(                \
            (arg("particle_def"), arg("medium"), arg("energy_cuts"), arg("shadow_effect"), arg("multiplier"))));

#define PHOTO_Q2_INTERPOL_DEF(cls)                                                                                     \
    class_<PhotoQ2Interpolant<Photo##cls>, boost::shared_ptr<PhotoQ2Interpolant<Photo##cls> >, bases<Photo##cls> >(    \
        #cls "Interpolant",                                                                                            \
        init<const ParticleDef&,                                                                                       \
             const Medium&,                                                                                            \
             const EnergyCutSettings&,                                                                                 \
             const ShadowEffect&,                                                                                      \
             double,                                                                                                   \
             InterpolationDef>((arg("particle_def"),                                                                   \
                                arg("medium"),                                                                         \
                                arg("energy_cuts"),                                                                    \
                                arg("shadow_effect"),                                                                  \
                                arg("multiplier"),                                                                     \
                                arg("interpolation_def"))));

using namespace PROPOSAL;

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

// using namespace PROPOSAL;
//
// // struct CrossSectionsWrap : CrossSections, boost::python::wrapper<CrossSections>
// // {
// //     CrossSectionsWrap(): CrossSections(), wrapper<CrossSections>(){}
// //
// //     // virtual double CalculatedEdx()
// //     // {
// //     //     return this->get_override("calculate_dEdx")();
// //     // }
// // };
//
//
// #<{(|*****************************************************************************
// *                             Register functions                              *
// *****************************************************************************|)}>#
//
// // ------------------------------------------------------------------------- //
// // Vectors
// // ------------------------------------------------------------------------- //
//
//
// struct CrossSectionToPython
// {
//     static PyObject* convert(std::vector<CrossSections*> const& vec)
//     {
//         boost::python::list python_list;
//
//         for (unsigned int i = 0; i < vec.size(); ++i)
//         {
//             switch (vec[i]->GetType())
//             {
//                 case ParticleType::Brems:
//                     python_list.append(boost::python::ptr((Bremsstrahlung*)vec[i]));
//                     break;
//                 case ParticleType::NuclInt:
//                     python_list.append(boost::python::ptr((Photonuclear*)vec[i]));
//                     break;
//                 case ParticleType::DeltaE:
//                     python_list.append(boost::python::ptr((Ionization*)vec[i]));
//                     break;
//                 case ParticleType::EPair:
//                     python_list.append(boost::python::ptr((Epairproduction*)vec[i]));
//                     break;
//                 default:
//                     boost::python::object obj(NULL);
//                     break;
//             }
//         }
//
//         PyObject* py = boost::python::incref(python_list.ptr());
//         return py;
//     }
// };
//

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

// // ------------------------------------------------------------------------- //
// // Pair
// // ------------------------------------------------------------------------- //
//
//
// template<typename T1, typename T2>
// struct PairToPythonList
// {
//     static PyObject* convert(std::pair<T1, T2> const& p)
//     {
//         boost::python::list python_list;
//
//         python_list.append(boost::python::object(p.first));
//         python_list.append(boost::python::object(p.second));
//         // typename std::pair<T1, T2>::const_iterator iter;
//
//         // for(iter = vec.begin(); iter != vec.end(); ++iter)
//         // {
//         //     python_list.append(boost::python::object(*iter));
//         // }
//
//         return boost::python::incref(python_list.ptr());
//     }
// };
//
//
// template<typename T1, typename T2>
// struct PairFromPythonList
// {
//
//     PairFromPythonList()
//     {
//         boost::python::converter::registry::push_back(&PairFromPythonList<T1, T2>::convertible,
//                             &PairFromPythonList<T1, T2>::construct,
//                             boost::python::type_id<std::pair<T1, T2> >());
//     }
//
//     // Determine if obj_ptr can be converted in a std::pair
//     static void* convertible(PyObject* obj_ptr)
//     {
//         if (!PyList_Check(obj_ptr))
//         {
//             return 0;
//         }
//
//         return obj_ptr;
//     }
//
//     // Convert obj_ptr into a std::pair<T>
//     static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
//     {
//         // Use borrowed to construct the object so that a reference
//         // count will be properly handled.
//         boost::python::list python_list(boost::python::handle<>(boost::python::borrowed(obj_ptr)));
//
//         // Grab pointer to memory into which to construct the new std::pair<T1, T2>
//         void* storage = reinterpret_cast<boost::python::converter::rvalue_from_python_storage<std::pair<T1, T2> >*>(data)->storage.bytes;
//
//         // in-place construct the new std::pair<T1, T2> using the character data
//         // extraced from the python object
//         std::pair<T1, T2>& p = *(new (storage) std::pair<T1, T2>());
//
//         assert(boost::python::len(python_list == 2));
//
//         // Populate pari from python list
//         p.first = boost::python::extract<T1>(python_list[0]);
//         p.second = boost::python::extract<T2>(python_list[1]);
//
//         // Stash the memory chunk pointer for later use by boost.python
//         data->convertible = storage;
//     }
// };
//
// #<{(|*****************************************************************************
// *                               Python Module                                 *
// *****************************************************************************|)}>#

// ------------------------------------------------------------------------- //
// Components sub module
// ------------------------------------------------------------------------- //

void export_components()
{
    using namespace boost::python;
    // map the Util namespace to a sub-module
    // make "from mypackage.Util import <whatever>" work
    object componentsModule(handle<>(borrowed(PyImport_AddModule("pyPROPOSAL.Components"))));
    // make "from mypackage import Util" work
    scope().attr("Components") = componentsModule;
    // set the current scope to the new sub-module
    scope components_scope = componentsModule;
    // export stuff in the Util namespace

    class_<Components::Component, boost::shared_ptr<Components::Component>, boost::noncopyable >("Component", no_init)

        .def(self_ns::str(self_ns::self))

        .add_property("name", &Components::Component::GetName)
        .add_property("nuclear_charge", &Components::Component::GetNucCharge)
        .add_property("atomic_number", &Components::Component::GetAtomicNum)
        .add_property("atoms_in_molecule", &Components::Component::GetAtomInMolecule)
        .add_property("log_constant", &Components::Component::GetLogConstant)
        .add_property("bprime", &Components::Component::GetBPrime)
        .add_property("average_nucleon_weight", &Components::Component::GetAverageNucleonWeight)
        .add_property("mn", &Components::Component::GetMN)
        .add_property("r0", &Components::Component::GetR0);

    // class_<Components::Hydrogen, boost::shared_ptr<Components::Hydrogen>, bases<Components::Component> >( "Hydrogen", init<double>());
    COMPONENT_DEF(Hydrogen)
    COMPONENT_DEF(Carbon)
    COMPONENT_DEF(Nitrogen)
    COMPONENT_DEF(Oxygen)
    COMPONENT_DEF(Sodium)
    COMPONENT_DEF(Magnesium)
    COMPONENT_DEF(Sulfur)
    COMPONENT_DEF(Argon)
    COMPONENT_DEF(Potassium)
    COMPONENT_DEF(Calcium)
    COMPONENT_DEF(Iron)
    COMPONENT_DEF(Copper)
    COMPONENT_DEF(Lead)
    COMPONENT_DEF(Uranium)
    COMPONENT_DEF(StandardRock)
    COMPONENT_DEF(FrejusRock)
}


void export_medium()
{
    using namespace boost::python;
    // map the Util namespace to a sub-module
    // make "from mypackage.Util import <whatever>" work
    object mediumModule(handle<>(borrowed(PyImport_AddModule("pyPROPOSAL.Medium"))));
    // make "from mypackage import Util" work
    scope().attr("Medium") = mediumModule;
    // set the current scope to the new sub-module
    scope medium_scope = mediumModule;
    // export stuff in the Util namespace

    class_<Medium, boost::shared_ptr<Medium>, boost::noncopyable >("Medium", no_init)

        .def(self_ns::str(self_ns::self))

        .add_property("name", &Medium::GetName);

    MEDIUM_DEF(Water)
    MEDIUM_DEF(Ice)
    MEDIUM_DEF(Salt)
    // MEDIUM_DEF(CalciumCarbonate) //TODO(mario): Whats going on here??? symbol not found at runtime Sun 2017/09/24
    MEDIUM_DEF(StandardRock)
    MEDIUM_DEF(FrejusRock)
    MEDIUM_DEF(Iron)
    MEDIUM_DEF(Hydrogen)
    MEDIUM_DEF(Lead)
    MEDIUM_DEF(Copper)
    MEDIUM_DEF(Uranium)
    MEDIUM_DEF(Air)
    MEDIUM_DEF(Paraffin)
    MEDIUM_DEF(AntaresWater)
}

void export_bremsstrahlung()
{
    using namespace boost::python;
    // map the Util namespace to a sub-module
    // make "from mypackage.Util import <whatever>" work
    object bremsstrahlungModule(handle<>(borrowed(PyImport_AddModule("pyPROPOSAL.Bremsstrahlung"))));
    // make "from mypackage import Util" work
    scope().attr("Bremsstrahlung") = bremsstrahlungModule;
    // set the current scope to the new sub-module
    scope bremsstrahlung_scope = bremsstrahlungModule;
    // export stuff in the Util namespace

    /**************************************************************************
    *                       Bind all Parametrizations                         *
    **************************************************************************/

    class_<Bremsstrahlung, boost::shared_ptr<Bremsstrahlung>, bases<Parametrization>, boost::noncopyable>("Bremsstrahlung", no_init);

    BREMS_DEF(KelnerKokoulinPetrukhin)
    BREMS_DEF(PetrukhinShestakov)
    BREMS_DEF(CompleteScreening)
    BREMS_DEF(AndreevBezrukovBugaev)
}

void export_photo()
{
    using namespace boost::python;
    // map the Util namespace to a sub-module
    // make "from mypackage.Util import <whatever>" work
    object photoModule(handle<>(borrowed(PyImport_AddModule("pyPROPOSAL.Photonuclear"))));
    // make "from mypackage import Util" work
    scope().attr("Photonuclear") = photoModule;
    // set the current scope to the new sub-module
    scope photo_scope = photoModule;
    // export stuff in the Util namespace

    /**************************************************************************
    *                       Bind all Parametrizations                         *
    **************************************************************************/

    // Shadow Effect
    class_<ShadowEffect, boost::shared_ptr<ShadowEffect>, bases<Parametrization>, boost::noncopyable>("ShadowEffect", no_init)
        .def("calculate_shadow_effect", &ShadowEffect::CalculateShadowEffect);

    class_<ShadowDutta, boost::shared_ptr<ShadowDutta>, bases<ShadowEffect> >("ShadowDutta", init<>());
    class_<ShadowButkevichMikhailov, boost::shared_ptr<ShadowButkevichMikhailov>, bases<ShadowEffect> >("ShadowButkevichMikhailov", init<>());

    // Real Photon
    class_<RealPhoton, boost::shared_ptr<RealPhoton>, boost::noncopyable>("RealPhoton", no_init)
        .def("calculate_hardbb", &RealPhoton::CalculateHardBB);

    class_<SoftBB, boost::shared_ptr<SoftBB>, bases<RealPhoton> >("SoftBB", init<>());
    class_<HardBB, boost::shared_ptr<HardBB>, bases<RealPhoton> >("HardBB", init<const ParticleDef&>((arg("particle_def"))));

    // Photnuclear
    class_<Photonuclear, boost::shared_ptr<Photonuclear>, bases<Parametrization>, boost::noncopyable>("Photonuclear", no_init);
    class_<PhotoRealPhotonAssumption, boost::shared_ptr<PhotoRealPhotonAssumption>, bases<Photonuclear>, boost::noncopyable>("PhotoRealPhotonAssumption", no_init);
    class_<PhotoQ2Integral, boost::shared_ptr<PhotoQ2Integral>, bases<Photonuclear>, boost::noncopyable>("PhotoQ2Integral", no_init);

    PHOTO_REAL_DEF(Zeus, RealPhotonAssumption)
    PHOTO_REAL_DEF(BezrukovBugaev, RealPhotonAssumption)
    PHOTO_REAL_DEF(Rhode, RealPhotonAssumption)
    PHOTO_REAL_DEF(Kokoulin, BezrukovBugaev) // Kokoulin derives from BezrukovBugaev

    PHOTO_Q2_DEF(AbramowiczLevinLevyMaor91)
    PHOTO_Q2_DEF(AbramowiczLevinLevyMaor97)
    PHOTO_Q2_DEF(ButkevichMikhailov)
    PHOTO_Q2_DEF(RenoSarcevicSu)

    PHOTO_Q2_INTERPOL_DEF(AbramowiczLevinLevyMaor91)
    PHOTO_Q2_INTERPOL_DEF(AbramowiczLevinLevyMaor97)
    PHOTO_Q2_INTERPOL_DEF(ButkevichMikhailov)
    PHOTO_Q2_INTERPOL_DEF(RenoSarcevicSu)
}

void export_epair()
{
    using namespace boost::python;
    // map the Util namespace to a sub-module
    // make "from mypackage.Util import <whatever>" work
    object epairModule(handle<>(borrowed(PyImport_AddModule("pyPROPOSAL.EpairProduction"))));
    // make "from mypackage import Util" work
    scope().attr("EpairProduction") = epairModule;
    // set the current scope to the new sub-module
    scope epair_scope = epairModule;
    // export stuff in the Util namespace

    /**************************************************************************
    *                       Bind all Parametrizations                         *
    **************************************************************************/

    class_<EpairProduction, boost::shared_ptr<EpairProduction>, bases<Parametrization>, boost::noncopyable>("EpairProduction", no_init);

    class_<EpairProductionRhoIntegral, boost::shared_ptr<EpairProductionRhoIntegral>, bases<EpairProduction> >(
        "EpairProductionRhoIntegral",
        init<const ParticleDef&, const Medium&, const EnergyCutSettings&, double, bool>());
            // (arg("particle_def"), arg("medium"), arg("energy_cuts"), arg("multiplier"), arg("lpm"))));

    class_<EpairProductionRhoInterpolant, boost::shared_ptr<EpairProductionRhoInterpolant>, bases<EpairProductionRhoIntegral> >(
        "EpairProductionRhoInterpolant",
        init<const ParticleDef&, const Medium&, const EnergyCutSettings&, double, bool, InterpolationDef>(
            (arg("particle_def"), arg("medium"), arg("energy_cuts"), arg("multiplier"), arg("lpm"), arg("interpolation_def"))));

}

void export_parametrizations()
{
    using namespace boost::python;
    // map the Util namespace to a sub-module
    // make "from mypackage.Util import <whatever>" work
    object parametrizationModule(handle<>(borrowed(PyImport_AddModule("pyPROPOSAL.Parametrization"))));
    // make "from mypackage import Util" work
    scope().attr("Parametrization") = parametrizationModule;
    // set the current scope to the new sub-module
    scope parametrization_scope = parametrizationModule;
    // export stuff in the Util namespace

    /**************************************************************************
    *                       Bind all Parametrizations                         *
    **************************************************************************/

    class_<Parametrization::IntegralLimits, boost::shared_ptr<Parametrization::IntegralLimits> >("IntegralLimits")

        .def_readwrite("v_max", &Parametrization::IntegralLimits::vMax)
        .def_readwrite("v_up", &Parametrization::IntegralLimits::vUp)
        .def_readwrite("v_min", &Parametrization::IntegralLimits::vMin);

    class_<Parametrization, boost::shared_ptr<Parametrization>, boost::noncopyable>("Parametrization", no_init)

        .def(self_ns::str(self_ns::self))

        .def("differntial_crosssection", &Parametrization::DifferentialCrossSection)
        .def("dEdx_integrand", &Parametrization::FunctionToDEdxIntegral)
        .def("dE2dx_integrand", &Parametrization::FunctionToDE2dxIntegral)
        .def("dNdx_integrand", &Parametrization::FunctionToDNdxIntegral)

        .add_property("name", make_function(&Parametrization::GetName, return_value_policy<copy_const_reference>()))
        .add_property("integral_limits", &Parametrization::GetIntegralLimits)
        .add_property("particle_def", make_function(&Parametrization::GetParticleDef, return_internal_reference<>()))
        .add_property("medium", make_function(&Parametrization::GetMedium, return_internal_reference<>()))
        .add_property("energy_cuts", make_function(&Parametrization::GetEnergyCuts, return_internal_reference<>()))
        // .add_property("medium", &Parametrization::GetMedium)
        // .add_property("energy_cuts", &Parametrization::GetEnergyCuts)
        .add_property("multiplier", &Parametrization::GetMultiplier)
        .add_property("hash", &Parametrization::GetHash);

    class_<Ionization, boost::shared_ptr<Ionization>, bases<Parametrization> >(
        "Ionization",
        init<const ParticleDef&, const Medium&, const EnergyCutSettings&, double>());

    export_epair();
    export_bremsstrahlung();
    export_photo();
}

void export_crosssections()
{
    using namespace boost::python;
    // map the Util namespace to a sub-module
    // make "from mypackage.Util import <whatever>" work
    object crosssectionModule(handle<>(borrowed(PyImport_AddModule("pyPROPOSAL.CrossSection"))));
    // make "from mypackage import Util" work
    scope().attr("Crosssection") = crosssectionModule;
    // set the current scope to the new sub-module
    scope crosssection_scope = crosssectionModule;
    // export stuff in the Util namespace


    /**************************************************************************
    *                        Bind all cross sections                          *
    **************************************************************************/

    double (CrossSection::*dndx)(double) = &CrossSection::CalculatedNdx;
    double (CrossSection::*stochastic)(double, double, double) = &CrossSection::CalculateStochasticLoss;

    class_<CrossSection, boost::shared_ptr<CrossSection>, boost::noncopyable>("CrossSection", no_init)

        .def(self_ns::str(self_ns::self))

        .def("calculate_dEdx", &CrossSection::CalculatedEdx)
        .def("calculate_dE2dx", &CrossSection::CalculatedE2dx)
        .def("calculate_dNdx", dndx)
        .def("calculate_stochastic_loss", stochastic)

        .add_property("id", &CrossSection::GetTypeId)
        .add_property("parametrization", make_function(&CrossSection::GetParametrization, return_internal_reference<>()));

    class_<CrossSectionIntegral, boost::shared_ptr<CrossSectionIntegral>, bases<CrossSection>, boost::noncopyable>("CrossSectionIntegral", no_init);
    class_<CrossSectionInterpolant, boost::shared_ptr<CrossSectionInterpolant>, bases<CrossSection>, boost::noncopyable>("CrossSectionInterpolant", no_init);

    class_<BremsIntegral, boost::shared_ptr<BremsIntegral>, bases<CrossSectionIntegral> >("BremsIntegral", init<const Bremsstrahlung&>((arg("parametrization"))));
    class_<EpairIntegral, boost::shared_ptr<EpairIntegral>, bases<CrossSectionIntegral> >("EpairIntegral", init<const EpairProduction&>((arg("parametrization"))));
    class_<PhotoIntegral, boost::shared_ptr<PhotoIntegral>, bases<CrossSectionIntegral> >("PhotoIntegral", init<const Photonuclear&>((arg("parametrization"))));
    class_<IonizIntegral, boost::shared_ptr<IonizIntegral>, bases<CrossSectionIntegral> >("IonizIntegral", init<const Ionization&>((arg("parametrization"))));

    class_<BremsInterpolant, boost::shared_ptr<BremsInterpolant>, bases<CrossSectionInterpolant> >("BremsInterpolant", init<const Bremsstrahlung&, InterpolationDef>((arg("parametrization"), arg("interpolation_def"))));
    class_<EpairInterpolant, boost::shared_ptr<EpairInterpolant>, bases<CrossSectionInterpolant> >("EpairInterpolant", init<const EpairProduction&, InterpolationDef>((arg("parametrization"), arg("interpolation_def"))));
    class_<PhotoInterpolant, boost::shared_ptr<PhotoInterpolant>, bases<CrossSectionInterpolant> >("PhotoInterpolant", init<const Photonuclear&, InterpolationDef>((arg("parametrization"), arg("interpolation_def"))));
    class_<IonizInterpolant, boost::shared_ptr<IonizInterpolant>, bases<CrossSectionInterpolant> >("IonizInterpolant", init<const Ionization&, InterpolationDef>((arg("parametrization"), arg("interpolation_def"))));
}

class PythonHardBBTables {};

BOOST_PYTHON_MODULE(pyPROPOSAL)
{
    using namespace boost::python;
    object package = scope();
    package.attr("__path__") = "pyPROPOSAL";

    export_components();
    export_medium();
    export_parametrizations();
    export_crosssections();

    docstring_options doc_options(true, true, false);

    // --------------------------------------------------------------------- //
    // Vector classes
    // --------------------------------------------------------------------- //

    // ----[ Register the to-python converter ]-------------- //

    to_python_converter< std::vector<double>, VectorToPythonList<double> >();
    to_python_converter< std::vector<std::vector<double> >, VectorToPythonList<std::vector<double> > > ();

    to_python_converter< std::vector<std::string>, VectorToPythonList<std::string> >();

    to_python_converter< std::vector<DynamicData*>, PVectorToPythonList<DynamicData*> >();
    to_python_converter< std::vector<PROPOSALParticle*>, PVectorToPythonList<PROPOSALParticle*> >();

    to_python_converter< std::vector<SectorFactory::Definition>, VectorToPythonList<SectorFactory::Definition> >();

    // ----[ Register the from-python converter ]------------ //

    VectorFromPythonList<double>();
    VectorFromPythonList<std::vector<double> >();
    VectorFromPythonList<std::string>();

    VectorFromPythonList<DynamicData*>();
    VectorFromPythonList<PROPOSALParticle*>();

    VectorFromPythonList<SectorFactory::Definition>();

    /**************************************************************************
    *                                Vector3D                                *
    **************************************************************************/

    class_<Vector3D, boost::shared_ptr<Vector3D> >("Vector3D", init<>())

        .def(init<double, double, double>((arg("x"), arg("y"), arg("z"))))
        .def(init<const Vector3D&>())

        .def(self_ns::str(self_ns::self))

        .add_property("x", &Vector3D::GetX)
        .add_property("y", &Vector3D::GetY)
        .add_property("z", &Vector3D::GetZ)
        .add_property("radius", &Vector3D::GetRadius)
        .add_property("phi", &Vector3D::GetPhi)
        .add_property("theta", &Vector3D::GetTheta)

        .def("set_cartesian_coordinates", &Vector3D::SetCartesianCoordinates)
        .def("set_spherical_coordinates", &Vector3D::SetSphericalCoordinates)
        .def("normalize", &Vector3D::normalise)
        .def("magnitude", &Vector3D::magnitude)
        .def("cartesian_from_spherical", &Vector3D::CalculateCartesianFromSpherical)
        .def("spherical_from_cartesian", &Vector3D::CalculateSphericalCoordinates)
    ;

    /**************************************************************************
    *                               EnergyCuts                                *
    **************************************************************************/

    class_<EnergyCutSettings, boost::shared_ptr<EnergyCutSettings> >("EnergyCutSettings", init<>())

        .def(init<double, double>((arg("ecut"), arg("vcut"))))
        .def(init<const EnergyCutSettings&>())

        .def(self_ns::str(self_ns::self))

        .add_property("ecut", &EnergyCutSettings::GetEcut, &EnergyCutSettings::SetEcut)
        .add_property("vcut", &EnergyCutSettings::GetVcut, &EnergyCutSettings::SetVcut)

        .def("get_cut", &EnergyCutSettings::GetCut, "Return the lower from E*v = e")
    ;

    /**************************************************************************
    *                              DecayChannel                               *
    **************************************************************************/

    class_<DecayChannel, boost::shared_ptr<DecayChannel>, boost::noncopyable >("DecayChannel", no_init)

        .def("decay", &DecayChannel::Decay, "Decay the given paritcle")
    ;

    class_<LeptonicDecayChannel, boost::shared_ptr<LeptonicDecayChannel>, bases<DecayChannel> >("LeptonicDecayChannel", init<>());

    class_<TwoBodyPhaseSpace, boost::shared_ptr<TwoBodyPhaseSpace>, bases<DecayChannel> >("TwoBodyPhaseSpace", init<double, double>());

    class_<StableChannel, boost::shared_ptr<StableChannel>, bases<DecayChannel> >("StableChannel", init<>());

    /**************************************************************************
    *                               DecayTable                                *
    **************************************************************************/

    enum_<DecayTable::Mode>("DecayMode")
        .value("LeptonicDecay", DecayTable::LeptonicDecay)
        .value("TwoBodyDecay", DecayTable::TwoBodyDecay)
        .value("Stable", DecayTable::Stable);

    class_<DecayTable, boost::shared_ptr<DecayTable> >("DecayTable", init<>())

        .def(init<const DecayTable&>())

        // .def(self_ns::str(self_ns::self))

        .def("add_channel", &DecayTable::addChannel, "Add an decay channel")
        .def("select_channel", make_function(&DecayTable::SelectChannel, return_internal_reference<>()), "Select an decay channel according to given branching ratios")
        .def("set_stable", &DecayTable::SetStable, "Define decay table for stable particles")
    ;

    /**************************************************************************
    *                              ParticleDef                                *
    **************************************************************************/

    class_<ParticleDef, boost::shared_ptr<ParticleDef> >("ParticleDef", init<>())

        .def(init<const ParticleDef&>())

        .def(self_ns::str(self_ns::self))

        .def_readwrite("name", &ParticleDef::name)
        .def_readwrite("mass", &ParticleDef::mass)
        .def_readwrite("low", &ParticleDef::low)
        .def_readwrite("charge", &ParticleDef::charge)
        .def_readwrite("decay_table", &ParticleDef::decay_table)
        // .def_readwrite("hardbb_table", &ParticleDef::hardbb_table) //TODO(mario): hardbb Fri 2017/09/22
    ;

    class_<MuMinusDef, boost::shared_ptr<MuMinusDef>, bases<ParticleDef>, boost::noncopyable >("MuMinusDef", no_init)

        .def("get", make_function(&MuMinusDef::Get, return_value_policy<reference_existing_object>()))
        .staticmethod("get")
        ;

    PARTICLE_DEF(MuPlus)

    PARTICLE_DEF(EMinus)
    PARTICLE_DEF(EPlus)

    PARTICLE_DEF(TauMinus)
    PARTICLE_DEF(TauPlus)

    PARTICLE_DEF(StauMinus)
    PARTICLE_DEF(StauPlus)

    PARTICLE_DEF(P0)
    PARTICLE_DEF(PiMinus)
    PARTICLE_DEF(PiPlus)

    PARTICLE_DEF(KMinus)
    PARTICLE_DEF(KPlus)

    PARTICLE_DEF(PMinus)
    PARTICLE_DEF(PPlus)

    PARTICLE_DEF(NuE)
    PARTICLE_DEF(NuEBar)

    PARTICLE_DEF(NuMu)
    PARTICLE_DEF(NuMuBar)

    PARTICLE_DEF(NuTau)
    PARTICLE_DEF(NuTauBar)

    PARTICLE_DEF(Monopole)
    PARTICLE_DEF(Gamma)
    PARTICLE_DEF(StableMassiveParticle)

    /**************************************************************************
    *                              Dynamic Data                               *
    **************************************************************************/


    enum_<DynamicData::Type>("Data")
        .value("None", DynamicData::None)
        .value("Particle", DynamicData::Particle)
        .value("Brems", DynamicData::Brems)
        .value("DeltaE", DynamicData::DeltaE)
        .value("Epair", DynamicData::Epair)
        .value("NuclInt", DynamicData::NuclInt)
        .value("MuPair", DynamicData::MuPair)
        .value("Hadrons", DynamicData::Hadrons)
        .value("ContinuousEnergyLoss", DynamicData::ContinuousEnergyLoss);

    class_<DynamicData, boost::shared_ptr<DynamicData> >("DynamicData", init<DynamicData::Type>())

        .def(init<const DynamicData&>())

        .def(self_ns::str(self_ns::self))

        .add_property("id", &DynamicData::GetTypeId)
        .add_property("position", &DynamicData::GetPosition, &DynamicData::SetPosition)
        .add_property("direction", &DynamicData::GetDirection, &DynamicData::SetDirection)
        .add_property("energy", &DynamicData::GetEnergy, &DynamicData::SetEnergy)
        .add_property("parent_particle_energy", &DynamicData::GetParentParticleEnergy, &DynamicData::SetParentParticleEnergy)
        .add_property("time", &DynamicData::GetTime, &DynamicData::SetTime)
        .add_property("propagated_distance", &DynamicData::GetPropagatedDistance, &DynamicData::SetPropagatedDistance)
    ;

    /**************************************************************************
    *                                Particle                                 *
    **************************************************************************/

    class_<PROPOSALParticle, boost::shared_ptr<PROPOSALParticle>, bases<DynamicData> >("Particle", init<>())

        .def(init<const ParticleDef&>())
        .def(init<const PROPOSALParticle&>())

        .add_property("momentum", &PROPOSALParticle::GetMomentum, &PROPOSALParticle::SetMomentum)
        .add_property("particle_def", make_function(&PROPOSALParticle::GetParticleDef, return_value_policy<reference_existing_object>()))

        .add_property("entry_point", &PROPOSALParticle::GetEntryPoint, &PROPOSALParticle::SetEntryPoint)
        .add_property("entry_time", &PROPOSALParticle::GetEntryTime, &PROPOSALParticle::SetEntryTime)
        .add_property("entry_energy", &PROPOSALParticle::GetEntryEnergy, &PROPOSALParticle::SetEntryEnergy)
        .add_property("exit_point", &PROPOSALParticle::GetExitPoint, &PROPOSALParticle::SetExitPoint)
        .add_property("exit_time", &PROPOSALParticle::GetExitTime, &PROPOSALParticle::SetExitTime)
        .add_property("exit_energy", &PROPOSALParticle::GetExitEnergy, &PROPOSALParticle::SetExitEnergy)
        .add_property("closet_approach_point", &PROPOSALParticle::GetClosestApproachPoint, &PROPOSALParticle::SetClosestApproachPoint)
        .add_property("closet_approach_time", &PROPOSALParticle::GetClosestApproachTime, &PROPOSALParticle::SetClosestApproachTime)
        .add_property("closet_approach_energy", &PROPOSALParticle::GetClosestApproachEnergy, &PROPOSALParticle::SetClosestApproachEnergy)

        .add_property("e_lost", &PROPOSALParticle::GetElost, &PROPOSALParticle::SetElost)
    ;

    /**************************************************************************
    *                                 Medium                                 *
    **************************************************************************/

    // --------------------------------------------------------------------- //
    // Medium Definition
    // --------------------------------------------------------------------- //

    enum_<MediumFactory::Enum>("MediumType")
        .value("Water", MediumFactory::Water)
        .value("Ice", MediumFactory::Ice)
        .value("Salt", MediumFactory::Salt)
        .value("StandardRock", MediumFactory::StandardRock)
        .value("FrejusRock", MediumFactory::FrejusRock)
        .value("Iron", MediumFactory::Iron)
        .value("Hydrogen", MediumFactory::Hydrogen)
        .value("Lead", MediumFactory::Lead)
        .value("Copper", MediumFactory::Copper)
        .value("Air", MediumFactory::Air)
        .value("Paraffin", MediumFactory::Paraffin)
        .value("AntaresWater", MediumFactory::AntaresWater);

    class_<MediumFactory::Definition, boost::shared_ptr<MediumFactory::Definition> >("MediumDefinition", init<>())

        .def_readwrite("type", &MediumFactory::Definition::type)
        .def_readwrite("density_correction", &MediumFactory::Definition::density_correction)
    ;

    /**************************************************************************
    *                                Geometry                                 *
    **************************************************************************/

    // --------------------------------------------------------------------- //
    // Geometry Definition
    // --------------------------------------------------------------------- //

    enum_<GeometryFactory::Enum>("Shape")
        .value("Sphere", GeometryFactory::Sphere)
        .value("Box", GeometryFactory::Box)
        .value("Cylinder", GeometryFactory::Cylinder);

    class_<GeometryFactory::Definition, boost::shared_ptr<GeometryFactory::Definition> >("GeometryDefinition", init<>())

        .def_readwrite("shape", &GeometryFactory::Definition::shape)
        .def_readwrite("position", &GeometryFactory::Definition::position)
        .def_readwrite("inner_radius", &GeometryFactory::Definition::inner_radius)
        .def_readwrite("radius", &GeometryFactory::Definition::radius)
        .def_readwrite("width", &GeometryFactory::Definition::width)
        .def_readwrite("height", &GeometryFactory::Definition::height)
        .def_readwrite("depth", &GeometryFactory::Definition::depth)
    ;

    // --------------------------------------------------------------------- //
    // Geometry
    // --------------------------------------------------------------------- //

    class_<Geometry, boost::shared_ptr<Geometry>, boost::noncopyable>("Geometry", no_init)

        .def(self_ns::str(self_ns::self))

        .def("is_infront", &Geometry::IsInfront)
        .def("is_inside", &Geometry::IsInside)
        .def("is_behind", &Geometry::IsBehind)
        .def("distance_to_border", &Geometry::DistanceToBorder)
        .def("distance_to_closet_approach", &Geometry::DistanceToClosestApproach)

        .add_property("name", &Geometry::GetName)
        .add_property("position", &Geometry::GetPosition, &Geometry::SetPosition)
        .add_property("hirarchy", &Geometry::GetHirarchy, &Geometry::SetHirarchy)
    ;

    class_<Sphere, boost::shared_ptr<Sphere>, bases<Geometry> >("Sphere", init<>())

        .def(init<Vector3D, double, double>())
        .def(init<const Sphere&>())

        .add_property("inner_radius", &Sphere::GetInnerRadius, &Sphere::SetInnerRadius)
        .add_property("radius", &Sphere::GetRadius, &Sphere::SetRadius)
    ;

    class_<Box, boost::shared_ptr<Box>, bases<Geometry> >("Box", init<>())

        .def(init<Vector3D, double, double, double>())
        .def(init<const Box&>())

        .add_property("width", &Box::GetX, &Box::SetX)
        .add_property("height", &Box::GetY, &Box::SetY)
        .add_property("depth", &Box::GetZ, &Box::SetZ)
    ;

    class_<Cylinder, boost::shared_ptr<Cylinder>, bases<Geometry> >("Cylinder", init<>())

        .def(init<Vector3D, double, double, double>())
        .def(init<const Cylinder&>())

        .add_property("inner_radius", &Cylinder::GetInnerRadius, &Cylinder::SetInnerRadius)
        .add_property("radius", &Cylinder::GetRadius, &Cylinder::SetRadius)
        .add_property("height", &Cylinder::GetZ, &Cylinder::SetZ)
    ;

    // --------------------------------------------------------------------- //
    // Scattering Definition
    // --------------------------------------------------------------------- //

    enum_<ScatteringFactory::Enum>("ScatteringModel")
        .value("default", ScatteringFactory::Default)
        .value("moliere", ScatteringFactory::Moliere)
        .value("moliere_first_order", ScatteringFactory::MoliereFirstOrder);

    // --------------------------------------------------------------------- //
    // Bremsstrahlung Definition
    // --------------------------------------------------------------------- //

    enum_<BremsstrahlungFactory::Enum>("BremsParametrization")
        .value("PetrukhinShestakov", BremsstrahlungFactory::PetrukhinShestakov)
        .value("KelnerKokoulinPetrukhin", BremsstrahlungFactory::KelnerKokoulinPetrukhin)
        .value("CompleteScreening", BremsstrahlungFactory::CompleteScreening)
        .value("AndreevBezrukovBugaev", BremsstrahlungFactory::AndreevBezrukovBugaev);

    class_<BremsstrahlungFactory::Definition, boost::shared_ptr<BremsstrahlungFactory::Definition> >("BremsDefinition", init<>())

        .def_readwrite("parametrization", &BremsstrahlungFactory::Definition::parametrization)
        .def_readwrite("lpm_effect", &BremsstrahlungFactory::Definition::lpm_effect)
        .def_readwrite("multiplier", &BremsstrahlungFactory::Definition::multiplier)
    ;

    // --------------------------------------------------------------------- //
    // Photonuclear  Definition
    // --------------------------------------------------------------------- //

    enum_<PhotonuclearFactory::Enum>("PhotoParametrization")
        .value("Zeus", PhotonuclearFactory::Zeus)
        .value("BezrukovBugaev", PhotonuclearFactory::BezrukovBugaev)
        .value("Rhode", PhotonuclearFactory::Rhode)
        .value("Kokoulin", PhotonuclearFactory::Kokoulin)
        .value("AbramowiczLevinLevyMaor91", PhotonuclearFactory::AbramowiczLevinLevyMaor91)
        .value("AbramowiczLevinLevyMaor97", PhotonuclearFactory::AbramowiczLevinLevyMaor97)
        .value("ButkevichMikhailov", PhotonuclearFactory::ButkevichMikhailov)
        .value("RenoSarcevicSu", PhotonuclearFactory::RenoSarcevicSu);

    enum_<PhotonuclearFactory::Shadow>("PhotoShadow")
        .value("Dutta", PhotonuclearFactory::ShadowDutta)
        .value("ButkevichMikhailov", PhotonuclearFactory::ShadowButkevichMikhailov);


    class_<PhotonuclearFactory::Definition, boost::shared_ptr<PhotonuclearFactory::Definition> >("PhotoDefinition", init<>())

        .def_readwrite("parametrization", &PhotonuclearFactory::Definition::parametrization)
        .def_readwrite("shadow", &PhotonuclearFactory::Definition::shadow)
        .def_readwrite("hardbb", &PhotonuclearFactory::Definition::hardbb)
        .def_readwrite("multiplier", &PhotonuclearFactory::Definition::multiplier)
    ;

    // --------------------------------------------------------------------- //
    // EpairProduction Definition
    // --------------------------------------------------------------------- //

    class_<EpairProductionFactory::Definition, boost::shared_ptr<EpairProductionFactory::Definition> >("EpairDefinition", init<>())

        .def_readwrite("lpm_effect", &EpairProductionFactory::Definition::lpm_effect)
        .def_readwrite("multiplier", &EpairProductionFactory::Definition::multiplier)
    ;

    // --------------------------------------------------------------------- //
    // Ionization Definition
    // --------------------------------------------------------------------- //

    class_<IonizationFactory::Definition, boost::shared_ptr<IonizationFactory::Definition> >("IonizationDefinition", init<>())

        .def_readwrite("multiplier", &IonizationFactory::Definition::multiplier)
    ;

    // --------------------------------------------------------------------- //
    // Utility Definition
    // --------------------------------------------------------------------- //


    class_<Utility::Definition, boost::shared_ptr<Utility::Definition> >("UtilityDefinition", init<>())

        .def_readwrite("brems_def", &Utility::Definition::brems_def)
        .def_readwrite("photo_def", &Utility::Definition::photo_def)
        .def_readwrite("epair_def", &Utility::Definition::epair_def)
        .def_readwrite("ioniz_def", &Utility::Definition::ioniz_def)
    ;

    // --------------------------------------------------------------------- //
    // Sector Definition
    // --------------------------------------------------------------------- //

    // ----[ Location ]-------------------------------------- //

    enum_<Sector::ParticleLocation::Enum>("ParticleLocation")
        .value("infront_detector", Sector::ParticleLocation::InfrontDetector)
        .value("inside_detector", Sector::ParticleLocation::InsideDetector)
        .value("behind_detector", Sector::ParticleLocation::BehindDetector);

    class_<SectorFactory::Definition, boost::shared_ptr<SectorFactory::Definition> >("SectorDefinition", init<>())

        .def_readwrite("e_cut", &SectorFactory::Definition::e_cut)
        .def_readwrite("v_cut", &SectorFactory::Definition::v_cut)
        .def_readwrite("medium_def", &SectorFactory::Definition::medium_def)
        .def_readwrite("geometry_def", &SectorFactory::Definition::geometry_def)
        .def_readwrite("do_weighting", &SectorFactory::Definition::do_weighting)
        .def_readwrite("weighting_order", &SectorFactory::Definition::weighting_order)
        .def_readwrite("do_continuous_randomization", &SectorFactory::Definition::do_continuous_randomization)
        .def_readwrite("do_exact_time_calculation", &SectorFactory::Definition::do_exact_time_calculation)
        .def_readwrite("scattering_model", &SectorFactory::Definition::scattering_model)
        .def_readwrite("particle_location", &SectorFactory::Definition::location)
        .def_readwrite("crosssection_defs", &SectorFactory::Definition::utility_def)
    ;

    // --------------------------------------------------------------------- //
    // Interpolation Definition
    // --------------------------------------------------------------------- //

    class_<InterpolationDef, boost::shared_ptr<InterpolationDef> >("InterpolationDef", init<>())

        // .def(self_ns::str(self_ns::self))

        .def_readwrite("order_of_interpolation", &InterpolationDef::order_of_interpolation)
        .def_readwrite("path_to_tables", &InterpolationDef::path_to_tables)
        .def_readwrite("raw", &InterpolationDef::raw)
    ;


    // --------------------------------------------------------------------- //
    // Propagator
    // --------------------------------------------------------------------- //

    class_<Propagator, boost::shared_ptr<Propagator> >(
        "Propagator", init<const ParticleDef&, const std::vector<SectorFactory::Definition>&, const Geometry&>((arg("partcle_def"), arg("sector_defs"), arg("detector"))))

        .def(init<const ParticleDef&,
                  const std::vector<SectorFactory::Definition>&,
                  const Geometry&,
                  const InterpolationDef&>((
            arg("particle_def"), arg("sector_defs"), arg("detector"), arg("interpolation_def"))))

        .def(init<const ParticleDef&,
                  const std::string&>((
            arg("particle_def"), arg("config_file"))))

        // .def(self_ns::str(self_ns::self))

        .def("propagate", &Propagator::Propagate, (arg("max_distance_cm") = 1e20))
        .add_property("particle", make_function(&Propagator::GetParticle, return_internal_reference<>()), "Get the internal created particle to modify its properties");

    // --------------------------------------------------------------------- //
    // PropagatorService
    // --------------------------------------------------------------------- //

    class_<PropagatorService, boost::shared_ptr<PropagatorService> >(
        "PropagatorService", init<>())

        .def("propagate", &PropagatorService::Propagate, (arg("particle")))
        .def("register_propagator", &PropagatorService::RegisterPropagator, (arg("propagator")));

    // --------------------------------------------------------------------- //
    // HardBBTable
    // --------------------------------------------------------------------- //

    boost::python::scope hardbb = class_<PythonHardBBTables>("HardBBTables");
    hardbb.attr("MuonTable") = HardBBTables::MuonTable;
    hardbb.attr("TauTable") = HardBBTables::TauTable;
}

#undef PARTICLE_DEF
#undef COMPONENT_DEF
#undef MEDIUM_DEF
#undef BREMS_DEF
#undef PHOTO_REAL_DEF
#undef PHOTO_Q2_DEF
#undef PHOTO_Q2_INTERPOL_DEF
