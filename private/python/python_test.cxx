#include <boost/python.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <string>

#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/EnergyCutSettings.h"


template<typename T>
struct VectorToPythonList
{
    // static PyObject* convert(std::vector<T> const& v)
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
        void* storage = reinterpret_cast<boost::python::converter::rvalue_from_python_storage<std::vector<T>>*>(data)->storage.bytes;

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

/******************************************************************************
*                               Python Module                                *
******************************************************************************/

BOOST_PYTHON_MODULE(pyPROPOSAL)
{

    using namespace boost::python;

    // --------------------------------------------------------------------- //
    // Vector classes
    // --------------------------------------------------------------------- //

    // register the to-python converter
    to_python_converter< std::vector<double>, VectorToPythonList<double> >();
    to_python_converter< std::vector<std::string>, VectorToPythonList<std::string> >();
    to_python_converter< std::vector<PROPOSALParticle*>, VectorToPythonList<PROPOSALParticle*> >();

    // register the from-python converter
    VectorFromPythonList<double>();
    VectorFromPythonList<std::string>();
    VectorFromPythonList<PROPOSALParticle*>();

    // class_<std::vector<double>>("DoubelVec")
    //     .def(vector_indexing_suite<std::vector<double>>())
    //     ;
    //
    // class_<std::vector<std::string>>("StringVec")
    //     .def(vector_indexing_suite<std::vector<std::string>>())
    //     ;
    //
    // class_<std::vector<PROPOSALParticle*>>("Secondarys")
    //     .def(vector_indexing_suite<std::vector<PROPOSALParticle*>>())
    //     ;


    // --------------------------------------------------------------------- //
    // ParticleType
    // --------------------------------------------------------------------- //


    enum_<PROPOSALParticle::ParticleType>("ParticleType")
            .value("EPlus",                PROPOSALParticle::ParticleType::EPlus)
            .value("EMinus",               PROPOSALParticle::ParticleType::EPlus)
            .value("MuPlus",               PROPOSALParticle::ParticleType::MuPlus)
            .value("MuMinus",              PROPOSALParticle::ParticleType::MuMinus)
            .value("TauPlus",              PROPOSALParticle::ParticleType::TauPlus)
            .value("TauMinus",             PROPOSALParticle::ParticleType::TauMinus)
            .value("NuE",                  PROPOSALParticle::ParticleType::NuE)
            .value("NuEBar",               PROPOSALParticle::ParticleType::NuEBar)
            .value("NuMu",                 PROPOSALParticle::ParticleType::NuMu)
            .value("NuMuBar",              PROPOSALParticle::ParticleType::NuMuBar)
            .value("NuTau",                PROPOSALParticle::ParticleType::NuTau)
            .value("NuTauBar",             PROPOSALParticle::ParticleType::NuTauBar)
            .value("Brems",                PROPOSALParticle::ParticleType::Brems)
            .value("DeltaE",               PROPOSALParticle::ParticleType::DeltaE)
            .value("EPair",                PROPOSALParticle::ParticleType::EPair)
            .value("NuclInt",              PROPOSALParticle::ParticleType::NuclInt)
            .value("MuPair",               PROPOSALParticle::ParticleType::MuPair)
            .value("Hadrons",              PROPOSALParticle::ParticleType::Hadrons)
            .value("ContinuousEnergyLoss", PROPOSALParticle::ParticleType::ContinuousEnergyLoss)
            .value("Monopole",             PROPOSALParticle::ParticleType::Monopole)
            .value("STauPlus",             PROPOSALParticle::ParticleType::STauPlus)
            .value("STauMinus",            PROPOSALParticle::ParticleType::STauMinus)
            ;

    // --------------------------------------------------------------------- //
    // Particle
    // --------------------------------------------------------------------- //

    std::string (PROPOSALParticle::*getNameParticle)() const = &PROPOSALParticle::GetName;

    class_<PROPOSALParticle, boost::shared_ptr<PROPOSALParticle>>("Particle",
                                                                  init<PROPOSALParticle::ParticleType>(
                                                                  (arg("type")=PROPOSALParticle::ParticleType::MuMinus)))

            .def(self_ns::str(self_ns::self))

            .add_property("energy", &PROPOSALParticle::GetEnergy, &PROPOSALParticle::SetEnergy)
            .add_property("propagated_distance", &PROPOSALParticle::GetPropagatedDistance, &PROPOSALParticle::SetPropagatedDistance)
            .add_property("X", &PROPOSALParticle::GetX, &PROPOSALParticle::SetX)
            .add_property("Y", &PROPOSALParticle::GetY, &PROPOSALParticle::SetY)
            .add_property("Z", &PROPOSALParticle::GetY, &PROPOSALParticle::SetZ)
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

    class_<EnergyCutSettings, boost::shared_ptr<EnergyCutSettings>>("EnergyCutSettings",
                                              init<double, double>(
                                              (arg("ecut"),
                                               arg("vcut"))))

            .def(self_ns::str(self_ns::self))

            .add_property("ecut", &EnergyCutSettings::GetEcut, &EnergyCutSettings::SetEcut)
            .add_property("vcut", &EnergyCutSettings::GetVcut, &EnergyCutSettings::SetVcut)

            .def("get_cut", &EnergyCutSettings::GetCut, "Return the lower from E*v = e")
        ;

    // --------------------------------------------------------------------- //
    // Medium
    // --------------------------------------------------------------------- //

    std::vector<std::string> (Medium::*getNameMed)() const = &Medium::GetElementName;

    class_<Medium, boost::shared_ptr<Medium>>("Medium",
                                              init<std::string, double>(
                                              (arg("type"),
                                               arg("rho")=1.0)))

            .def(self_ns::str(self_ns::self))

            .add_property("num_components", &Medium::GetNumComponents, &Medium::SetNumComponents)
            .add_property("nuc_charge", &Medium::GetNucCharge, &Medium::SetNucCharge)
            .add_property("atomic_num", &Medium::GetAtomicNum, &Medium::SetAtomicNum)
            .add_property("atom_in_molecule", &Medium::GetAtomInMolecule, &Medium::SetAtomInMolecule)
            .add_property("sum_charge", &Medium::GetSumCharge, &Medium::SetSumCharge)
            .add_property("ZA", &Medium::GetZA, &Medium::SetZA)
            .add_property("I", &Medium::GetI, &Medium::SetI)
            .add_property("C1", &Medium::GetC1, &Medium::SetC1)
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
            .add_property("MN", &Medium::GetMN, &Medium::SetMN)
            .add_property("MM", &Medium::GetMM, &Medium::SetMM)
            .add_property("sum_nucleons", &Medium::GetSumNucleons, &Medium::SetSumNucleons)
            .add_property("R0", &Medium::GetR0, &Medium::SetR0)
        ;

    // --------------------------------------------------------------------- //
    // Propagator
    // --------------------------------------------------------------------- //

    class_<Propagator, boost::shared_ptr<Propagator>>("Propagator",
                                                      init<std::string, PROPOSALParticle*, bool>(
                                                      (arg("config"),
                                                       arg("particle"),
                                                       arg("applyoptions")=true)))

            .def("propagate", &Propagator::propagate, (arg("max_distance_cm") = 1e20))
            .def("apply_options", &Propagator::ApplyOptions)
            .def("reset_particle", &Propagator::ResetParticle)

            // .add_property("particle", &Propagator::SetParticle)
            .add_property("particle", make_function(&Propagator::GetParticle, return_value_policy<reference_existing_object>()), &Propagator::SetParticle)
            .add_property("seed",&Propagator::GetSeed ,&Propagator::SetSeed)
            .add_property("brems",&Propagator::GetBrems ,&Propagator::SetBrems)
            .add_property("photo",&Propagator::GetPhoto ,&Propagator::SetPhoto)
            .add_property("path_to_tables",&Propagator::GetPath_to_tables ,&Propagator::SetPath_to_tables)
            .add_property("stopping_decay",&Propagator::GetStopping_decay ,&Propagator::SetStopping_decay)
        ;
}
