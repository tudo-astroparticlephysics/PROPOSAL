#include <I3Test.h>

#include <PROPOSAL-icetray/I3PropagatorServicePROPOSAL.h>
#include <PROPOSAL-icetray/SimplePropagator.h>
#include <dataclasses/physics/I3Particle.h>
#include <phys-services/I3SPRNGRandomService.h>

#include <boost/make_shared.hpp>

TEST_GROUP(Repeatablility);

static I3Particle make_particle()
{
    I3Particle p;
    p.SetPos(0, 0, 2e3);
    p.SetDir(0, 0);
    p.SetTime(0);
    p.SetType(I3Particle::MuMinus);
    p.SetLocationType(I3Particle::InIce);
    p.SetEnergy(1e6);

    return p;
}

TEST(SimplePropagator)
{
    I3RandomServicePtr rng(new I3SPRNGRandomService(1, 10000, 2));

    boost::shared_ptr<PROPOSAL::SimplePropagator> prop(
        new PROPOSAL::SimplePropagator(I3Particle::MuMinus, "ice", 5e2, -1));
    prop->SetRandomNumberGenerator(rng);
    double distance = 1e2;

    std::vector<std::vector<I3Particle> > daughters;
    std::vector<I3Particle> primaries;
    std::vector<I3FrameObjectPtr> states;
    for (int i = 0; i < 100; i++)
    {
        I3Particle p = make_particle();
        p.SetEnergy(std::pow(10, rng->Uniform(3, 6)));
        primaries.push_back(p);
        states.push_back(rng->GetState());
        boost::shared_ptr<std::vector<I3Particle> > d(new std::vector<I3Particle>);
        prop->propagate(p, distance, d);
        daughters.push_back(*d);
    }

    // Now, replay the simulation and ensure we get the same thing

    for (size_t i = 0; i < daughters.size(); i++)
    {
        // Throw out a random subset of events. The remainder should
        // be reproducible given the same sequence of random numbers,
        // unless of course the propagator keeps secret state.
        if (rng->Uniform() < 0.5)
            continue;
        rng->RestoreState(states[i]);
        boost::shared_ptr<std::vector<I3Particle> > d(new std::vector<I3Particle>);
        prop->propagate(primaries[i], distance, d);
        ENSURE_EQUAL(daughters[i].size(), d->size(), "Same number of particles must be produced");
        for (size_t j = 0; j < daughters[i].size(); j++)
        {
            I3Particle& p1 = daughters[i][j];
            I3Particle& p2 = (*d)[j];
            ENSURE_EQUAL(p1.GetType(), p2.GetType(), "Secondaries must have the same type");
            ENSURE_EQUAL(p1.GetTime(), p2.GetTime(), "Times must be completely identical");
            ENSURE_EQUAL(p1.GetEnergy(), p2.GetEnergy(), "Energies must completely identical");
        }
    }
}

TEST(PropagatorService)
{
    I3RandomServicePtr rng(new I3SPRNGRandomService(1, 10000, 1));
    I3FrameObjectPtr state = rng->GetState();

    PROPOSAL::I3PropagatorServicePROPOSALPtr prop(new PROPOSAL::I3PropagatorServicePROPOSAL);
    prop->SetRandomNumberGenerator(rng);

    I3PropagatorService::DiagnosticMapPtr frame(new I3PropagatorService::DiagnosticMap);
    // the dummy I3Frame makes compiler happy, but won't be used
    I3FramePtr dummy(new I3Frame()); 

    std::vector<std::vector<I3Particle> > daughters;
    for (int i = 0; i < 2; i++)
    {
        I3Particle p              = make_particle();
        std::vector<I3Particle> d = prop->Propagate(p, frame, dummy);
        daughters.push_back(d);
    }

    rng->RestoreState(state);

    for (size_t i = 0; i < daughters.size(); i++)
    {
        I3Particle p              = make_particle();
        std::vector<I3Particle> d = prop->Propagate(p, frame, dummy);
        ENSURE_EQUAL(daughters[i].size(), d.size());
        for (size_t j = 0; j < daughters[i].size(); j++)
        {
            I3Particle& p1 = daughters[i][j];
            I3Particle& p2 = d[j];
            ENSURE_EQUAL(p1.GetType(), p2.GetType());
            ENSURE_EQUAL(p1.GetTime(), p2.GetTime());
            ENSURE_EQUAL(p1.GetEnergy(), p2.GetEnergy());
        }
    }
}
