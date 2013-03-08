#include <I3Test.h>

#include <dataclasses/physics/I3Particle.h>
#include <phys-services/I3SPRNGRandomService.h>
#include <PROPOSAL-icetray/I3PropagatorServicePROPOSAL.h>
#include <PROPOSAL-icetray/SimplePropagator.h>

#include <boost/make_shared.hpp>

TEST_GROUP(Repeatablility);

static I3Particle make_particle()
{
	I3Particle p;
	p.SetPos(0,0,2e3);
	p.SetDir(0,0);
	p.SetTime(0);
	p.SetType(I3Particle::MuMinus);
	p.SetLocationType(I3Particle::InIce);
	p.SetEnergy(1e6);
	
	return p;
}

TEST(SimplePropagator)
{
	I3RandomServicePtr rng(new I3SPRNGRandomService(1,10000,1));
	I3FrameObjectPtr state = rng->GetState();
	
	boost::shared_ptr<PROPOSAL::SimplePropagator> prop(new PROPOSAL::SimplePropagator("ice", 5e2, -1));
	prop->SetRandomNumberGenerator(rng);
	double distance = 1e1;
	
	std::vector<std::vector<I3Particle> > daughters;
	for (int i=0; i<1; i++) {
		I3Particle p = make_particle();
		boost::shared_ptr<std::vector<I3Particle> > d(new std::vector<I3Particle>);
		prop->propagate(p, distance, d);
		daughters.push_back(*d);
	}
	
	rng->RestoreState(state);
	// prop = boost::make_shared<PROPOSAL::SimplePropagator>("ice", 5e2, -1);
	// prop->SetRandomNumberGenerator(rng);
	
	for (int i=0; i<daughters.size(); i++) {
		I3Particle p = make_particle();
		boost::shared_ptr<std::vector<I3Particle> > d(new std::vector<I3Particle>);
		prop->propagate(p, distance, d);
		ENSURE_EQUAL(daughters[i].size(), d->size());
		for (int j=0; j < daughters[i].size(); j++) {
			I3Particle &p1 = daughters[i][j];
			I3Particle &p2 = (*d)[j];
			ENSURE_EQUAL(p1.GetType(), p2.GetType());
			ENSURE_EQUAL(p1.GetTime(), p2.GetTime());
			ENSURE_EQUAL(p1.GetEnergy(), p2.GetEnergy());
		}
	}
	
}

TEST(PropagatorService)
{
	I3RandomServicePtr rng(new I3SPRNGRandomService(1,10000,1));
	I3FrameObjectPtr state = rng->GetState();
	
	I3PropagatorServicePROPOSALPtr prop(new I3PropagatorServicePROPOSAL);
	prop->SetRandomNumberGenerator(rng);
	
	I3FramePtr frame(new I3Frame);
	
	std::vector<std::vector<I3Particle> > daughters;
	for (int i=0; i<1; i++) {
		I3Particle p = make_particle();
		std::vector<I3Particle> d = prop->Propagate(p, frame);
		daughters.push_back(d);
	}
	
	rng->RestoreState(state);
	// prop = I3PropagatorServicePROPOSALPtr(new I3PropagatorServicePROPOSAL);
	// prop->SetRandomNumberGenerator(rng);
	
	for (int i=0; i<daughters.size(); i++) {
		I3Particle p = make_particle();
		std::vector<I3Particle> d = prop->Propagate(p, frame);
		ENSURE_EQUAL(daughters[i].size(), d.size());
		for (int j=0; j < daughters[i].size(); j++) {
			I3Particle &p1 = daughters[i][j];
			I3Particle &p2 = d[j];
			ENSURE_EQUAL(p1.GetType(), p2.GetType());
			ENSURE_EQUAL(p1.GetTime(), p2.GetTime());
			ENSURE_EQUAL(p1.GetEnergy(), p2.GetEnergy());
		}
	}
	
}