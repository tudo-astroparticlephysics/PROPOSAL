import pyPROPOSAL

ptype = pyPROPOSAL.ParticleType.MuMinus
mu = pyPROPOSAL.Particle(ptype)
mu.energy = 1e10

p = pyPROPOSAL.Propagator(
    config="resources/configuration",
    particle=mu
)

mu_length = list()
n_daughters = list()

for i in range(1000):
    p.reset_particle()
    d = p.propagate()

    mu_length.append(mu.propagated_distance)
    # mu_length.append(p.GetParticle().GetPropagatedDistance())
    n_daughters.append(len(d))


try:
    import matplotlib as mpl
    mpl.use('Agg')
    import pylab

    pylab.figure()
    pylab.title("Mu Lengths")
    pylab.hist(mu_length, histtype="step", log=True, bins=100)
    pylab.xlabel(r'$l_{\mu}(\rm{m})$')
    pylab.savefig('MuonLenghts.pdf')

    pylab.figure()
    pylab.title("N Daughters")
    pylab.hist(n_daughters, histtype="step", log=True, bins=100)
    pylab.xlabel('N')
    pylab.savefig('Daughters.pdf')
except ImportError :
    print("pylab not installed.  no plots for you.")
