import pyPROPOSAL

energy = 1e8
statistics = 10000

ptype = pyPROPOSAL.ParticleType.MuMinus

med = pyPROPOSAL.Medium("standard_rock")
E = pyPROPOSAL.EnergyCutSettings()
p = pyPROPOSAL.Propagator(
    med,
    E,
    ptype,
    "resources/tables",
    moliere=False,
    scattering_model=2
)

# mu = pyPROPOSAL.Particle(ptype)
# mu.energy = energy
# mu.theta = 180
# p = pyPROPOSAL.Propagator(
#     config="resources/configuration",
#     particle=mu
# )

mu_length = list()
n_daughters = list()

for i in range(statistics):
    p.reset_particle()
    p.particle.energy = energy
    d = p.propagate()

    mu_length.append(p.particle.propagated_distance / 100)
    # mu_length.append(p.GetParticle().GetPropagatedDistance())
    n_daughters.append(len(d))


try:
    import matplotlib as mpl
    mpl.use('Agg')
    import pylab

    pylab.figure()
    pylab.title("{} muons with energy {} TeV in {}".format(
        statistics,
        energy / 1e6,
        p.collections[0].medium.name.lower()
    ))
    pylab.hist(mu_length, histtype="step", log=True, bins=100)
    pylab.xlabel(r'range / m')
    pylab.ylabel(r'count')
    pylab.savefig('MuonLenghts.pdf')

    pylab.figure()
    pylab.title("N Daughters")
    pylab.hist(n_daughters, histtype="step", log=True, bins=100)
    pylab.xlabel('N')
    pylab.savefig('Daughters.pdf', bbox_inches='tight', pad_inches=0.02)
except ImportError :
    print("pylab not installed.  no plots for you.")
