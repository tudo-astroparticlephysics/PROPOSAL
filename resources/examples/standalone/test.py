import pyPROPOSAL

med = pyPROPOSAL.Medium("ice")
E = pyPROPOSAL.EnergyCutSettings(500, 1)
ptype = pyPROPOSAL.ParticleType.MuMinus

p = pyPROPOSAL.Propagator(med, E, ptype, "")
p.particle.energy = 1e12

d = p.propagate()
print(len(d))

# a = pr.cross_sections

# p = pyPROPOSAL.Propagator(med, E, ptype, "resources/tables")
# print(len(p.current_collection.cross_sections))

# mu_length = list()
# n_daughters = list()
#
# for i in range(10):
#     p.reset_particle()
#     p.particle.energy = 1e10
#
#     d = p.propagate()
#
#     mu_length.append(p.particle.propagated_distance)
#     n_daughters.append(len(d))

# try:
#     import matplotlib as mpl
#     mpl.use('Agg')
#     import pylab
#
#     pylab.figure()
#     pylab.title("Mu Lengths")
#     pylab.hist(mu_length, histtype="step", log=True, bins=100)
#     pylab.xlabel(r'$l_{\mu}(\rm{m})$')
#     pylab.savefig('MuonLenghts.pdf')
#
#     pylab.figure()
#     pylab.title("N Daughters")
#     pylab.hist(n_daughters, histtype="step", log=True, bins=100)
#     pylab.xlabel('N')
#     pylab.savefig('Daughters.pdf')
# except ImportError :
#     print("pylab not installed.  no plots for you.")
