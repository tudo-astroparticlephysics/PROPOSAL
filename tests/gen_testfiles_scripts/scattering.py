import pyPROPOSAL as pp
import numpy as np

parametrizations = [
    pp.Scattering.Moliere,
    pp.Scattering.HighlandIntegral,
    pp.Scattering.Highland,
    pp.Scattering.NoScattering
]

scat_names = [
    "Moliere",
    "HighlandIntegral",
    "Highland",
    "NoScattering"
]

particle_defs = [
    pp.MuMinusDef.get(),
    pp.TauMinusDef.get(),
    pp.EMinusDef.get()
]

mediums = [
    pp.Medium.Ice(1.0),
    pp.Medium.Hydrogen(1.0),
    pp.Medium.Uranium(1.0)
]

cuts = [
    pp.EnergyCutSettings(-1, -1),
    pp.EnergyCutSettings(500, -1),
    pp.EnergyCutSettings(-1, 0.05),
    pp.EnergyCutSettings(500, 0.05)
]

energies = np.logspace(4, 13, num=10) # MeV
distance = 1000 # cm
position_init = pp.Vector3D(0, 0, 0)
direction_init = pp.Vector3D(1, 0, 0)
direction_init.spherical_from_cartesian()

interpoldef = pp.InterpolationDef()
pp.RandomGenerator.get().set_seed(1234)


def create_table_scatter():

    with open("TestFiles/Scat_scatter.txt", "a") as f:

        for particle_def in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for idx, parametrization in enumerate(parametrizations):
                        particle = pp.Particle(particle_def)

                        if scat_names[idx] == "HighlandIntegral":
                            util = pp.Utility(particle_def, medium, cut, pp.UtilityDefinition(), pp.InterpolationDef())
                            scattering = parametrization(particle, util, pp.InterpolationDef())
                        else:
                            scattering = parametrization(particle, medium)

                        for energy in energies:
                            # TODO wrap UtilityDecorator
                            particle.energy = energy
                            particle.position = position_init
                            particle.direction = direction_init

                            energy_out = energy - 1000
                            scattering.scatter(distance, energy, energy_out)

                            buf = [""]

                            buf.append(particle_def.name)
                            buf.append(medium.name)
                            buf.append(scat_names[idx])
                            buf.append(str(cut.ecut))
                            buf.append(str(cut.vcut))
                            buf.append(str(energy))
                            buf.append(str(energy_out))
                            buf.append(str(distance))
                            buf.append(str(particle.position.x))
                            buf.append(str(particle.position.y))
                            buf.append(str(particle.position.z))
                            buf.append(str(particle.direction.radius))
                            buf.append(str(particle.direction.phi))
                            buf.append(str(particle.direction.theta))
                            buf.append("\n")

                            f.write("\t".join(buf))
#

if __name__ == "__main__":
    create_table_scatter()
