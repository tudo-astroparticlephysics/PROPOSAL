import pyPROPOSAL as pp
import numpy as np

parametrizations = [
    pp.scattering.Moliere,
    pp.scattering.HighlandIntegral,
    pp.scattering.Highland,
    pp.scattering.NoScattering
]

scat_names = [
    "Moliere",
    "HighlandIntegral",
    "Highland",
    "NoScattering"
]

particle_defs = [
    pp.particle.MuMinusDef.get(),
    pp.particle.TauMinusDef.get(),
    pp.particle.EMinusDef.get()
]

mediums = [
    pp.medium.Ice(1.0),
    pp.medium.Hydrogen(1.0),
    pp.medium.Uranium(1.0)
]

cuts = [
    pp.EnergyCutSettings(-1, -1),
    pp.EnergyCutSettings(500, -1),
    pp.EnergyCutSettings(-1, 0.05),
    pp.EnergyCutSettings(500, 0.05)
]

energies = np.logspace(4, 13, num=10)  # MeV
energies_out = np.logspace(3, 12, num=10)
distance = 1000  # cm
position_init = pp.Vector3D(0, 0, 0)
direction_init = pp.Vector3D(1, 0, 0)
direction_init.spherical_from_cartesian()

interpoldef = pp.InterpolationDef()
pp.RandomGenerator.get().set_seed(1234)


def create_table_scatter(dir_name):

    with open(dir_name + "Scattering_scatter.txt", "a") as file:

        for particle_def in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for idx, parametrization in enumerate(parametrizations):
                        particle = pp.particle.Particle(particle_def)

                        if scat_names[idx] == "HighlandIntegral":
                            util = pp.Utility(particle_def, medium, cut, pp.UtilityDefinition(), pp.InterpolationDef())
                            scattering = parametrization(particle, util, pp.InterpolationDef())
                        else:
                            scattering = parametrization(particle, medium)

                        buf = [""]
                        for jdx, energy in enumerate(energies):
                            # TODO wrap UtilityDecorator
                            particle.energy = energy
                            particle.position = position_init
                            particle.direction = direction_init

                            scattering.scatter(distance, energy, energies_out[jdx])


                            buf.append(particle_def.name)
                            buf.append(medium.name)
                            buf.append(scat_names[idx])
                            buf.append(str(cut.ecut))
                            buf.append(str(cut.vcut))
                            buf.append(str(energy))
                            buf.append(str(energies_out[jdx]))
                            buf.append(str(distance))
                            buf.append(str(particle.position.x))
                            buf.append(str(particle.position.y))
                            buf.append(str(particle.position.z))
                            buf.append(str(particle.direction.radius))
                            buf.append(str(particle.direction.phi))
                            buf.append(str(particle.direction.theta))
                            buf.append("\n")

                        file.write("\t".join(buf))


def main(dir_name):
    create_table_scatter(dir_name)


if __name__ == "__main__":

    import os

    dir_name = "TestFiles/"

    if os.path.isdir(dir_name):
        print("Directory {} already exists".format(dir_name))
    else:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))

    main(dir_name)
