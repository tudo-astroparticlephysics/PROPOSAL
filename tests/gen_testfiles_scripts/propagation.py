import os
import proposal as pp


def create_table(dir_name):
    """TODO: Docstring for create_table.
    Returns: TODO

    """

    statistics = 10

    prop = pp.Propagator(
        pp.particle.MuMinusDef(),
        "examples/config_minimal.json"
    )

    init_particle = pp.particle.ParticleState()

    init_particle.energy = 1e8
    init_particle.propagated_distance = 0
    init_particle.position = pp.Cartesian3D(0, 0, 0)
    init_particle.direction = pp.Cartesian3D(0, 0, -1)
    init_particle.time = 0
    pp.RandomGenerator.get().set_seed(1234)

    with open(dir_name + "Propagator_propagation.txt", "w") as file:

        buf = [""]
        buf.append("name")
        buf.append("length")
        buf.append("energy")
        buf.append("x")
        buf.append("y")
        buf.append("z")
        buf.append("dx")
        buf.append("dy")
        buf.append("dz")
        buf.append("\n")
        buf.append(str(statistics))
        buf.append(str(init_particle.energy))
        buf.append("\n")

        file.write("\t".join(buf))

        for i in range(statistics):
            daughters = prop.propagate(init_particle)

            buf = [""]
            for d in daughters.stochastic_losses():
                buf.append(str(d.type))
                buf.append(str(d.propagated_distance))
                buf.append(str(d.energy))
                buf.append(str(d.position.x))
                buf.append(str(d.position.y))
                buf.append(str(d.position.z))
                buf.append(str(d.direction.x))
                buf.append(str(d.direction.y))
                buf.append(str(d.direction.z))
                buf.append("\n")

            file.write("\t".join(buf))


def main(dir_name):
    create_table(dir_name)


if __name__ == "__main__":

    dir_name = "TestFiles/"

    if os.path.isdir(dir_name):
        print("Directory {} already exists".format(dir_name))
    else:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))

    main(dir_name)
