import pyPROPOSAL as pp


def create_table(dir_name):
    """TODO: Docstring for create_table.
    Returns: TODO

    """

    statistics = 10

    prop = pp.Propagator(
        pp.particle.MuMinusDef.get(),
        "resources/config_ice.json"
    )

    mu = prop.particle

    mu.energy = 1e8
    mu.propagated_distance = 0
    mu.position = pp.Vector3D(0, 0, 0)
    mu.direction = pp.Vector3D(0, 0, -1)
    pp.RandomGenerator.get().set_seed(1234)

    with open(dir_name + "Propagator_propagation.txt", "a") as file:

        buf = [""]
        buf.append("name")
        buf.append("lenght")
        buf.append("energy")
        buf.append("x")
        buf.append("y")
        buf.append("z")
        buf.append("dx")
        buf.append("dy")
        buf.append("dz")
        buf.append("\n")
        buf.append(str(statistics))
        buf.append(str(mu.energy))
        buf.append("\n")

        file.write("\t".join(buf))

        for i in range(statistics):
            mu.energy = 1e8
            mu.propagated_distance = 0
            mu.position = pp.Vector3D(0, 0, 0)
            mu.direction = pp.Vector3D(0, 0, -1)

            daughters = prop.propagate()

            buf = [""]
            for d in daughters:
                if d.id == pp.particle.Data.Particle:
                    buf.append(d.particle_def.name)
                else:
                    buf.append(str(d.id).split(".")[1])
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

    import os

    dir_name = "TestFiles/"

    if os.path.isdir(dir_name):
        print("Directory {} already exists".format(dir_name))
    else:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))

    main(dir_name)
