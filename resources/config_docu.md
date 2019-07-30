
# Documentation of the Configuration file #

In the configuration file you can specify all parameters for your propagation.
It's written in a .json file. That's why we have this (nice) Docu. ;)

The configuration file is structured into three parts:
- In the "global" part the general parameters are defined like the path for the interpolation tables.
- In the "sectors" part the different media with their geometry are defined. The energy-cuts, which means how accurate the propagation process should be calculated, can be set for each sector, which then overrides the global configurations.
- In the "detector" part the geometry of the detector is defined.

## The `global` configurations ##
There is the option to set a seed, if the internal random number generator is used.
The seed has no effect, if you use an external random number generator.

If the Output of the Secondaries should not only include the Stochastic energy losses (and the Particles produced in a decay), but also the continuous energy losses, this can be set.

When the Output should just contain the secondaries (energy losses or decay products), that accured inside the detector volume, and not the ones outside of the detector, this can also be set.

| Keyword                     | Type    | Default   | Description |
| --------------------------- | ------- | --------- | ----------- |
| `seed`                      | Integer | `0`       | seed for the internal random number generator|
| `continous_loss_output`     | Bool    | `False`   | Decides, whether continuous losses should be emitted in the Output of Secondaries|
| `only_loss_inside_detector` | Bool    | `False`   | Decides, whether the secondaries outside of the detector should be included in the output, or not.|

### Interpolation parameters ###
The `interpolation` parameter is an own json-object.
This object can contain multiple parameters dealing with the interpolation tables.

There are two kind of paths, that can be set; a path, where the program only requires reading permission and a path, where it requires writing permission.
The argument can be a String, or a list of Strings.
If the given String or the first String in the list fulfills the readable/writable requirement, it is set.If the first item in the list don't satisfy the requirement, the program is looking for the next and so on until it reaches the end of the list.
If none of the given elements in the list or the given String is a valid path, then the program writes the tables in the memory.

If both paths are valid, PROPOSAL first looks at the reading path.
If the interpolation table file doesn't exist or is empty (e.g. if two instances are called in parallel and one process is already building the table but hasn't finished yet), PROPOSAL uses writing path.
There it again looks, if the interpolation table file already exist or is not empty.
If the tables have been already build PROPOSAL just uses them.
If there are no tables (or at least not ones with the desired particle or medium properties) PROPOSAL builds the tables in the folder.
If the String is empty, the Folder doesn't exists or PROPOSAL has no permission, the tables are stored in the cache.
Note: The tables differ in the parameters given below, that are stored in the file name. For not too long file names, these values are hashed.

There is the option that just the readonly path should be used. So if there is not the required tables prebuild in the readonly path the Initialization/program wil break and not try to look or write at the `path_to_tables` or in the memory.

The tables can be stored as a binary file or as a text file.

It's also possible to do the calculations with integrations, but that increases the amount of time by around 4 orders of magnitude !!! This should just be used for tests and comparisons. The Interpolations are accurate enough. Note: The reason why integration is that much slower is that the interpolation used for propagation is near the surface, using already calculated numbers from 'lower' interpolations. When integrating, these integrals going 'deep' (up to 4 layers).

The upper energy limit can be modified up to the maximum possible primary particle energy, 
to prevent energies greater than the maximum energy from being extrapolated.
If particles are propagated whose primary energy is greater than the energy of the maximum node, the error increases rapidly. 
This should be avoided.

If the error of the interpolation becomes too large, the number of sampling points can be increased. 
This increases the runtime of the program.

| Keyword                   | Type   | Default | Description |
| ------------------------- | ------ | ------- | ----------- |
| `do_interpolation`              | Bool   | `True`  | Decides, whether to calculate with interpolations or integrations |
| `path_to_tables`                | String | `""`    | path pointing to the folder with the interpolation tables |
| `path_to_tables_readonly`       | String | `""`    | path pointing to the folder with the interpolation tables |
| `do_binary_tables`              | Bool   | `True`  | Decides, whether the tables are stored in binary format or in human readable text format |
| `just_use_readonly_path`        | Bool   | `False` | Decides, if just the readonly path should be used. |
| `max_node_energy`               | Double | `1.e14` | Energy in MeV up to which the interpolation tables are built |
| `nodes_cross_section`           | Integer| `100`   | Number of interpolation points for the interpolation of the cross section integral |
| `nodes_continous_randomization` | Integer| `200`   | Number of interpolation points for the interpolation of the continous randomization integral |
| `nodes_propagate`               | Integer| `1000`  | Number of interpolation points for the interpolation of the propagation integral |

### Accuracy parameters ###
There are several parameters with which the precision or speed for advancing the particles can be adjusted.

For the time calculation the energy dependence can be included via the dE/dx.
Or the time is calculated by simply dividing distance by the speed of light.

There are four different multiple scattering parametrizations describing the deviation of the propagation direction.
  - `"Moliere"` Parametrization of [Moliere](http://zfn.mpdl.mpg.de/data/Reihe_A/3/ZNA-1948-3a-0078.pdf)
  - `"Highland"` First Order approximation (a gaussian function) of Moliere's theory derived by [Highland](https://doi.org/10.1016/0029-554X(75)90743-0) and corrected by [Lynch/Dahl](https://doi.org/10.1016/0168-583X(91)95671-Y).
  - `"HighlandIntegral"` From the old PROPOSAL version (which had just this scattering mode). Also using Highland approximation (corrected by Lynch/Dahl), but taking into account the energy dependence of the propagated distance
  - `"NoScattering"` Here the particle always propagate in the initial direction without deviation from the propagation axes.

If a particle gets below a lower threshold energy (this is a property of a particle that should be defined in the Particle Definitions and is by default the same as the particle mass) but hasn't decayed, it can get forced to stop and to decay.

| Keyword          | Type   | Default     | Description |
| ---------------- | ------ | ----------- | ----------- |
| `exact_time`     | Bool   | `True`      | Decide, whether the energy dependence of the time is included |
| `scattering`     | String | `"Moliere"` | multiple scattering parametrization describing the displacement from the initial propagation direction |
| `stopping_decay` | Bool   | `True`      | Decide, whether the particle gets stopped and forced to decay, if its energy is too low. |

### Cross section parameters ###
There are several parametrizations defining the cross sections and further options which simply scales them.

The cross section multiplier, available for each cross section, scales this cross section by its factor.

For **Ionization**, there is just one parametrization, but for pair production, bremsstrahlung and nuclear interaction, it's possible to choose between multiple parametrizations.

The **electron pair production** parametrizations are:
  - `"EpairKelnerKokoulinPetrukhin"` (Proc. 12th ICCR (1971), 2436) with corrections for the interaction with atomic electrons (Phys. Atom. Nucl. 61 (1998), 448)
  - `"EpairSandrockSoedingreksoRhode"` 

The **bremsstrahlung** parametrizations are:
  - `"BremsKelnerKokoulinPetrukhin"` ([Preprint MEPhI (1995) no. 024-95](http://cds.cern.ch/record/288828)) and (Phys. Atom. Nucl. 62 (1999), 272)
  - `"BremsAndreevBezrukovBugaev"` (Phys. Atom. Nucl. 57 (1994), 2066)
  - `"BremsPetrukhinShestakov"` (Canad. J. Phys. 46 (1968), 377)
  - `"BremsCompleteScreening"` taken from Tsai [Rev. Mod. Phys. 46 (1974), 815](https://doi.org/10.1103/RevModPhys.46.815)
  - `"BremsSandrockSoedingreksoRhode"`

There are two different approaches to parametrise the **nuclear interaction**:

Either with the approximation where a real photon scatters inelastically with a nucleus with an additional factor to make the photon virtual.
The available parametrizations using this approach are:
- `"PhotoKokoulin"`
- `"PhotoRhode"`
- `"PhotoBezrukovBugaev"` (Sov. J. Nucl. Phys. 33 (1981), 635) with corrections for particles with higher mass like taus ([Phys. Rev. D67 (2003), 034027](https://doi.org/10.1103/PhysRevD.67.034027))
- `"PhotoZeus"` ([Eur. Phys. J. C 7 (1999), 609](https://doi.org/10.1007/s100529901084))

For these parametrizations the hard component can additionally be considered (only affecting the parametrizations using the real photon approximation).

Or using data from electron proton scattering experiments and extrapolate to low momentum of the virtual photon transfered to the nucleus ($Q^2$).
There the cross section is first integrated over the photon momentum.
The available parametrizations using this approach are:
- `"PhotoAbramowiczLevinLevyMaor91"` by Abramowicz Levin Levy Maor [Phys. Let. B269 (1991), 465](https://doi.org/10.1016/0370-2693(91)90202-2)
- `"PhotoAbramowiczLevinLevyMaor97"` by Abramowicz Levin Levy Maor [arXiv::hep-ph/9712415](https://arxiv.org/abs/hep-ph/9712415)
- `"PhotoButkevichMikhailov"` [JETP 95 (2002), 11](https://doi.org/10.1134/1.1499897)
- `"PhotoRenoSarcevicSu"` [Astrop. Phys. 24 (2005), 107](https://doi.org/10.1016/j.astropartphys.2005.06.002) (which uses the parametrization of ALLM97, but with corrections for spin 0 particles and should therefore only be selected for spin 0 particles like sTau)

For these parametrizations the parametrization of the shadowing factor can be chosen from one of the following parametrizations (only affecting the nuclear interaction parametrizations with momentum integration):
- `"ShadowButkevichMikhailov"` from their calculation of nuclear interaction
- `"ShadowDuttaRenoSarcevicSeckel"` by Dutta, Reno, Sarcevic, Seckel [Phys. Rev. D63 (2001), 094020](https://doi.org/10.1103/PhysRevD.63.094020)

The **LPM effect** (Landau-Pomeranschuk-Migdal), suppressing the bremsstrahlung and the pair production at high energies and the Ter-Mikaelian effect, suppress low bremsstrahlung energy losses, can also be incorporated.

The **muon pair production** (which is an optional process and per default disabled) parametrizations are:
  - `"MupairKelnerKokoulinPetrukhin"` (Phys. Atom. Nucl. Vol. 63, No.9 (2000),  pp. 1603-1611, DOI: 10.1134/1.1312894)

| Keyword                 | Type   | Default    | Description |
| ----------------------  | ------ | ---------- | ----------- |
| `brems_multiplier`      | Double | `1.0`        | scales the bremsstrahlung |
| `epair_multiplier`      | Double | `1.0`        | scales the electron pair production |
| `ioniz_multiplier`      | Double | `1.0`        | scales the ionization |
| `photo_multiplier`      | Double | `1.0`        | scales the nuclear interaction |
| `mupair_multiplier`     | Double | `1.0`        | scales the muon pair production |
| `epair`                 | String | `"EpairKelnerKokoulinPetrukhin"` | electron pair production parametrization |
| `epair`                 | String | `"MupairKelnerKokoulinPetrukhin"` | muon pair production parametrization |
| `brems`                 | String | `"BremsKelnerKokoulinPetrukhin"` | Bremsstrahlung parametrization |
| `photo`                 | String | `"PhotoAbramowiczLevinLevyMaor97"` | nuclear interaction parametrization |
| `photo_hard_component`  | Bool   | `True`     | including the hard components |
| `photo_shadow`          | String | `"ShadowButkevichMikhailov"` | shadowing parametrization |
| `lpm`                   | Bool   | `True`     | Incorporate the LPM-effect and TM-effect |
| `mupair_enable`         | Bool   | `False`     | Include production of muon pairs in the calculation of energy losses |
| `mupair_particle_output`| Bool   | `True`     | Produced muon pairs are treated as Particles with corresponding energies in the Output of Secondaries (and not as DynamicData objects) |


### Energy-cut parameters ###
The energy cut settings and the continuous randomization option are separated between `cuts_inside`, `cuts_infront` and `cuts_behind` the detector, which are again own json-objects.
**ecut** describes the cut in the total energy loss in MeV and **vcut** describes the cut in the energy loss relative to the particle energy.
Above this cut, the energy losses are calculated stochastically and below this cut, the energy loss is calculated via the average energy loss (continuously).
- If both values are set, the minimum of the total energy cut and the relative cut times the particle energy is applied.
- If just one of the energy cuts should be used, the other value should be set to a negative value (usually it's -1), because negative values for the energy cuts are neglected.
- If both values are set to a negative value, only continuous losses (average energy loss dE/dx) are taking into account.

The continuous randomization randomizes the continuous losses, which affects, that the continuous losses are not always the same for the same particle energy, because they are randomised a little bit.
This is just a useful for small energy cuts.

Note: The energy cuts and the continuous randomization settings can also be specified for each sector.
Then the global settings will be overwritten.

For the `cuts_inside` option, the default values are

| Keyword     | Type   | Default   | Description |
| ----------- | ------ | --------- | ----------- |
| `e_cut`     | Double | `500.0`   | total energy loss cut inside the detector |
| `v_cut`     | Double | `-1.0`    | relative energy loss cut inside the detector |
| `cont_rand` | Bool   | `True`    | includes the continuous randomization inside the detector |

For the `cuts_infront` option, the default values are

| Keyword     | Type   | Default   | Description |
| ----------- | ------ | --------- | ----------- |
| `e_cut`     | Double | `-1.0`    | total energy loss cut in front the detector |
| `v_cut`     | Double | `0.001`   | relative energy loss cut in front the detector |
| `cont_rand` | Bool   | `True`    | includes the continuous randomization in front the detector |

For the `cuts_behind` option, the default values are

| Keyword     | Type   | Default   | Description |
| ----------- | ------ | --------- | ----------- |
| `e_cut`     | Double | `-1.0`    | total energy loss cut behind the detector |
| `v_cut`     | Double | `-1.0`    | relative energy loss cut behind the detector |
| `cont_rand` | Bool   | `False`   | includes the continuous randomization behind the detector |


## The `sectors` configurations ##

Every sector needs a medium and a geometry, while it's priority, compared to other sectors, can be scaled with the hierarchy option.
If multiple sectors are defined and if these sectors overlap, the decision which sector should be preferred is made with the hierarchy option. The bigger hierarchy is preferred. If the overlapping sectors have the same hierarchy, then the one with the more dense medium is chosen. If both hierarchy and Massdensity is equal, then the first in the sector list is chosen.

| Keyword     | Type    | Default | Description |
| ----------- | ------- | ------- | ----------- |
| `hierarchy` | Integer | `1`     | decides, which sector is favored, when they overlap |

### Medium parameters ###
The following media are implemented:
- `"Hydrogen"`
- `"Iron"`
- `"Copper"`
- `"Lead"`
- `"Uranium"`
- `"Air"`
- `"Water"`
- `"Ice"`
- `"StrandardRock"`
- `"FrejusRock"`
- `"Salt"`
- `"Paraffin"`
- `"AntaresWater"`

A density correction factor for the medium, specifying the density of the medium for the experiment, can optionally be set.

| Keyword              | Type   | Default | Description |
| -------------------- | ------ | ------- | ----------- |
| `medium`             | String | `"Water"` | Medium in this sector |
| `density_correction` | Double | `1.0`   | density correction factor of the medium |

### Geometry parameters ###

Every Geometry needs an origin as a vector and a shape.

The following three shapes are available:
- `"Sphere"`
- `"Cylinder"`
- `"Box"`

| Keyword  | Type                     | Default           | Description |
| -------- | ------------------------ | ----------------- | ----------- |
| `shape`  | String                   |  `-`              | Shape of the Geometry |
| `origin` | [Double, Double, Double] | `[0.0, 0.0 ,0.0]` | Center of the Geometry |

Every shape has its own specifications:
- a sphere needs an inner and an outer radius between which the sector is defined
- a cylinder needs also an inner and an outer radius between which the sector is defined and a height (there is just one height, not an inner and an outer one)
- a box needs a length, a width and a height. Despite the other two geometries, where the sector is defined between two shapes, the box has just one shape and the sector is inside this.

| Keyword        | Type   | Default | Description |
| -------------- | ------ | ------- | ----------- |
| `outer_radius` | Double | `-`     | outer radius of the sphere or the cylinder |
| `inner_radius` | Double | `-`     | inner radius of the sphere or the cylinder |
| `height`       | Double | `-`     | height of the cylinder or the box (length in z-direction) |
| `length`       | Double | `-`     | length of the box (length in x-direction) |
| `width`        | Double | `-`     | width of th box (length in y-direction) |

### Energy cut parameters ###

The energy cut parameters can be specified for every sector, which then overwrite the globally defined cut settings.
The cuts uses the same specifications like the global ones and for a description see the previous section.

## The `detector` configurations ##

The Detector is defined by a Geometry, which parameters can be specified like the geometry of a sector described in the previous section.
