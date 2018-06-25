#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Creating test files for the unit tests of PROPOSAL

This internal package is thought to be used by developers and
maintainers of PROPOSAL only.

You can execute every module separate to create test files
for the specific tests or run __init__.py to execute all modules.
Therfore you must be in the tests/gen_testfiles directory, because
propagation.py uses the config_ice.json located in the
resources directory which is hard coded in that module.

"""

import bremsstrahlung
import continous_randomization
import epairproduction
import ionization
import photonuclear
import propagation
import scattering
import sector


def main():
    bremsstrahlung.main()
    continous_randomization.main()
    epair_production.main()
    ionization.main()
    photonuclear.main()
    propagation.main()
    scattering.main()
    sector.main()

if __name__ == "__main__":
    main()
