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

This script uses multiprocessing. You can specify the number of
process as fist argument of this script.

"""


import multiprocessing
import sys
import warnings

import bremsstrahlung
import continous_randomization
import epairproduction
import ionization
import photonuclear
import propagation
import scattering
import sector


def main():
    print("There are %d CPUs on this machine" % multiprocessing.cpu_count())

    number_processes = 2

    try:
        number_processes = int(sys.argv[1])
    except IndexError:
        pass
    except ValueError:
        warnings.warn("first argument must be an integer. (Number of processes to ues)")
        pass

    pool = multiprocessing.Pool(number_processes)

    pool.apply_async(bremsstrahlung.main)
    pool.apply_async(continous_randomization.main)
    pool.apply_async(epairproduction.main)
    pool.apply_async(ionization.main)
    pool.apply_async(photonuclear.main)
    pool.apply_async(propagation.main)
    pool.apply_async(scattering.main)
    pool.apply_async(sector.main)

    pool.close()
    pool.join()

if __name__ == "__main__":
    main()
