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

import multiprocessing
import sys
import os
import warnings
import subprocess

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
    dir_name = "TestFiles/"
    tar_name = "TestFiles.tar.gz"

    try:
        number_processes = int(sys.argv[1])
    except IndexError:
        pass
    except ValueError:
        warnings.warn("first argument must be an integer. (Number of processes to ues)")
        pass

    try:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))
    except OSError:
        print("Directory {} already exists".format(dir_name))

    pool = multiprocessing.Pool(number_processes)
    results = []

    pool.apply_async(bremsstrahlung.main, (dir_name, ))
    pool.apply_async(continous_randomization.main, (dir_name, ))
    pool.apply_async(epairproduction.main, (dir_name, ))
    pool.apply_async(ionization.main, (dir_name, ))
    pool.apply_async(photonuclear.main, (dir_name, ))
    pool.apply_async(propagation.main, (dir_name, ))
    pool.apply_async(scattering.main, (dir_name, ))
    pool.apply_async(sector.main, (dir_name, ))

    pool.close()
    pool.join()

    print("all threads are joined")

    p = subprocess.Popen(['tar', '-czf', tar_name, dir_name])

    p.communicate()
    print("compressed test files {}".format(tar_name))

    p = subprocess.Popen(['rm', '-r', dir_name])

    p.communicate()
    print("Directory {} removed".format(dir_name))

if __name__ == "__main__":
    main()
