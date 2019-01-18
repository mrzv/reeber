#!/usr/bin/env python3

import os.path
import subprocess as sp

# General comment regarding format: cannot use full-scale string interpolation
# because Nersc does not have the latest python3.
# Otherwise clumsy format(name1=name1, name2=name2, ..) wouldn't be necessary.

def read_diagram_file(fname):
    dgm = set()
    with open(fname, 'r') as f:
        for line in f:
            b, d = line.split()
            b = float(b)
            d = float(d)
            dgm.add((b,d))
    return dgm


def read_integral_file(fname):
    integral = {}
    with open(fname, 'r') as f:
        for line in f:
            x,y,z,v = line.split()
            x = int(x)
            y = int(y)
            z = int(z)
            v = float(v)
            integral[(x,y,z)] = v
    return integral


def compare_integrals(fname1, fname2):
    integral1 = read_integral_file(fname1)
    integral2 = read_integral_file(fname2)
    if (len(integral1) != len(integral2)):
        return False
    for vertex in integral1:
        v1 = integral1[vertex]
        try:
            v2 = integral2[vertex]
            if (v1 - v2) / v1  > 0.01:
                return False
        except KeyError:
            return False
    return True


def compare_diagrams(fname1, fname2):
    dgm1 = read_diagram_file(fname1)
    dgm2 = read_diagram_file(fname2)
    return dgm1 == dgm2


if __name__ == "__main__":
    rho = "81.66"
    theta = "90.0"

    n_failed = 0
    n_passed = 0

    simple_dir = "/scratch2/scratchdirs/greynarn/output/my_answer/final/simple"
    cc_dir = "/scratch2/scratchdirs/greynarn/output/my_answer/final/cc"

    # amr_cc_z2-1024-f4.npy_blocks_1024_rho_81.66_diagram.txt
    # amr_mt_simple_z2-1024-f4.npy_blocks_1024_rho_81.66_diagram.txt
    # amr_cc_z2-1024-f4.npy_blocks_2048_rho_81.66_diagram.txt


    for data_size in ["512", "1024", "2048"]:
        for nblocks in ["512", "1024", "2048", "4096", "8192"]:
            fname = "z2-{}-f4.npy".format(data_size)
            integral_simple_fname = "{simple_dir}/amr_mt_simple_{fname}_blocks_{nblocks}_rho_{rho}_theta_{theta}_integral.txt".format(nblocks=nblocks,
                    simple_dir=simple_dir,fname=fname, rho=rho, theta=theta)
            integral_cc_fname = "{cc_dir}/amr_cc_{fname}_blocks_{nblocks}_rho_{rho}_theta_{theta}_integral.txt".format(nblocks=nblocks,
                    cc_dir=cc_dir,fname=fname, rho=rho, theta=theta)
            if not compare_integrals(integral_simple_fname, integral_cc_fname):
                print("Test failed: {} {} {}-{}".format(integral_simple_fname, integral_cc_fname, rho, theta))
                n_failed += 1
                continue
            else:
                print("Integrals match: {} {} {}-{}".format(integral_simple_fname, integral_cc_fname, rho, theta))

            diagram_simple_fname = "{simple_dir}/amr_mt_simple_{fname}_blocks_{nblocks}_rho_{rho}_diagram.txt".format(nblocks=nblocks,
                    simple_dir=simple_dir, fname=fname, rho=rho)
            diagram_cc_fname = "{cc_dir}/amr_cc_{fname}_blocks_{nblocks}_rho_{rho}_diagram.txt".format(nblocks=nblocks,
                    cc_dir=cc_dir, fname=fname, rho=rho)

            if not compare_diagrams(diagram_simple_fname, diagram_cc_fname):
                print("Test failed: {} {}, rho = {}".format(diagram_simple_fname, diagram_cc_fname, rho))
                n_failed += 1
            else:
                print("Diagrams match: {} {}, rho = {}".format(diagram_simple_fname, diagram_cc_fname, rho))
                n_passed += 1


    if n_failed == 0:
        print("All {} tests passed".format(n_passed))
    else:
        print("Failed tests: {}".format(n_failed))
