#!/usr/bin/env python3

import os.path
import subprocess as sp
import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))

build_dir = "./../../../../build_master"

execname_si = f"{build_dir}/examples/amr-merge-tree/amr_merge_tree_simple_float"
execname_sc = f"{build_dir}/examples/amr-connected-components/amr_connected_components_float"
execname_mt = f"{build_dir}/examples/local-global/mt-lg-ghosts-float"
execname_pi = f"{build_dir}/examples/local-global/persistent-integral-lg-float"

nblocks = "64"

executables_found = (os.path.exists(execname_si) and
                     os.path.exists(execname_sc) and
                     os.path.exists(execname_mt) and
                     os.path.exists(execname_pi))

if not executables_found:
    raise FileNotFoundError("One of the executables is missing. Rebuild the project or edit paths in script")


beta_params = [(0.5, 0.5), (5, 1), (1, 3), (2,2), (2,5)]

fnames = []

def generate_data():
    global fnames
    fnames = []
    np.random.seed(1)
    fname = "a-64-64-64-normal-float.npy"
    fnames.append(fname)
    if not os.path.exists(fname):
        a = np.random.randn(64,64,64)
        np.save(fname, a)
    for (alpha, beta) in beta_params:
        fname = "alpha-{}-beta-{}-size-64-float.npy".format(alpha, beta)
        fnames.append(fname)
        if not os.path.exists(fname):
            beta = np.random.beta(alpha, beta, size=(64,64,64))
            np.save(fname, beta)


def save_correct_tree(input_npy_fname, tree_fname, log_file):
    if not os.path.exists(tree_fname):
        sp.call([execname_mt, "-b", nblocks, "-w", "-n", input_npy_fname, tree_fname], stdout=log_file, stderr=log_file)


def save_correct_integral(input_npy_fname, out_int_fname, rho, log_file):
    tree_fname = current_dir + "/" + input_npy_fname[:-4] + ".tree"
    save_correct_tree(input_npy_fname, tree_fname, log_file)
    if not os.path.exists(out_int_fname):
        exec_str = f"./get_intl.sh {tree_fname} {out_int_fname} {rho} {rho} {execname_pi}"
        sp.call(exec_str, shell = True)


def read_integral_file(fname):
    integral = {}
    if os.path.exists(fname):
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

n_failed = 0
n_passed = 0
rhos = ["0.2", "0.3", "0.5", "0.6", "0.9"]
theta = "1"

generate_data()

for fname in fnames:
    for rho in rhos:

        out_tree_simple_fname = "none"
        out_diagram_simple_fname = "none"
        out_integral_simple_fname = f"{current_dir}/{fname}-b-{nblocks}-simple-integral-rho-{rho}-theta-{theta}.txt"

        out_tree_sc_fname = "none"
        out_diagram_sc_fname = "none"
        out_integral_sc_fname = f"{current_dir}/{fname}-b-{nblocks}-sc-integral-rho-{rho}-theta-{theta}.txt"

        out_integral_correct_fname = f"{current_dir}/{fname}-correct-integral-rho-{rho}-theta-{theta}.txt"

        with open("compare-integrals-log.txt", 'a') as log_file:

            save_correct_integral(fname, out_integral_correct_fname, rho, log_file)

            retcode = sp.call([execname_si, "-b", nblocks, "-w", "-n", "-a", "-i", rho, "-x",
                theta, "-f", "density", fname, out_tree_simple_fname, out_diagram_simple_fname, out_integral_simple_fname],
                stdout=log_file, stderr=log_file)
            if retcode != 0:
                print("Execution failed: {} {} {}-{}".format(execname_si, fname, rho, theta))
                n_failed += 1

            # retcode = sp.call([execname_sc, "-b", nblocks, "-w", "-n", "-a", "-i", rho, "-x", theta, fname, 
            #     out_tree_sc_fname, out_diagram_sc_fname, out_integral_sc_fname], 
            #     stdout=log_file, stderr=log_file)
            # if retcode != 0:
            #     print("Execution failed: {} {} {}-{}".format(execname_sc, fname, rho, theta))
            #     n_failed += 1
            # else:
            #     sp.call(["./assemble_intl.sh", out_integral_sc_fname])


        if not compare_integrals(out_integral_simple_fname, out_integral_correct_fname):
            print("Test failed: {} {} {}-{}".format("simple", fname, rho, theta))
            n_failed += 1
        else:
            n_passed += 1

        # if not compare_integrals(out_integral_sc_fname, out_integral_correct_fname):
        #     print("Test failed: {} {} {}-{}".format("send_components", fname, rho, theta))
        #     n_failed += 1
        # else:
        #     n_passed += 1

        print("Processed {}".format(fname))

if n_failed == 0:
    print("All {} tests passed".format(n_passed))
else:
    print("Failed tests: {}".format(n_failed))
