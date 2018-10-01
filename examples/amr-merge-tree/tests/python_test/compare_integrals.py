#!/usr/bin/env python3

import os.path
import subprocess as sp
import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))

execname_si = "./../../../../build/examples/amr-merge-tree/amr_merge_tree_simple_double"
execname_sc = "./../../../../build/examples/amr-merge-tree/amr_merge_tree_send_components_double"
execname_mt = "./../../../../build/examples/local-global/mt-lg-ghosts-double"
execname_pi = "./../../../../build/examples/local-global/persistent-integral-lg-double"
nblocks = "8"


def generate_data():
    a = np.random.randn(64,64,64)
    np.save("a-64-64-64-normal-double", a)


def write_correct_answer(input_npy_fname, out_int_fname, rho, theta, log_file):
    tree_fname = current_dir + "/" + input_npy_fname[:-4] + ".tree"
    sp.call([execname_mt, "-b", nblocks, "-w", "-n", input_npy_fname, tree_fname], stdout=log_file, stderr=log_file)
    exec_str = "./get_intl.sh {tree_fname} {out_int_fname} {rho} {theta} {execname_pi}".format(tree_fname=tree_fname,
            out_int_fname=out_int_fname, rho=rho, theta=theta, execname_pi=execname_pi)
    sp.call(exec_str, shell = True)
    sp.call(["./get_intl.sh",  tree_fname, out_int_fname, rho, theta, execname_pi], stdout=log_file, stderr=log_file)


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


fnames = ["a-64-64-64-normal-double.npy"]
bounds = [("0.5", "0.8"), ("0.2", "0.4"), ("0.5", "0.7")]
# bounds = [("0.5", "0.8")]
for fname in fnames:
    if not os.path.exists(fname):
        generate_data()
    for (rho, theta) in bounds:
        out_tree_simple_fname = "none"
        out_diagram_simple_fname = "{current_dir}/{fname}-simple-rho_{rho}-diagram.txt".format(current_dir=current_dir,fname=fname, rho=rho)
        out_integral_simple_fname = "{current_dir}/{fname}-simple-integral-rho-{rho}-theta-{theta}.txt".format(current_dir=current_dir,fname=fname, rho=rho, theta=theta)
        out_integral_correct_fname = "{current_dir}/{fname}-correct-integral-rho-{rho}-theta-{theta}.txt".format(current_dir=current_dir,fname=fname, rho=rho, theta=theta)
        with open("compare-integrals-log.txt", 'a') as log_file:
            write_correct_answer(fname, out_integral_correct_fname, rho, theta, log_file)
        with open("compare-integrals-log.txt", 'a') as log_file:
            retcode = sp.call([execname_si, "-b", nblocks, "-w", "-n", "-a", "-i", rho, "-x", theta, fname, out_tree_simple_fname, out_diagram_simple_fname, out_integral_simple_fname], stdout=log_file, stderr=log_file)
            # retcode = sp.call([execname_si, "-b", nblocks, "-w", "-n", "-a", "-i", rho, "-x", theta, fname, out_tree_simple_fname, out_diagram_simple_fname, out_integral_simple_fname])
        result = compare_integrals(out_integral_simple_fname, out_integral_correct_fname)
        print("{}: {}".format(os.path.basename(out_integral_simple_fname), result))
        assert(result)

