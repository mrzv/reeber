#!/usr/bin/env python3

qname = "debug"
cores_per_node = 24

def write_sl(file_name, z, data_size, n_cores, max_minutes):
    n_nodes = 1 + n_cores // 24
    rho = "81.66"
    theta = "90.0"

    npy_fname = "z{0}-{1}-f4.npy".format(z, data_size)
    dgm_fname = "$SCRATCH/output/my_answers/amr_cc_{0}_blocks_{1}_rho_{2}_diagram.txt".format(npy_fname, n_cores, rho)
    intl_fname = "$SCRATCH/output/my_answers/amr_cc_{0}_blocks_{1}_rho_{2}_theta_{3}_integral.txt".format(npy_fname, n_cores, rho, theta)
    npy_full_fname = "$SCRATCH/striped_dataset/{0}".format(npy_fname)
    mt_fname = "$SCRATCH/output/my_answers/amr_cc_z_{0}_size_{1}_blocks_{2}.mt".format(z, data_size, n_cores)
    log_fname = "log_cc_size_{0}_blocks_{1}_z_{2}_npy.txt".format(data_size, n_cores, z)

    script_template = """#!/bin/bash -l

#SBATCH -q debug
#SBATCH -N {0}
#SBATCH -t 00:{1}:00
#SBATCH -L SCRATCH
#SBATCH -o {2}
#SBATCH -e {2}

srun -n {3} ~/code/reeber2-git/build_edison_release/examples/amr-connected-components/amr_connected_components_float -a -i {4} -n -w -x {5} {6} {7} {8} {9}"""

    with open(file_name, 'w') as sl_file:
        sl_file.write(script_template.format(n_nodes, max_minutes, log_fname, n_cores, rho, theta, npy_full_fname, mt_fname, dgm_fname, intl_fname))


if __name__ == "__main__":
    max_minutes = 15
    sbatch_file_name = "queue_all.sh"
    with open(sbatch_file_name, 'w') as sbatch_file:
        sbatch_file.write('#!/bin/sh\n\n')
        # for z in ['2', '6']:
        for z in ['2']:
            for data_size in [512, 2048]:
                for n_cores in [512, 1024, 2048, 4096, 8192]:
                    file_name = "run_cc_{0}_{1}_{2}_npy.sl".format(data_size, n_cores, z)
                    write_sl(file_name, z, data_size, n_cores, max_minutes)
                    sbatch_file.write("sbatch {0}\n".format(file_name))
