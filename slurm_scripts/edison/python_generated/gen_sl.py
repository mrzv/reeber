#!/usr/bin/env python3

qname = "regular"
cores_per_node = 24

def write_sl(file_name, z, data_size, n_cores, max_minutes):
    n_nodes = 1 + n_cores // 24

    script_template = """#!/bin/bash -l

    #SBATCH -q debug
    #SBATCH -N {0}
    #SBATCH -t 00:{1}:00
    #SBATCH -L SCRATCH
    #SBATCH -o log_simple_size_{2}_blocks_{3}_z_{4}_npy.txt
    #SBATCH -e log_simple_size_{2}_blocks_{3}_z_{4}_npy.txt

    srun -n {3} ~/code/reeber2-git/build_edison/examples/amr-merge-tree/amr_merge_tree_simple_float -t 81.66 -n $SCRATCH/dataset/z{4}-{2}-f4.npy none"""

    with open(file_name, 'w') as sl_file:
        sl_file.write(script_template.format(n_nodes, max_minutes, data_size, n_cores, z))


if __name__ == "__main__":
    max_minutes = 25
    sbatch_file_name = "queue_all.sh"
    with open(sbatch_file_name, 'w') as sbatch_file:
        sbatch_file.write('#!/bin/sh\n\n')
        for z in ['2', '6']:
            for data_size in [512, 2048]:
                for n_cores in [512, 1024, 2048, 4096, 8192]:
                    file_name = "run_simple_{0}_{1}_{2}_npy.sl".format(data_size, n_cores, z)
                    write_sl(file_name, z, data_size, n_cores, max_minutes)
                    sbatch_file.write("sbatch {0}\n".format(file_name))
