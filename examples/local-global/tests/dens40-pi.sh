#!/bin/bash

b=8; ../mt-lg-ghosts-double dens40.npy -b $b -n dens40-b$b-n.lg
b=8; ../persistent-integral-lg-double -a -v -x 7.5e9 -i 7e9 dens40-b$b-n.lg dens40-pi-b$b
b=8; diff dens40-b8-n.pi <(./sort.sh dens40-pi-b$b-b*) || exit 1

b=8; ../mt-lg-ghosts-double dens40.npy -b $b -n -w dens40-b$b-n-w.lg
b=8; ../persistent-integral-lg-double -a -v -x 7.5e9 -i 7e9 dens40-b$b-n-w.lg dens40-pi-w-b$b
b=8; diff dens40-b8-n-w.pi <(./sort.sh dens40-pi-w-b$b-b*) || exit 1

b=10; ../mt-lg-ghosts-double dens40.npy -b $b -n dens40-b$b-n.lg
b=10; ../persistent-integral-lg-double -a -v -x 7.5e9 -i 7e9 dens40-b$b-n.lg dens40-pi-b$b
b=10; diff dens40-b8-n.pi <(./sort.sh dens40-pi-b$b-b*) || exit 1

# No wrap with 10 since 10 has only two prime factors, and wrap doesn't work with a single division on a side

b=20; ../mt-lg-ghosts-double dens40.npy -b $b -n dens40-b$b-n.lg
b=20; ../persistent-integral-lg-double -a -v -x 7.5e9 -i 7e9 dens40-b$b-n.lg dens40-pi-b$b
b=20; diff dens40-b8-n.pi <(./sort.sh dens40-pi-b$b-b*) || exit 1

b=20; ../mt-lg-ghosts-double dens40.npy -b $b -n -w dens40-b$b-n-w.lg
b=20; ../persistent-integral-lg-double -a -v -x 7.5e9 -i 7e9 dens40-b$b-n-w.lg dens40-pi-w-b$b
b=20; diff dens40-b8-n-w.pi <(./sort.sh dens40-pi-w-b$b-b*) || exit 1

b=32; ../mt-lg-ghosts-double dens40.npy -b $b -n dens40-b$b-n.lg
b=32; ../persistent-integral-lg-double -a -v -x 7.5e9 -i 7e9 dens40-b$b-n.lg dens40-pi-b$b
b=32; diff dens40-b8-n.pi <(./sort.sh dens40-pi-b$b-b*) || exit 1

b=32; ../mt-lg-ghosts-double dens40.npy -b $b -n -w dens40-b$b-n-w.lg
b=32; ../persistent-integral-lg-double -a -v -x 7.5e9 -i 7e9 dens40-b$b-n-w.lg dens40-pi-w-b$b
b=32; diff dens40-b8-n-w.pi <(./sort.sh dens40-pi-w-b$b-b*) || exit 1

# TMT

b=8; ../tmt-lg-ghosts-double dens40.npy -b $b -n dens40-b$b-n.tmt     || exit 1
b=8; ../triplet-persistence-lg-double dens40-b$b-n.tmt dens40-tmt-pd-n-b$b
b=8; diff dens40-tmt-pd-n.dgm <(./sort.sh dens40-tmt-pd-n-b$b-b*) || exit 1

b=16; ../tmt-lg-ghosts-double dens40.npy -b $b -n dens40-b$b-n.tmt     || exit 1
b=16; ../triplet-persistence-lg-double dens40-b$b-n.tmt dens40-tmt-pd-n-b$b
b=16; diff dens40-tmt-pd-n.dgm <(./sort.sh dens40-tmt-pd-n-b$b-b*) || exit 1

# -w

b=64; ../tmt-lg-ghosts-double dens40.npy -b $b -n -w dens40-b$b-n-w.tmt     || exit 1
b=64; ../triplet-persistence-lg-double dens40-b$b-n-w.tmt dens40-tmt-pd-n-w-b$b
b=64; diff dens40-tmt-pd-n-w.dgm <(./sort.sh dens40-tmt-pd-n-w-b$b-b*) || exit 1

b=128; ../tmt-lg-ghosts-double dens40.npy -b $b -n -w dens40-b$b-n-w.tmt     || exit 1
b=128; ../triplet-persistence-lg-double dens40-b$b-n-w.tmt dens40-tmt-pd-n-w-b$b
b=128; diff dens40-tmt-pd-n-w.dgm <(./sort.sh dens40-tmt-pd-n-w-b$b-b*) || exit 1
