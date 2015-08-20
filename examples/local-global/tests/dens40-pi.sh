#!/bin/bash

../mt-lg-ghosts-double dens40.npy -b 8 -n dens40-b8-n.lg
diff dens40-b8-n.pi <(../persistent-integral-lg-double -a -v -x 7.5e9 -i 7e9 dens40-b8-n.lg dens40-pi && cat dens40-pi-b* | sort) || exit 1

../mt-lg-ghosts-double dens40.npy -b 8 -n -w dens40-b8-n-w.lg
diff dens40-b8-n-w.pi <(../persistent-integral-lg-double -a -v -x 7.5e9 -i 7e9 dens40-b8-n-w.lg dens40-pi-w && cat dens40-pi-w-b* | sort) || exit 1
