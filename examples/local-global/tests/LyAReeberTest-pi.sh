#!/bin/bash

fn=$1

../mt-lg-ghosts-double $fn -b 8 -n LyAReeberTest_plt00242_128x128x128_density-b8-n.lg
../persistent-integral-lg-double -v -x 200 -i 81.66 LyAReeberTest_plt00242_128x128x128_density-b8-n.lg LyAReeberTest_plt00242_128x128x128_density-pi
diff LyAReeberTest_plt00242_128x128x128_density-b8-n.pi <(./sort.sh LyAReeberTest_plt00242_128x128x128_density-pi-b*) || exit 1

../mt-lg-ghosts-double $fn -b 8 -n -w LyAReeberTest_plt00242_128x128x128_density-b8-n-w.lg
../persistent-integral-lg-double -v -x 200 -i 81.66 LyAReeberTest_plt00242_128x128x128_density-b8-n-w.lg LyAReeberTest_plt00242_128x128x128_density-pi-w
diff LyAReeberTest_plt00242_128x128x128_density-b8-n-w.pi <(./sort.sh LyAReeberTest_plt00242_128x128x128_density-pi-w-b*) || exit 1
