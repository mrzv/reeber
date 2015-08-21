#!/bin/bash

fn=$1

../mt-lg-ghosts-double $fn -b 8 -n LyAReeberTest_plt00242_128x128x128_density-b8-n.lg
diff LyAReeberTest_plt00242_128x128x128_density-b8-n.pi <(../persistent-integral-lg-double -v -x 200 -i 81.66 LyAReeberTest_plt00242_128x128x128_density-b8-n.lg LyAReeberTest_plt00242_128x128x128_density-pi && cat LyAReeberTest_plt00242_128x128x128_density-pi-b* | ./sort.sh) || exit 1

../mt-lg-ghosts-double $fn -b 8 -n -w LyAReeberTest_plt00242_128x128x128_density-b8-n-w.lg
diff LyAReeberTest_plt00242_128x128x128_density-b8-n-w.pi <(../persistent-integral-lg-double -v -x 200 -i 81.66 LyAReeberTest_plt00242_128x128x128_density-b8-n-w.lg LyAReeberTest_plt00242_128x128x128_density-pi-w && cat LyAReeberTest_plt00242_128x128x128_density-pi-w-b* | ./sort.sh) || exit 1
