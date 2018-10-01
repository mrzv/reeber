#!/bin/bash

if [ "$#" -ne 5 ]; then
    echo "Incorrect number of parameters"
    echo "$#"
    exit 1
fi

infname=$1
full_integral_name=$2
#full_integral_name="${infname%.tree}.npy-correct-integral-rho-${rho}-theta-${theta}.txt"
rho=$3
theta=$4
execname=$5

"$execname" -a -i "$rho" -x "$theta" -l fatal "$infname" "${infname}-intl" 2>&1
cat ${infname}-intl-b* > ${full_integral_name}
rm ${infname}-intl-b*
sort -o "${full_integral_name}" "${full_integral_name}"
