#!/bin/bash

scanID=$1
HV=$2

mkdir ScanId_${scanID}/
mkdir ScanId_${scanID}/HV${HV}/
cp index.php ScanId_${scanID}/
cp index.php ScanId_${scanID}/HV${HV}/

root -l -b <<EOF
.L run_ana_FEBv2.C
run_ana_FEBv2(${HV},${scanID},"/eos/cms/store/user/mthiel/roottrees/")
.q
EOF

#root -l -b <<EOF
#.L dqm_beamdump.C
#dqm_beamdump(${HV},${scanID},5000,"/eos/cms/store/user/mthiel/roottrees/")
#.q
#EOF

