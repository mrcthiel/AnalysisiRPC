#!/bin/bash

scanID=$1
HV=$2

mkdir plots_904
mkdir outputs
mkdir plots_gif
mkdir plots_gif/DQM_plots
cp index.php plots_904/
cp index.php plots_gif/
cp index.php plots_gif/DQM_plots/

root -l -b <<EOF
.L run_ana_FEBv2.C
run_ana_FEBv2(${HV},${scanID},"./data/")
.q
EOF

#root -l -b <<EOF
#.L dqm_beamdump.C
#dqm_beamdump(${HV},${scanID},"/eos/cms/store/user/mthiel/roottrees/")
#.q
#EOF

mkdir ScanId_${scanID}/
mkdir ScanId_${scanID}/HV${HV}/
cp index.php ScanId_${scanID}/
cp index.php ScanId_${scanID}/HV${HV}/
mv plots_904 ScanId_${scanID}/HV${HV}/
mv plots_gif/DQM_plots ScanId_${scanID}/HV${HV}/

