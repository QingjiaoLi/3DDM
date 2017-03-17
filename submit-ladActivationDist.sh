#!/bin/bash
#calculate the activation distance for NE interacting beads
python_path="/home/cmb-panasas2/same/Software/anaconda2/bin"

die () {
 echo >&2 "$@"
 exit 1
}


maindir=`pwd`
model='domain1169'
runName=".DamID_HiC"
inputDir=${maindir}"/input"

prog="$maindir/LAD.peffective.nextActivationDistance.py"
lastpym=$1  #this is dual: pxxx.lxxx, the previous file pym to read
nextfreq=$2 #this is the hi-c pxxx
nextladp=$3 #this is the lad lxxx
ns=$4 # the population size
pfile=$inputDir/${model}.from40k.prob 
beadfile=$inputDir/${model}.ind
pymdir=$maindir
outputdir=$maindir

[ $# -gt 3 ] || die "4 argument required, $# provided: lastpym nextfreq nextladp i.e.:  p006i.l040 p006j l020 10000"

${python_path}/python $prog $nextfreq $lastpym $pymdir $ns $beadfile $outputdir $runName $nextladp

