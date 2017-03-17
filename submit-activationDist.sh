#!/bin/bash
#FREQ LASTCOOR PROBFILE PYMDIR NSTRUCT LSTART
#arg1: target probability
#arg2: last pym file prefix
#arg3: probability file
#arg4: path to pymfile
#arg5: number of structures
#arg6: start reading line
#arg7: bead index file
#arg8: blocksize
#arg9: output directory

python_path="/home/cmb-panasas2/same/Software/anaconda2/bin"

die () {
 echo >&2 "$@"
 exit 1
}

maindir=`pwd`
model='domain1169'
runName=".DamID_HiC"
inputDir=${maindir}"/input"

prog="$maindir/DM.peffective.nextActivationDistance.incr.py"
lastpym=$1 #name like "p010.l040" the hi-c & lamin prefixes
nextfreq=$2 #next pxxx for hi-c
nextladp=$3 #next lxxx for lad, not used for freq here, just a name
ns=$4 # the population size
pfile=$inputDir/${model}.from40k.prob 
beadfile=$inputDir/${model}.ind
pymdir=$maindir
outputdir=$maindir

[ $# -gt 3 ] || die "4 argument required, $# provided: lastpym nextfreq nextladp i.e.:  p006i.l040 p006j l020 10000"

np=`wc -l $pfile |awk '{print $1}'`  # 650000
echo $np

i=1
blocksize=2500 #submit a job for a block
for l in `seq 0 $blocksize $np`
do
   echo $i $l
   pbsfile="$i.peffectiveActDist.pbs"
#   sleep 1
cat << ENDHERE > $pbsfile
#!/bin/bash
#PBS -q cmb
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -l mem=1777mb
#PBS -l pmem=1777mb
#PBS -l vmem=1777mb
#PBS -l walltime=10:00:00
#PBS -m a
source /usr/usc/python/default/setup.sh

cd \$PBS_O_WORKDIR

${python_path}/python $prog $nextfreq $lastpym $pfile $pymdir $ns $l $beadfile $blocksize $outputdir $runName $nextladp

ENDHERE

   qsub $pbsfile
   let i=i+1
done
