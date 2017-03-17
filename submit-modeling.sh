#!/bin/bash
python_path="/home/cmb-panasas2/same/Software/anaconda2/bin"

die () {
 echo >&2 "$@"
 exit 1
}

[ $# -gt 3 ] || die "4 argument required, $# provided: current_probability  last_probability lad_p lad_lastp, i.e: p006j p006i l020 l040"
maindir=`pwd`
prog='DM.simult_DamID_HiC'

CODEDIR=$maindir
modelName='domain1169'
inputDir=${maindir}"/input"
countname='fly'
currentfb=$1 
lastfb=$2
lad_p=$3
lad_lastp=$4
ns=$5

ensemble=$maindir/'ensemble'

function job {
    jobnum=$3
    copyStart=$1
    copyEnd=$2
    pbs="$jobnum.$countname.pbs"
cat << ENDHERE > $pbs
#!/bin/bash
#PBS -q cmb
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -l mem=950mb
#PBS -l pmem=990mb
#PBS -l vmem=990mb
#PBS -l walltime=20:00:00
#PBS -m a

cd $ensemble

LD_LIBRARY_PATH='/home/cmb-04/fa/htjong/local/imp-1.0/build/lib:/home/cmb-04/fa/htjong/myBin/lib/:/home/cmb-04/xjz/qingjial/tools/localpython/lib'
export LD_LIBRARY_PATH

PYTHONPATH='/home/cmb-04/fa/htjong/local/imp-1.0/build/lib:/home/cmb-04/xjz/qingjial/tools/localpython/bin'
export PYTHONPATH

IMP_BUILD_ROOT='/home/cmb-04/fa/htjong/local/imp-1.0'
export IMP_BUILD_ROOT

for i in {$copyStart..$copyEnd}
do
    newdir="copy\$i"
    pymdir="$ensemble/\$newdir"
    if [ ! -d \$newdir ]; then
        mkdir \$newdir
    fi
    cd \$newdir
        ${python_path}/python $CODEDIR/$prog.py $maindir $inputDir $modelName $currentfb $lastfb \$pymdir $lad_p $lad_lastp > ${prog}_${currentfb}_${lad_p}.log
    cd ..

done


ENDHERE

qsub $pbs
sleep 1
}

#--------------

i=1 #pbs index
struct=1 # number of structures per job
for aw in `seq 0 $struct $ns`
do
    let ak=aw+struct-1
    job $aw $ak $i
    let i=i+1
done

