#step 0: setup the parameters & prepare the folders
ns=10
let ns_less=ns-1
runName='DamID_HiC'

mkdir ensemble
cd ensemble
for i in {0..$ns_less};
do
    mkdir copy$i;
done

#step 1: inital structures with all the constraints with prob=1
currentp='p100'
lastp='p100'
lad_p='l100'
lad_lastp='l100'
./submit-modeling.sh $currentp $lastp $lad_p $lad_lastp $ns_less

#step 2: Activation distance
currentp='p100' #'p006s'
nextp='p070' #'p006t'
lad_currentp='l100' #'l020'
lad_nextp='l020' #'l010'
./submit-activationDist.sh $currentp.$lad_currentp $nextp $lad_nextp $ns
cat ${nextp}.${lad_nextp}.${runName}.peffect_ActDist.* > ${nextp}.${lad_nextp}.${runName}.peffect_ActDist
rm -f ${nextp}.${lad_nextp}.${runName}.peffect_ActDist.*
./submit-ladActivationDist.sh $currentp.$lad_currentp $nextp $lad_nextp $ns

#step 3: Modeling
currentp='p070' #'p006t'
lastp='p100' #'p006s'
lad_p='l020' #'l010'
lad_lastp='l100' #'l020'
./submit-modeling.sh $currentp $lastp $lad_p $lad_lastp $ns_less

###repeat step 2 and step 3 until all the restraints are satisified within certain tolerance 
