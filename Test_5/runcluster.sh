#!/bin/bash -l
#module load python3/recommended

nA=100000
mpinum=25

for DENS in 0.1 0.2
do
for REAL in 100 #101 102 103 104 105 106 107 #108 109 110 111 112 113 114 115
do
for EXTTORQUE in 0 0.1 0.3 1 3 5 10 15 30 
do
for qA in 5
do
for EA in 20
do
for EP in -10
do
for RA in 0.15
do
RP=${RA}
python3 make_simulation.py -nA ${nA} -ea ${EA} -ep ${EP} -ra ${RA} -rp ${RP} -dens ${DENS} --real ${REAL} -qA ${qA} -eT ${EXTTORQUE} --MPInum ${mpinum} --RunSteps 3000000 --dumpevery 50000 --MaxRunTime 200 --clusterevery 100
runscriptfile=$(ls -rt runscript_* | tail -n 1)
echo $runscriptfile
sbatch $runscriptfile
mv $runscriptfile runscriptFiles/ 
sleep 10
done
done
done
done
done
done
done
