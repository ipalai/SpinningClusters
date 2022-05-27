#!/bin/bash -l
#module load python3/recommended

nA=10000
mpinum=9
DENS=0.1


for REAL in 100 #101 102 103 104 105 106 107 #108 109 110 111 112 113 114 115
do
for EXTTORQUE in 0 0.1 0.3 1 3 10 30 100
do
for qA in 5 #6 8
do
for EA in 20 #5 6 8 10 #6 7 8 #7 8 9 10 1
do
for EP in -10 #-5 5 10 20 #0.5 1.0 1.5 2.0 2.5 3.0 #1 1.5 2.5 3  #2 3
do
for RA in 0.15 #0.5 1.0 2.0 #1 1.5 2
do
RP=${RA}
python3 make_simulation.py -nA ${nA} -ea ${EA} -ep ${EP} -ra ${RA} -rp ${RP} -dens ${DENS} --real ${REAL} -qA ${qA} -eT ${EXTTORQUE} --MPInum ${mpinum} --RunSteps 5000000 --dumpevery 50000 --MaxRunTime 200
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
