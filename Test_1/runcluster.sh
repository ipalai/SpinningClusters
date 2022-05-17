#!/bin/bash -l
#module load python3/recommended

nA=200
mpinum=1
DENS=0.2
RP=0.3
EXTTORQUE=1.0

for REAL in 100 
do
for qA in 5 6 8
do
for EA in 20 #5 6 8 10 #6 7 8 #7 8 9 10 1
do
for EP in -20 -10 -5 5 10 20 #0.5 1.0 1.5 2.0 2.5 3.0 #1 1.5 2.5 3  #2 3
do
for RA in 0.05 0.2 0.5 1 #0.5 1.0 2.0 #1 1.5 2
do
python3 make_simulation.py -nA ${nA} -ea ${EA} -ep ${EP} -ra ${RA} -rp ${RP} -dens ${DENS} --real ${REAL} -qA ${qA} -eT ${EXTTORQUE} --MPInum ${mpinum} --RunSteps 50000000 --dumpevery 100000
runscriptfile=$(ls -rt runscript_* | tail -n 1)
echo $runscriptfile
sbatch $runscriptfile
mv $runscriptfile runscriptFiles/ 
sleep 1
done
done
done
done
done
