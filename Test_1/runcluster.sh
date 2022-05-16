#!/bin/bash -l
module load python3/recommended

nA=200
mpinum=1

for REAL in 311 312 313 314 315 #301 302 303 304 305 #1 2 3 #51 52 
do
for qA in 5
do
for EA in 5 #5 6 8 10 #6 7 8 #7 8 9 10 1
do
for EP in 1.5 #0.5 1.0 1.5 2.0 2.5 3.0 #1 1.5 2.5 3  #2 3
do
for RA in 2.0 #0.5 1.0 2.0 #1 1.5 2
do
python3 make_simulation.py -nA ${nA} -ea ${EA} -ep ${EP} -ra ${RA} --density_A_target 0.03 --real ${REAL} -qA ${qA} -qB ${qB} --mpinum ${mpinum} --runsteps 50000000 --ligand_range 0.2 --dumpstep 1000000
runscriptfile=$(ls -rt runscript_* | tail -n 1)
echo $runscriptfile
qsub $runscriptfile
mv $runscriptfile runscriptFiles/ 
sleep 1
done
done
done
done
done
done
done
