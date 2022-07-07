#!/bin/bash -l
#module load python3/recommended

mpinum=16
qA=5
EA=20
EP=-10
RA=0.15
RP=${RA}
DENS=0.40
ETCLU=0.00

#for configfile in Input/Configurations/ConfigCluFrom_qA${qA}_dp0.50_dens${DENS}_eT${ETCLU}_nA10000_rp0.15_ra0.15_ep-10.0_ea20.0_T1.0_*/Config_C*.dat
for configfile in Input/Configurations/ConfigCluFrom_*/Config_C*.dat
do

	if [[ ! -e "$configfile" ]]
	then
		 continue
	fi

	echo $configfile

	for REAL in 100 #101 102 103 104 105 106 107 #108 109 110 111 112 113 114 115
	do

		filepattern=$(echo $configfile | sed 's/.*Configurations\/Config_//' | sed  's/\/Config.*//' | sed "s/_ts[0-9]*//")
		if ls Results*/TrajLog_${filepattern}_${REAL}.xyz 1> /dev/null 2>&1
		then
			echo "TrajLog with same filepattern already exists! File was skipped."
   			continue
		fi

		for EXTTORQUE in 0 #1 3 10 #0.1 0.3 1 3 5 10 15
		do
			python3 make_simulation_extractedcluster.py -ea ${EA} -ep ${EP} -ra ${RA} -rp ${RP} --real ${REAL} -eT ${EXTTORQUE} --MPInum ${mpinum} --MaxRunTime 240 --submitflag sbatch --configfile ${configfile}
			sleep 5

		done
	done
done
