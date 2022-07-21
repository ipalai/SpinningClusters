#!/bin/bash -l
#module load python3/recommended

mpinum=1
EA=20
EP=-10
RA=0.15
RP=${RA}

DENS=0.40


for configfile in Input/Configurations/ConfigCluFrom_qA5_dp0.50_dens${DENS}_eT*_nA10000_rp0.15_ra0.15_ep-10.0_ea20.0_T1.0_201_ts10000000/Config_C*.dat
#for configfile in Input/Configurations/ConfigCluFrom_*/Config_C*.dat
do

	if [[ ! -e "$configfile" ]]
	then
		 continue
	fi

	echo $configfile

	for REAL in 95 #101 102 103 104 105 106 107 #108 109 110 111 112 113 114 115
	do

		filepattern=$(echo $configfile | sed 's/.*Configurations\/Config_//' | sed  's/\/Config.*//' | sed "s/_ts[0-9]*//")
		if ls Results*/Traj*_${filepattern}_${REAL}.xyz 1> /dev/null 2>&1
		then
			echo "Traj with same filepattern already exists! File was skipped."
   			continue
		fi

		EXTTORQUE=0
		for DEFDIR in x
		do
			python3 make_simulation_extractedcluster.py -ea ${EA} -ep ${EP} -ra ${RA} -rp ${RP} --real ${REAL} -eT ${EXTTORQUE} --MPInum ${mpinum} --MaxRunTime 240 --submitflag no --configfile ${configfile} --DeformDirection ${DEFDIR} --DeformDelta 1.0 --dumpevery 100
			sleep 5

		done
	done
done
