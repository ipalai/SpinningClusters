#!/bin/bash -l
#module load python3/recommended

mpinum=1
RunSteps=300000
qA=5
EA=20
EP=-10
RA=0.15
RP=${RA}

for configfile in Input/Configurations/ConfigCluFrom_*/Config_C{1..99}.dat
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
   			continue
		fi

		for EXTTORQUE in 0.1 0.3 1 3 5 10 15
		do

			python3 make_simulation_extractedcluster.py -ea ${EA} -ep ${EP} -ra ${RA} -rp ${RP} --real ${REAL} -eT ${EXTTORQUE} --MPInum ${mpinum} --RunSteps ${RunSteps} --dumpevery 50000 --MaxRunTime 24 --submitflag run --configfile ${configfile}
			sleep 10

		done
	done
done
