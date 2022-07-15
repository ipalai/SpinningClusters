#! bin/bash

for folder in Results_*/
do
	echo $folder
	cd $folder
	for traj0 in TrajLog_qA5_dp0.50_*_20[1-9].xyz
	do
		if [ -f $traj0 ]
		then
			filepattern=$(echo $traj0 | sed 's/TrajLog_//')
			filemerged=$(echo $traj0 | sed 's/TrajLog_/TrajLog_merge_/')
			echo $filepattern
			echo $traj0
			if [ -f $filemerged ]
			then
				rm $filemerged
			fi
			cat $traj0 > $filemerged
			for i in {2..20}
			do	
				traji=TrajLog_${i}_${filepattern}
				if [ -f $traji ]
				then
					echo $traji
					cat $traji >> $filemerged
				fi
			done
			echo $filemerged
			echo "***"
		fi
	done
	cd ..
	echo "------"
done



