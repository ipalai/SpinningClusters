#!bin/bash/

for q in 5 #4
do
for EL in 7 #5 #6 7 8 9 10 12 14 20
do
for EA in 1.0 #1.5 #0.5 1.0 1.5 2.0 2.5 3.0 3.5
do
for RA in 2.0 #2.0 #0.5 1.0 1.5 2.0 2.5 3.0
do
for nB in 50 100 150 200 300 500 800 1200 1600
do

qA=$q
qB=$q

tmpfile=tmpfileClusterTimeEvol.dat
findfile=(ClusterAnalysisStatisticsFiles/ClusterAnalysisStatistics_ljcossq_ep10.0_sp2.0_ea${EA}_ra${RA}_el${EL}*_nA200_qA${qA}_nB${nB}_qB${qB}_*_dA0.03*_ts*.dat)

if [ -s ${findfile[0]} ]
then

	OutputFile=$(echo ${findfile} | awk '{print $1}' | sed 's/ClusterAnalysisStatistics_/ClusterTimeEvol_/g' | sed 's/_ts.*/.dat/g')
	if [ -s $OutputFile ]
	then
		rm $OutputFile
	fi

	echo "# 0) ts  1) NumberofClustersBigger +StDev  3) NumberOfAClusters +StDev  5) AvgClustersSizeBigger +StDev  7) AvgAConnectivity +StDev  9) AvgBConnectivity +StDev  11) CoalescenceLikelihood1 +StDev  13) CoalescenceLikelihood2 +StDev  15) FractionOfAInClustersOfSizeUpTo1 +StDev  17) FractionOfAInClustersOfSizeUpTo5 +StDev  19) FractionOfAInClustersOfSizeUpTo10 +StDev  21) FractionOfFreeB +StDev 23) AvgAClusterSizeBigger +StDev 25) NumberOfTerminalNodes +StDev 27) NumberOfUnboundNodes +StDev 29) Compactness +StDev 31) Roundness +StDev" >> $tmpfile
	for file in ${findfile[@]}  #$(find ./ -maxdepth 1 -name '${findfile}')
	do
        if [ -s ${file} ]
        then
		ts=$(echo $file | sed 's/C.*_ts//' | sed 's/.dat//')
		Arg1=$(grep -A1 '^NumberOfClustersBigger' $file | tail -n 1)
		Arg3=$(grep -A1 '^NumberOfAClusters' $file | tail -n 1)
		Arg5=$(grep -A1 '^AvgClusterSizeBigger' $file | tail -n 1)
		Arg7=$(grep -A1 '^AvgAConnectivity' $file | tail -n 1)
		Arg9=$(grep -A1 '^AvgBConnectivity' $file | tail -n 1)
		CoalescenceLikelihoodOne=$(grep -A1 '^CoalescenceLikelihood1' $file | tail -n 1)
		CoalescenceLikelihoodTwo=$(grep -A1 '^CoalescenceLikelihood2' $file | tail -n 1)
		FractionOfAInClustersOfSizeUpTo1=$(grep -A1 '^FractionOfAInClustersOfSizeUpTo1 ' $file | tail -n 1)
		FractionOfAInClustersOfSizeUpTo5=$(grep -A1 '^FractionOfAInClustersOfSizeUpTo5 ' $file | tail -n 1)
		FractionOfAInClustersOfSizeUpTo10=$(grep -A1 '^FractionOfAInClustersOfSizeUpTo10 ' $file | tail -n 1)
		FractionOfFreeB=$(grep -A1 '^FractionOfFreeB ' $file | tail -n 1)
		Arg23=$(grep -A1 '^AvgAClusterSizeBigger ' $file | tail -n 1)
		Arg25=$(grep -A1 '^NumberOfTerminalNodes ' $file | tail -n 1)
		Arg27=$(grep -A1 '^NumberOfUnboundNodes ' $file | tail -n 1)		
		Arg29=$(grep -A1 '^Compactness' $file | tail -n 1)		
		Arg31=$(grep -A1 '^Roundness' $file | tail -n 1)

		echo $ts $Arg1 $Arg3 $Arg5 $Arg7 $Arg9 $CoalescenceLikelihoodOne $CoalescenceLikelihoodTwo $FractionOfAInClustersOfSizeUpTo1 $FractionOfAInClustersOfSizeUpTo5 $FractionOfAInClustersOfSizeUpTo10 $FractionOfFreeB $Arg23 $Arg25 $Arg27 $Arg29 $Arg31 '#'$file >> $tmpfile
		echo $file
		echo $OutputFile
	fi
	done
	
	(head -n 1 $tmpfile && tail -n +2 $tmpfile | sort -n) > $OutputFile
	rm $tmpfile

fi
done
done
done
done
done

mv ClusterAnalysisStatisticsFiles/ClusterTimeEvol_* ClusterTimeEvolFiles/
