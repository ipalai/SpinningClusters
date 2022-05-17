#!bin/bash/

cd ClusterAnalysisFiles/
executable=ClusterAnalysisStatisticsFluidljcossq1.1.exe
logfile=log_ClusterAnalysisStatistics_Fluidljcossq.dat
cp ../${executable} ./
echo ${logfile}
date > ${logfile} 
echo ${executable} >> ${logfile}

for ts in 50000000 #50000000 #200000 400000 500000 600000 800000 1000000 2000000 5000000 10000000 20000000 50000000 #1000000 5000000 10000000 50000000
do
for EL in 5.0 #4.0 4.5 5.0 5.5 6 7 8 9 10 
do
for EA in 1.5 #0.5 1.0 1.5 2.0 2.5 3.0
do
for RA in 2.0 #0.5 1.0 1.5 2.0 #2.5 3.0
do
for q in 5
do
for nB in 50 100 150 200 300 400 500 800 1200 1600
do	

qA=$q
qB=$q

echo '' >> ${logfile}
findfile=(ClusterAnalysis_ljcossq_ep10.0_sp2.0_ea${EA}_ra${RA}_el${EL}*_sl0.1_rl0.2_d0.475_nA200_qA${qA}_nB${nB}_qB${qB}_T1_Tdamp1.0_dA0.030*_ts${ts}.dat)

if [ -s ${findfile[0]} ]
then
	echo $findfile
	OutputFile=$(echo ${findfile} | awk '{print $1}' | sed 's/ClusterAnalysis_/ClusterAnalysisStatistics_/g' | sed 's/_\?[0-9]\?[0-9]\?[0-9]_ts/_ts/')
	./${executable} ${findfile[@]} ${OutputFile} >> ${logfile}
fi
done
done
done
done
done
done

mv ClusterAnalysisStatistics_* ../ClusterAnalysisStatisticsFiles/
mv ${logfile} ../
rm ${executable}
