#!/bin/bash -l

cd ..
executable=ClusterAnalysisFluidljcossq1.2.exe

for ts in 50000000 #200000 400000 500000 600000 800000 1000000 2000000 5000000 10000000 20000000 50000000 # 1000000 5000000 10000000 50000000 
do
for q in 5 
do
for DIR in Results_el*_qA${q}_qB${q}
do

echo $ts $q $DIR
cp Analysis/$executable $DIR/
cd $DIR

for FILE in Movie_2d_ljcossq_ep10.0_sp2.0_ea*_ra2.0_el*_sl0.1_rl0.2_d0.475_nA200_qA${q}_nB*_qB${q}_T1_Tdamp1.0_dA0.0300_3*.xyz
do
ClusterAnalysisFile=$(echo $FILE | sed 's/Movie_2d_/ClusterAnalysis_/' | sed 's/.xyz/_ts'${ts}'.dat/' )
#echo $ClusterAnalysisFile

if [ ! -f ../Analysis/ClusterAnalysisFiles/${ClusterAnalysisFile} ]; then
echo '******************'
./$executable $FILE $ts $ClusterAnalysisFile
else
echo $ClusterAnalysisFile ' already exists'
fi
echo ' '

done
cd ..
done
done
done
mv Results_*/ClusterAnalysis_*.dat Analysis/ClusterAnalysisFiles
