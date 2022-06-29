#!/bin/bash -l

rsync -av --progress --human-readable --include={'TrajLog_*.xyz','Out_*.dat','Log_*.dat','RotationStatsMolecule*.dat','*.sh','Input*.in'} --include='*/' --exclude='*' ipalaia@marcin80.ista.local:/nfs/scistore15/saricgrp/ipalaia/SpinningClusters/Simulations/SpinningClusters/* ./
#rsync -av --progress -e "ssh -A -J ucapipa@socrates.ucl.ac.uk:22" ucapipa@myriad.rc.ucl.ac.uk:/lustre/home/ucapipa/Scratch/TheoryClustersNowca/* ./

