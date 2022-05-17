#!/bin/bash -l

rsync -av --progress --include={'*.xyz','*.dat','*.txt','*.sh'} --include='*/' --exclude='*' ipalaia@marcin80.ista.local:/nfs/scistore15/saricgrp/ipalaia/SpinningClusters/Simulations/SpinningClusters/* ./
#rsync -av --progress -e "ssh -A -J ucapipa@socrates.ucl.ac.uk:22" ucapipa@myriad.rc.ucl.ac.uk:/lustre/home/ucapipa/Scratch/TheoryClustersNowca/* ./

