Realisations 100, 101,... 115 are the core simulations used for most of the paper for density 0.20.

Realisations 200, 201, ... are much longer simulations, run at various densities to build the phase diagram (Fig. 4 of the draft). 200 has no printed Restart files, 201 and 202 do. The restarted simulations (ran through runcluster_restart.sh or restart_simulation.py) have file name patterns that look like
 Traj_2_*, Traj_3_*, Traj_4_*, ...
for, respectively, the 2nd restart, the 3rd, the 4th...

The script mergeTrajLog.sh is used to merge files from consecutive restarts of the same simulation, for realisations 201, 202, ... . Notice that this is a brute merge, so in the merged file frames are not necessary in chronological order as simulations may restart from sightly before the last printed time step.