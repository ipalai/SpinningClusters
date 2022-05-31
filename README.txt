Test_0 was used to develop the LAMMPS script that gived spinning and to test a bit the external torque parameter.

Test_1 contains tests about the shape of clusters and some sweep of parameter spaces (patch interaction, central interaction, ranges...)
)

Test_2 contains simulations at density=0.2 or 0.1, nA~10000, RunSteps~5e7, to get some more statistics and check the fractal dimension. Dump of trajectories is both linear and logarithmic, though linear dump is not downloaded.

Test_3 contains tests for oriented binding. I tried two strategies and a mix of them:
1) Central atoms form bonds beetween themselves. When a bond is formed an angle is also created, so that the bond is directional. This requires very strong angle and bond energy to oppose rotation when clusters are big. The timestep was decreased from 0.008 to 0.002, and still rotation partially happens. 
2) Clusters are treated as rigid bodies: namely, every "clusterevery" timesteps (e.g. 100 ts), a compute cluster/atom is performed and the fix rigid/nve updated with the new rigid bodies. This seem to work very ewell.

Test_4 contains simulation using the latter mechanism (on-the-go clustering). This allows to suppress folding, thus isolating the effect of pure rotation. However, clustering does not work properly: clusters across boundaries do not rotate. Maybe a bug in LAMMPS, or a problem with images of particles inside clusters.

Test_5 contains simulations as before (rigid clustering on the go), but with fixed boundary conditions and walls. This avoids the problem with clusters across boundaries not moving. Clustering cutoff reduced to sigma+InteractionRange*2/3=0.10 (it was 0.15).

Test_6 contains simulations where I isolate preformed clusters from simulations in Test_5, and allow for folding.
