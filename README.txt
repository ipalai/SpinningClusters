Test_0 was used to develop the LAMMPS script that gived spinning and to test a bit the external torque parameter.

Test_1 contains tests about the shape of clusters and some sweep of parameter spaces (patch interaction, central interaction, ranges...)
)

Test_2 contains simulations at density=0.2 or 0.1, nA~10000, RunSteps~5e7, to get some more statistics and check the fractal dimension. Dump of trajectories is both linear and logarithmic, though linear dump is not downloaded.

Test_3 contains tests for oriented binding. Central atoms form bonds beetween themselves. When a bond is formed an angle is also created, so that the bond is directional. This allows to disentangle the action of rotation from that of folding. The role of patches is the same as before, i.e. just create friction among colloids and transmit torque. 
