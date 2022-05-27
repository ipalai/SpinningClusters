import numpy as np
import pandas as pd
import subprocess

# Look for all Traj files

dens=0.20
filestringList = subprocess.check_output("ls /Users/ivan/Dropbox/Mac/Documents/SpinningClusters/Simulations/Test_2/Results_qA5_dp0.50_dens{:.3f}_eT*/TrajLog_*_1[0-9][0-9].xyz | sed 's/.*Test_2\///' ".format(dens), shell=True)
#test filestringList = subprocess.check_output("ls /Users/ivan/Dropbox/Mac/Documents/SpinningClusters/Simulations/Test_2/Results_qA5_dp0.50_dens{:.3f}_eT30.00/TrajLog_*_114.xyz | sed 's/.*Test_2\///' ".format(dens), shell=True)
filestringList=filestringList.split()
filestringList=[a.decode("utf-8") for a in filestringList]
print(filestringList)

folder='/Users/ivan/Dropbox/Mac/Documents/SpinningClusters/Simulations/Test_2'

#filestringList=['Lc24_dens0.010_c_nC800_Rcyl100.0_cbP0.010_104']


from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *
from ovito.vis import *
import PySide6.QtCore
import os.path

nA = 10000

for filestring in filestringList:
    print(filestring)
    filepattern = subprocess.check_output(
        " echo {:s} | sed 's/.*\/TrajLog_//' | sed 's/.xyz//'".format(folder + '/' + filestring), shell=True)
    filepattern = filepattern.split()[0]
    filepattern = filepattern.decode("utf-8")
    ClusterStatisticsFile = '{}/Analysis/ClustersFiles/ClustersStatistics_{}.dat'.format(folder, filepattern)
    laststring = subprocess.check_output("tail -n5 {}/Results_*/Log_{}.dat".format(folder, filepattern), shell=True)

    # if (os.path.isfile(ClusterStatisticsFile)==False) & ("Total wall time" in laststring.decode("utf-8")):
    # if ("Total wall time" in laststring.decode("utf-8")):        # use this instead of previous line to overwrite file
    if True:  # rifugium peccatorum, needs to run again once sims are finished

        # Data import:
        pipeline = import_file('{}/{}'.format(folder, filestring), multiple_frames=True)

        # Select type:
        pipeline.modifiers.append(SelectTypeModifier(types={2, 3}))

        # Delete selected:
        pipeline.modifiers.append(DeleteSelectedModifier())

        # Cluster analysis:
        pipeline.modifiers.append(ClusterAnalysisModifier(
            cutoff=1.15,
            sort_by_size=True,
            unwrap_particles=True,  # needs to stay true, otherwise omega&AngMom per cluster will be wrong
            compute_com=True,
            compute_gyration=True,
            cluster_coloring=True))

        # Export Cluster analysis
        AvgClusterSize = []
        StdClusterSize = []
        for ThisFrame in range(pipeline.source.num_frames):
            data = pipeline.compute(frame=ThisFrame)
            print(data.attributes['Timestep'])
            if data.particles.count != nA:
                print("Missing atoms (only {:d} present) in timestep {:d}, file\n {:s}".format(data.particles.count, data.attributes['Timestep'], filestring))
                break
            #data.particles_.positions[:, 2] = 0
            ClustersListFile = '{}/Analysis/ClustersFiles/ClustersList_{}_ts{:d}.dat'.format(folder, filepattern, data.attributes['Timestep'])
            # export_file(data.tables['clusters'], '{}/Analysis/ClustersFiles/ClustersList_{}_ts{:d}.dat'.format(folder,filepattern,data.attributes['Timestep']),'txt/table', precision=3)
            ClustersOvito = data.tables['clusters']
            assert (ClustersOvito.x[:] == np.arange(1, len(ClustersOvito.x[:]) + 1, 1)).all(), "ClusterIDs are weird... Add additional check in the ClustersRotation computation (AngMom and Omega) "
            ClustersDataDf = pd.DataFrame(ClustersOvito.xy(), columns=['ClusterID', 'Size'])
            ClustersDataDf['Xcm'] = ClustersOvito['Center of Mass'][:, 0]
            ClustersDataDf['Ycm'] = ClustersOvito['Center of Mass'][:, 1]
            ClustersDataDf['Gxx'] = ClustersOvito['Gyration Tensor'][:, 0]  # the gyration tensor is computed as in lammps, for instance: Gxx =  1/M sum(mi xi^2), where xi is the position of particle i wrt the CM of its cluster
            ClustersDataDf['Gyy'] = ClustersOvito['Gyration Tensor'][:, 1]
            ClustersDataDf['Gxy'] = ClustersOvito['Gyration Tensor'][:, 3]
            ClustersDataDf['Iz'] = (ClustersOvito['Gyration Tensor'][:, 0] + ClustersOvito['Gyration Tensor'][:, 0]) * ClustersOvito['Cluster Size'][:]
            ClustersDataDf['Rg'] = np.sqrt(ClustersDataDf['Gxx'].values + ClustersDataDf['Gyy'].values)
            ClustersRotation = np.zeros((len(ClustersOvito.x[:]), 2))  # column 0 is ang mom, 1 is omega
            for i in range(0, data.particles.count):
                AtID = data.particles.identifiers[i]
                ClID = data.particles.cluster[i]
                rowID = ClID - 1
                ClData = ClustersDataDf[ClustersDataDf['ClusterID'] == ClID].iloc[0]
                RelativePos = data.particles.positions[i] - np.array([ClData['Xcm'], ClData['Ycm'], 0])  # z coordinate is meaningless here
                Velocity = data.particles.velocities[i]
                ClustersRotation[rowID, 0] += 1.0 * (RelativePos[0] * Velocity[1] - RelativePos[1] * Velocity[0])
                ClustersRotation[rowID, 1] += (RelativePos[0] * Velocity[1] - RelativePos[1] * Velocity[0]) / (RelativePos[0] ** 2 + RelativePos[1] ** 2 + 1e-100)
            for rowID in range(0, len(ClustersOvito.x[:])):
                ClID = rowID + 1
                ClSize = ClustersDataDf[ClustersDataDf['ClusterID'] == ClID].iloc[0].Size
                ClustersRotation[rowID, 1] /= ClSize
            ClustersDataDf['AngMomz'] = ClustersRotation[:, 0]
            ClustersDataDf['Omegaz'] = ClustersRotation[:, 1]
            ClustersDataDf['Omegaz'] = ClustersDataDf['Omegaz'].map(lambda x: '{:.3e}'.format(x))
            ClustersDataDf.to_csv(ClustersListFile, index=False, header=True, sep=' ', float_format='{:.3f}'.format)

            # data.attributes['AvgClusterSize']=data.tables['clusters'].xy()[:,1].mean()
            # data.attributes['StdClusterSize']=data.tables['clusters'].xy()[:,1].std()

            AvgClusterSize.append(data.tables['clusters'].xy()[:, 1].mean())
            StdClusterSize.append(data.tables['clusters'].xy()[:, 1].std())

        export_file(pipeline, ClusterStatisticsFile,
                    'txt/attr', precision=3, columns=["Timestep", "ClusterAnalysis.largest_size"], multiple_frames=True)

        with open(ClusterStatisticsFile, 'r') as txt:
            lines = txt.readlines()
        lines[0] = lines[0][:-1] + ' "AvgClusterSize" "StdClusterSize"\n'
        for i in range(1, len(lines)):
            lines[i] = lines[i][:-2] + ' {:.3f} {:+.3f}\n'.format(AvgClusterSize[i - 1], StdClusterSize[i - 1])
        with open(ClusterStatisticsFile, 'w') as f:
            f.writelines(lines)
        print("Written " + ClusterStatisticsFile)