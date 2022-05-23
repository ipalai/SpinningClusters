
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt


# Look for all Traj files
import subprocess

from ovito.plugins.ParticlesPython import ClusterAnalysisModifier
from ovito.plugins.StdModPython import SelectTypeModifier, DeleteSelectedModifier

filestringList = subprocess.check_output("ls /Users/ivan/Dropbox/Mac/Documents/SpinningClusters/Simulations/Test_2/Results_*/TrajLog_*_1[0-9][0-9].xyz | sed 's/.*Test_2\///' ", shell=True)
filestringList=filestringList.split()
filestringList=[a.decode("utf-8") for a in filestringList]
print(filestringList)


from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *
from ovito.vis import *
import PySide6.QtCore
import os.path

for filestring in filestringList:
    
    folder='/Users/ivan/Dropbox/Mac/Documents/SpinningClusters/Simulations/Test_2'
    filepattern = subprocess.check_output(" echo {:s} | sed 's/.*\/TrajLog_//' | sed 's/.xyz//'".format(folder+'/'+filestring), shell=True)
    filepattern=filepattern.split()[0]
    filepattern = filepattern.decode("utf-8")
    ClusterStatisticsFile='{}/Analysis/ClustersFiles/ClustersStatistics_{}.dat'.format(folder,filepattern)
    laststring = subprocess.check_output("tail -n5 {}/Results_*/Log_{}.dat".format(folder,filepattern), shell=True)
    
    #if (os.path.isfile(ClusterStatisticsFile)==False) & ("Total wall time" in laststring.decode("utf-8")):        
    #if ("Total wall time" in laststring.decode("utf-8")):        # use this instead of previous line to overwrite file
    if True:                                                      # rifugium peccatorum, needs to run again once sims are finished
        
        # Data import:
        pipeline = import_file('{}/{}'.format(folder,filestring), multiple_frames = True)

        # Select type:
        pipeline.modifiers.append(SelectTypeModifier(types = {2, 3}))

        # Delete selected:
        pipeline.modifiers.append(DeleteSelectedModifier())

        # Cluster analysis:
        pipeline.modifiers.append(ClusterAnalysisModifier(
            cutoff = 1.1, 
            sort_by_size = True, 
            unwrap_particles = True, 
            compute_com = True, 
            compute_gyration = True, 
            cluster_coloring = True))

        #Export Cluster analysis
        AvgClusterSize=[]
        StdClusterSize=[]
        for ThisFrame in range(1,50):
            data=pipeline.compute(frame=ThisFrame)
            ClustersListFile='{}/Analysis/ClustersFiles/ClustersList_{}_ts{:d}.dat'.format(folder,filepattern,data.attributes['Timestep'])
            #export_file(data.tables['clusters'], '{}/Analysis/ClustersFiles/ClustersList_{}_ts{:d}.dat'.format(folder,filepattern,data.attributes['Timestep']),'txt/table', precision=3)
            ClustersData=[ data.tables['clusters'].xy [:,[0,1,2,3,5,6,8]]
            ClustersDataDf=pd.DataFrame(ClustersData, columns=['ClusterID', 'Size', 'Xcm', 'Ycm', 'Ixx', 'Iyy', 'Ixy'])
            #data.attributes['AvgClusterSize']=data.tables['clusters'].xy()[:,1].mean()
            #data.attributes['StdClusterSize']=data.tables['clusters'].xy()[:,1].std()
            ClustersDataDf['Rg']=np.sqrt(ClustersDataDf['Ixx'].values+ClustersDataDf['Iyy'].values)
            ClustersDataDf.to_csv(sep=' ')
            
            AvgClusterSize.append(data.tables['clusters'].xy()[:,1].mean())
            StdClusterSize.append(data.tables['clusters'].xy()[:,1].std())
            
        export_file(pipeline, ClusterStatisticsFile,
                    'txt/attr', precision=3, columns=["Timestep", "ClusterAnalysis.largest_size"], multiple_frames=True)

        with open(ClusterStatisticsFile, 'r') as txt:
            lines = txt.readlines()
        lines[0]= lines[0][:-1]+' "AvgClusterSize" "StdClusterSize"\n'
        for i in range(1,len(lines)):
            lines[i]= lines[i][:-2]+' {:.3f} {:+.3f}\n'.format(AvgClusterSize[i-1],StdClusterSize[i-1])
        with open(ClusterStatisticsFile, 'w') as f:
            f.writelines(lines)
        print("Written "+ClusterStatisticsFile)


# In[5]:


# Erase files created by mistake for simulations that didn't complete time
# (it shouldn't happen anymore, because now there's a check on this in the previous cell)
for filestring in filestringList:
    folder='/Users/ivan/Dropbox/Mac/Documents/CurvatureLAT/Simulations/CylPla'
    laststring = subprocess.check_output("tail -n5 {}/Results/Log_{}.dat".format(folder,filestring), shell=True)
    ClusterStatisticsFile='{}/Analysis/ClustersFiles/ClustersStatistics_{}.dat'.format(folder,filestring)
    if ("Total wall time" not in laststring.decode("utf-8")) and (os.path.isfile(ClusterStatisticsFile)==True):
        print("Removed "+filestring)
        subprocess.check_output("rm {}/Analysis/ClustersFiles/Clusters*_{}*".format(folder,filestring), shell=True)    


# ## Plot cluster analysis and bond analysis

# ### Cython definitions

# In[6]:


get_ipython().run_line_magic('load_ext', 'Cython')
#Don't merge this cell with the next. The next contains cython code, as it starts by %%cython


# In[9]:


get_ipython().run_cell_magic('cython', '', '\n# -*- coding: utf-8 -*-\n# This program is free software: you can redistribute it and/or modify\n# it under the terms of the GNU General Public License as published by\n# the Free Software Foundation, either version 3 of the License, or\n# (at your option) any later version.\n#\n# This program is distributed in the hope that it will be useful,\n# but WITHOUT ANY WARRANTY; without even the implied warranty of\n# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n# GNU General Public License for more details.\n#\n# You should have received a copy of the GNU General Public License\n# along with this program.  If not, see <http://www.gnu.org/licenses/>.\n# Copyright Guglielmo Saggiorato 2018\n""" An "efficient" reader of LAMMPS dump files in cython.\nUsage:\n```import pyximport\npyximport.install()\nfrom load_lammps import read_lammps\nfor time, column_names,data in read_lammps(filepath):\n  pass\n```\n"""\n__author__ = "Guglielmo Saggiorato" #Modified by Ivan Palaia\n__copyright__ = "Copyright 2018, Guglielmo Saggiorato"\n__credits__ = ["Guglielmo Saggiorato",]\n__license__ = "GPLv3"\n__version__ = "1.0"\n__maintainer__ = "Guglielmo Saggiorato"\n__email__ = "astyonax@gmail.com"\n__status__ = "Production"\n\nimport numpy as np\ncimport numpy as np\nimport pandas as pd\ncimport cython\nfrom libc.stdio cimport FILE, fopen, fwrite, fscanf, fclose, fprintf, fseek, ftell, SEEK_SET, rewind, fread\n\n@cython.boundscheck(False)\n@cython.wraparound(False)\n@cython.nonecheck(False)\ndef read_traj(str fname):\n    """ Reads a LAMMPS dump file frame by frame yielding the data as a pandas DataFrame\nThis function automatically reads the fields/columns names and sets the pandas dataframe columns accordingly\nFor maximum performance, the bulk of the data is loaded without using any python object\nIn: fname [str]-- a file name as strings\nOut: frame time [int], columns name [list], data [pandas.DataFrame], boxbounds [pandas.DataFrame]\nUsage:\n```import pyximport\npyximport.install()\nfrom load_lammps import read_lammps\nfor time, column_names,data in read_lammps(filepath):\n  pass\n    """\n    cdef str line,values\n    cdef int i,N,T,j,cln,\n    cdef list toarr,columns,\n    cdef list columnsbox,\n    cdef int Nbox,clnbox\n    cdef double[:,:] data\n    cdef double[:,:] box\n    cdef double tmp\n    cdef FILE * ptr_r\n\n    fin = open(fname,\'r\')\n    ptr_r = fopen(bytes(fname.encode(\'utf-8\')), "r")\n    line = fin.readline()\n    while line:\n        if \'ITEM: TIMESTEP\' in line:\n            # begin new timestep\n            T = int(fin.readline())\n            N = 0\n            columns = []\n        if \'ITEM: NUMBER OF ATOMS\' in line:\n            N = int(fin.readline())\n        if \'ITEM: BOX BOUNDS\' in line:\n            Nbox=3\n            columnsbox = [\'lo\',\'hi\']\n            clnbox = len(columnsbox)\n            if not (Nbox and clnbox):\n                raise StopIteration\n            box = np.zeros((Nbox,clnbox),dtype=\'float64\')#,dtype=[(j,\'float64\') for j in columns])\n\n            fseek(ptr_r,int(fin.tell()),SEEK_SET)\n            # loop over x, y and z coordinate\n            for i in range(Nbox):\n                for j in range(clnbox):\n                    fscanf(ptr_r,"%le",&tmp)\n                    box[i,j] = tmp#toarr[j]\n\n            fin.seek(ftell(ptr_r))\n            qbox  = pd.DataFrame(np.asarray(box),columns=columnsbox)\n        if \'ITEM: ATOMS\' in line:\n            columns = line.split()[2:]\n            cln = len(columns)\n            if not (N and cln):\n                raise StopIteration\n            data = np.zeros((N,cln),dtype=\'float64\')#,dtype=[(j,\'float64\') for j in columns])\n\n            fseek(ptr_r,int(fin.tell()),SEEK_SET)\n            # loop over particles\n            for i in range(N):\n                for j in range(cln):\n                    fscanf(ptr_r,"%le",&tmp)\n                    data[i,j] = tmp#toarr[j]\n\n            fin.seek(ftell(ptr_r))\n            q  = pd.DataFrame(np.asarray(data),columns=columns)\n            yield T,columns,q,qbox\n\n        line = fin.readline()')


# In[10]:


get_ipython().run_cell_magic('cython', '', '\nimport numpy as np\ncimport numpy as np\nimport pandas as pd\ncimport cython\nfrom libc.stdio cimport FILE, fopen, fwrite, fscanf, fclose, fprintf, fseek, ftell, SEEK_SET, rewind, fread\n\n@cython.boundscheck(False)\n@cython.wraparound(False)\n@cython.nonecheck(False)\ndef read_bonds(str fname):\n    """ Reads a Bonds file as generated by my custom dump command for CuvatureLAT\n    """\n    cdef str line,values\n    cdef int i,N,T,j,cln,\n    cdef list toarr,columns,\n    cdef list columnsbox,\n    cdef int Nbox,clnbox\n    cdef double[:,:] data\n    cdef double[:,:] box\n    cdef double tmp\n    cdef FILE * ptr_r\n\n    fin = open(fname,\'r\')\n    ptr_r = fopen(bytes(fname.encode(\'utf-8\')), "r")\n    line = fin.readline()\n    while line:\n        if \'ITEM: TIMESTEP\' in line:\n            # begin new timestep\n            T = int(fin.readline())\n            N = 0\n            columns = []\n        if \'ITEM: NUMBER OF ENTRIES\' in line:\n            N = int(fin.readline())\n        if \'ITEM: BOX BOUNDS\' in line:\n            Nbox=3\n            columnsbox = [\'lo\',\'hi\',\'0\']\n            clnbox = len(columnsbox)\n            if not (Nbox and clnbox):\n                raise StopIteration\n            box = np.zeros((Nbox,clnbox),dtype=\'float64\')#,dtype=[(j,\'float64\') for j in columns])\n\n            fseek(ptr_r,int(fin.tell()),SEEK_SET)\n            # loop over x, y and z coordinate\n            for i in range(Nbox):\n                for j in range(clnbox):\n                    fscanf(ptr_r,"%le",&tmp)\n                    box[i,j] = tmp#toarr[j]\n\n            fin.seek(ftell(ptr_r))\n            qbox  = pd.DataFrame(np.asarray(box),columns=columnsbox)\n        if \'ITEM: ENTRIES\' in line:\n            #columns = line.split()[2:] \n            columns =[\'ind\',\'AtomID1\',\'AtomID2\',\'BondType\',\'BondLength\', \'BondEnergy\']\n            cln = len(columns)\n            if not (N and cln):\n                line = fin.readline()\n                q=[]\n                yield T,columns,q,qbox,N\n                continue\n                #raise StopIteration\n            data = np.zeros((N,cln),dtype=\'float64\')#,dtype=[(j,\'float64\') for j in columns])\n\n            fseek(ptr_r,int(fin.tell()),SEEK_SET)\n            # loop over bonds\n            for i in range(N):\n                for j in range(cln):\n                    fscanf(ptr_r,"%le",&tmp)\n                    data[i,j] = tmp#toarr[j]\n\n            fin.seek(ftell(ptr_r))\n            q  = pd.DataFrame(np.asarray(data),columns=columns)\n            yield T,columns,q,qbox,N\n\n        line = fin.readline()')


# ### Python functions definitions

# In[20]:


def get_densdatadic(densList, RcylList, realrange, Lc, nC, cbP):

    densdatadic={}

    for dens in densList:
        
        densStr = densitystring(dens)
        datadic={}

        for Rcyl in RcylList:

            times = []
            Nbonds = []
            NbondsVar = []
            ClustersTmpDfs = []
            realcounter=0

            for real in realrange:
                if Rcyl==0:
                    filestring='Lc{:d}_dens{:s}_p_nC{:d}_cbP{:.3f}_{:d}'.format(Lc, densStr, nC, cbP, real)
                else:
                    filestring='Lc{:d}_dens{:s}_c_nC{:d}_Rcyl{:.1f}_cbP{:.3f}_{:d}'.format(Lc, densStr, nC, Rcyl, cbP, real)            
                bonds_file='../Results/Bonds_{}.dat'.format(filestring)
                clustersstatistics_file='ClustersFiles/ClustersStatistics_{}.dat'.format(filestring)
                if (os.path.isfile(clustersstatistics_file)==False):
                    print("Inexistant {}".format(filestring))
                    continue
                realcounter+=1
                print(filestring)
                j=0
                for time, columns, data, box, N in read_bonds(bonds_file):
                    if realcounter==1:
                        times.append(time)
                        Nbonds.append(N)
                        NbondsVar.append(N**2)
                    else:
                        assert time==times[j], "ERROR: {:d}-th time ({:d}, {:d}) does not correspond for file {:s}".format(j,time,times[j],bonds_file)
                        Nbonds[j] += N
                        NbondsVar[j] += N**2
                    j+=1
                ClustersTmpDfs.append(pd.read_csv(clustersstatistics_file, names=['ts','LargestSize_{:d}'.format(real),'AvgSize_{:d}'.format(real)],skiprows=1, sep=' ' ))

            Nbonds=np.array(Nbonds)/realcounter
            NbondsStd=np.sqrt(np.array(NbondsVar)/realcounter - Nbonds**2)
            NbondsDf = pd.DataFrame(np.array([times,Nbonds,NbondsStd]).transpose(), columns=['ts','mean','std'])

            LargestSizeDf = ClustersTmpDfs[0].iloc[:,[0,1]]
            AvgSizeDf = ClustersTmpDfs[0].iloc[:,[0,2]]
            for df in ClustersTmpDfs[1:]:
                assert (LargestSizeDf['ts'].values==df['ts']).all(), print(df['ts'])
                LargestSizeDf = pd.concat((LargestSizeDf, df.iloc[:,[1]]) ,axis=1)
                AvgSizeDf = pd.concat((AvgSizeDf, df.iloc[:,[2]]) ,axis=1)
            LargestSizeDf.iloc[:,1:] /= Lc+3
            tmpmean = LargestSizeDf.iloc[:,1:].mean(axis=1)
            tmpstd = LargestSizeDf.iloc[:,1:].std(axis=1)
            LargestSizeDf['mean'] = tmpmean
            LargestSizeDf['std'] = tmpstd
            AvgSizeDf.iloc[:,1:] /= Lc+3
            tmpmean = AvgSizeDf.iloc[:,1:].mean(axis=1)
            tmpstd = AvgSizeDf.iloc[:,1:].std(axis=1)
            AvgSizeDf['mean'] = tmpmean
            AvgSizeDf['std'] = tmpstd

            for real in realrange:
                if 'LargestSize_{:d}'.format(real) in LargestSizeDf.columns:
                    del LargestSizeDf['LargestSize_{:d}'.format(real)]
                    del AvgSizeDf['AvgSize_{:d}'.format(real)]

            datadic['{:.1f}'.format(Rcyl)]={}
            datadic['{:.1f}'.format(Rcyl)]['Nbonds']=NbondsDf
            datadic['{:.1f}'.format(Rcyl)]['LargestSize']=LargestSizeDf
            datadic['{:.1f}'.format(Rcyl)]['AvgSize']=AvgSizeDf

        densdatadic[densStr] = datadic
    
    return densdatadic


#  ### Preliminary plots of number of bonds (Nbonds), largest cluster size (LargestSize), and average cluster size (AvgSize)

# In[21]:


## Import Bonds files into a DataFrame, for all configurations (planar and cylindrical, all Rcyl) 

densList=[0.01, 0.005, 0.001]
#densStringList = ["{:.3f}".format(dens) for dens in densList]

RcylList=[0, 12.5, 15, 20, 25, 50, 100] # 0 stands for plane
#RcylStringList = ["{:.1f}".format(Rcyl) for Rcyl in RcylList]

realrange=range(100,105)

Lc, nC, cbP = 24, 200, 0.001


# In[23]:


densdatadic=get_densdatadic(densList, RcylList, realrange, Lc, nC, cbP)


# In[51]:


## Plot Nbonds, LargestSize and AvgSize, at fixed density

dens = 0.001
densString='{:.3f}'.format(dens)
datadic=densdatadic[densString]

labeldic={}
colordic={}
VariableList=['Nbonds','LargestSize','AvgSize']
for i in range(0,len(RcylStringList)):
    colordic[RcylStringList[i]]=[
        'black','limegreen','mediumspringgreen','aquamarine','darkturquoise','cornflowerblue','mediumblue'][i]
    
# safety check on times
for key in VariableList:
    times=datadic['0.0'][key]['ts'].values
    for Rcyl in RcylStringList:
        assert (datadic[Rcyl][key]['ts'].values==times).all(), "Different times for {}, {}".format(Rcyl,key)
        if Rcyl=='0.0':
            labeldic[Rcyl]='Plane'
        else:
            labeldic[Rcyl]='Cylinder {}'.format(Rcyl)

# plot everything
for key in VariableList:
    plt.figure()
    for Rcyl in RcylStringList:
        df=datadic[Rcyl][key]
        plt.errorbar(df['ts'], df['mean'], yerr=df['std'], label=labeldic[Rcyl], color=colordic[Rcyl])
    plt.legend(frameon=False)
    plt.ylabel(key)
    plt.xlabel('time')
    plt.title('Lc{:d}, dens{:.3f}, nC{:d}, cbP{:.3f}'.format(Lc, dens, nC, cbP))
    plt.gcf().savefig('Graphs/{:s}CylPla_Lc{:d}_dens{:.3f}_nC{:d}_cbP{:.3f}.pdf'.format(key, Lc, dens, nC, cbP))
    plt.show()
    plt.close()


# In[52]:


## Plot Nbonds, LargestSize and AvgSize, at different densities

for dens in densList:
    densString='{:.3f}'.format(dens)
    datadic=densdatadic[densString]

    labeldic={}
    colordic={}
    VariableList=['Nbonds','LargestSize','AvgSize']
    for i in range(0,len(RcylStringList)):
        colordic[RcylStringList[i]]=[
            'black','limegreen','mediumspringgreen','aquamarine','darkturquoise','cornflowerblue','mediumblue'][i]

    # safety check on times
    for key in VariableList:
        times=datadic['0.0'][key]['ts'].values
        for Rcyl in RcylStringList:
            assert (datadic[Rcyl][key]['ts'].values==times).all(), "Different times for {}, {}".format(Rcyl,key)
            if Rcyl=='0.0':
                labeldic[Rcyl]='Plane'
            else:
                labeldic[Rcyl]='Cylinder {}'.format(Rcyl)

    # plot everything
    for i in range(0,len(VariableList)):
        print(i)
        plt.figure(i)
        key=VariableList[i]
        for Rcyl in RcylStringList:
            df=datadic[Rcyl][key]
            #plt.errorbar(df['ts'], df['mean'], yerr=df['std'], label=labeldic[Rcyl], color=colordic[Rcyl])
            if dens==densList[0]:
                plt.errorbar(df['ts'], df['mean'], label=labeldic[Rcyl], color=colordic[Rcyl])
            else:
                plt.errorbar(df['ts'], df['mean'], color=colordic[Rcyl])               
        plt.legend(frameon=False)
        plt.ylabel(key)
        plt.xlabel('time')
        plt.title('Lc{:d}, nC{:d}, cbP{:.3f}'.format(Lc, nC, cbP))
        plt.gcf().savefig('Graphs/{:s}CylPla_Lc{:d}_nC{:d}_cbP{:.3f}.pdf'.format(key,Lc, nC, cbP))
        #plt.show()


# ### Find density of linkers (and map to effective density)

# In[11]:


from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *


# Find average posizion of linkers

# 1) for planes

def get_linkers_z(filesstring):

    LenHist=14.0
    NBins=100

    # Data import:
    pipeline = import_file(filesstring, multiple_frames = True)

    # Manual modifications of the imported data objects:
    def modify_pipeline_input(frame: int, data: DataCollection):
        data.particles_.particle_types_.type_by_id_(1).color = (0.23056382085908295, 0.5206073090714886, 0.8764324406805524)
        data.particles_.particle_types_.type_by_id_(1).radius = 0.5
        data.particles_.particle_types_.type_by_id_(2).color = (0.0, 0.5898069733730068, 1.0)
        data.particles_.particle_types_.type_by_id_(2).radius = 0.5
        data.particles_.particle_types_.type_by_id_(3).color = (0.9969787136644541, 0.8339513237201496, 0.28714427405203324)
        data.particles_.particle_types_.type_by_id_(3).radius = 0.5
        data.particles_.particle_types_.type_by_id_(4).radius = 0.45
        #data.particles_.particle_types_.type_by_id_(7).radius = 0.55
    pipeline.modifiers.append(modify_pipeline_input)

    # Select type:
    pipeline.modifiers.append(SelectTypeModifier(types = {3}))

    # Compute property:
    pipeline.modifiers.append(ComputePropertyModifier(
        expressions = ('Position.Z',), 
        output_property = 'MyHeight', 
        only_selected = True))

    # Histogram:
    pipeline.modifiers.append(HistogramModifier(
        property = 'MyHeight', 
        bin_count = NBins, 
        fix_xrange = True, 
        xrange_end = LenHist, 
        only_selected = True))

    # Time averaging:
    pipeline.modifiers.append(TimeAveragingModifier(operate_on = 'table:histogram[MyHeight]'))

    # Get value
    data = pipeline.compute()
    histogram = data.tables['histogram[MyHeight][average]'].xy()
    bins = histogram[:,0]
    histovalue = histogram[:,1]
    binLen = LenHist/NBins
    norm = sum(binLen*histovalue)
    avg = sum(binLen*histovalue*bins) / norm
    std = np.sqrt( sum(binLen*histovalue*bins**2) / norm - avg**2 )
    #print([bins, histovalue, avg, std])
    #plt.plot(bins,histovalue)
    return [avg, std]




# 2) for cylinders

def get_linkers_deltaR(filesstring, Rcyl):

    LenHist=14.0
    NBins=100

    # Data import:
    pipeline = import_file(filesstring, multiple_frames = True)

    # Manual modifications of the imported data objects:
    def modify_pipeline_input(frame: int, data: DataCollection):
        data.particles_.particle_types_.type_by_id_(1).color = (0.23056382085908295, 0.5206073090714886, 0.8764324406805524)
        data.particles_.particle_types_.type_by_id_(1).radius = 0.5
        data.particles_.particle_types_.type_by_id_(2).color = (0.0, 0.5898069733730068, 1.0)
        data.particles_.particle_types_.type_by_id_(2).radius = 0.5
        data.particles_.particle_types_.type_by_id_(3).color = (0.9969787136644541, 0.8339513237201496, 0.28714427405203324)
        data.particles_.particle_types_.type_by_id_(3).radius = 0.5
        data.particles_.particle_types_.type_by_id_(4).radius = 0.45
        #data.particles_.particle_types_.type_by_id_(7).radius = 0.55
    pipeline.modifiers.append(modify_pipeline_input)

    # Select type:
    pipeline.modifiers.append(SelectTypeModifier(types = {3}))

    # Compute property:
    pipeline.modifiers.append(ComputePropertyModifier(
        expressions = ('{:.3f}-sqrt(Position.X^2 + Position.Y^2)'.format(Rcyl),), 
        output_property = 'MyDistance',  
        only_selected = True))

    # Histogram:
    pipeline.modifiers.append(HistogramModifier(
        property = 'MyDistance', 
        bin_count = NBins, 
        fix_xrange = True, 
        xrange_end = LenHist, 
        only_selected = True))

    # Time averaging:
    pipeline.modifiers.append(TimeAveragingModifier(operate_on = 'table:histogram[MyDistance]'))

    # Get value
    data = pipeline.compute()
    histogram = data.tables['histogram[MyDistance][average]'].xy()
    bins = histogram[:,0]
    histovalue = histogram[:,1]
    binLen = LenHist/NBins
    norm = sum(binLen*histovalue)
    avg = sum(binLen*histovalue*bins) / norm
    std = np.sqrt( sum(binLen*histovalue*bins**2) / norm - avg**2 )
    #print([bins, histovalue, avg, std])
    #plt.plot(bins,histovalue)
    return [avg, std]


# In[16]:


# Una tantum

Table=[]
cbP=0.010
for dens in densList:
    for Rcyl in RcylList:
        if Rcyl==0:
            filesstring = '../Results/Traj_Lc24_dens{:.3f}_p_nC800_cbP{:.3f}_10*.xyz'.format(dens,cbP)
            avg, std = get_linkers_z(filesstring)
            EquivalentPlaneDensity = dens 
        else:
            filesstring = '../Results/Traj_Lc24_dens{:.3f}_c_nC800_Rcyl{:.1f}_cbP{:.3f}_10*.xyz'.format(dens, Rcyl,cbP)
            avg, std = get_linkers_deltaR(filesstring, Rcyl)
            EquivalentPlaneDensity = dens*Rcyl/(Rcyl-avg)    # density correction
        print([dens, Rcyl, avg, std, EquivalentPlaneDensity])
        Table.append([dens, Rcyl, avg, std, EquivalentPlaneDensity])
DistanceLinkersWallDf = pd.DataFrame(Table, columns=['dens','Rcyl','avg','std','EquivalentPlaneDensity'])
DistanceLinkersWallDf.to_csv('DistanceLinkersWallDf_Lc{:d}_cbP{:.3f}.dat'.format(Lc,cbP), index = False, header=True, sep=' ')   


# In[29]:


cbP=0.010
DistanceLinkersWallDf = pd.read_csv('DistanceLinkersWallDf_Lc{:d}_cbP{:.3f}.dat'.format(Lc,cbP), sep=' ')
equivdens = DistanceLinkersWallDf['EquivalentPlaneDensity'].values.round(6)
DistanceLinkersWallDf


# ### Compare cylinders with planes at equal effective density

# In[12]:





# In[128]:


densList=[0.01    , 0.023105, 0.018282, 0.014678, 0.013384, 0.011418,
       0.010658, 0.005   , 0.01075 , 0.008676, 0.007199, 0.006594,
       0.005683, 0.00532 , 0.001   , 0.001878, 0.001637, 0.00141 ,
       0.001304, 0.001134, 0.001063]
RcylList = [0]
nC=200
cbP=0.001
datadensdicEqPla = get_densdatadic(densList, RcylList, range(100,105), Lc, nC, cbP)

densList=[0.01, 0.005 , 0.001]
RcylList = [0.0, 12.5, 15, 20, 25, 50, 100]
nC=200
cbP=0.001
datadensdicCyl = get_densdatadic(densList, RcylList, range(100,105), Lc, nC, cbP)


# In[125]:


colors=['k','b','r','m']
i=0
RcylList = [12.5]
VariableList=['Nbonds','LargestSize','AvgSize']
for VarStr in VariableList:
    plt.figure()
    i=0
    for dens in densList:
        densStr=densitystring(dens)
        for Rcyl in RcylList:
            RcylStr='{:.1f}'.format(Rcyl)
            df= datadensdicCyl[densStr][RcylStr][VarStr]
            plt.errorbar(df['ts'].values, df['mean'].values, yerr=df['std'].values, c=colors[i], lw=1, fmt='-',label='dens = {:s}, Rcyl = {:s}'.format(densStr,RcylStr))
            df= datadensdicCyl[densStr]['0.0'][VarStr]
            plt.errorbar(df['ts'].values, df['mean'].values, c=colors[i], lw=1, fmt=':')
            densEqPla = DistanceLinkersWallDf[(DistanceLinkersWallDf['dens']==dens) & (DistanceLinkersWallDf['Rcyl']==Rcyl)]['EquivalentPlaneDensity'].iloc[0]
            df= datadensdicEqPla[densitystring(densEqPla)]['0.0'][VarStr]
            plt.errorbar(df['ts'].values, df['mean'].values, fmt=':', lw=3, c=colors[i])
            i+=1
    plt.legend(frameon=False)
    plt.xlabel('time')
    plt.ylabel(VarStr)
    plt.savefig('Graphs/{:s}CylPla_EqDensPla_nC{:d}_cbP{:.3f}.pdf'.format(VarStr,nC,cbP))
    plt.show()


# In[66]:


df['std'].values


# In[ ]:




