import numpy as np
import subprocess
import os
import re
import pandas as pd
import argparse



if __name__ == "__main__":

    TrajLogFilePath = ''
    timestep = -1

    # Get input
    parser = argparse.ArgumentParser(description="Script to extract single clusters from a pre-run simulations, as config files for a new simulation.")
    parser.add_argument('--TrajLogFilePath', '--trajfile', dest='TrajLogFilePath', action='store', type=str, default=TrajLogFilePath, help='Trajectory file to extract clusters from.')
    parser.add_argument('--timestep', '-ts', dest='timestep', action='store', type=int, default=timestep, help='Timestep at which clusters are to be extracted.')

    args = parser.parse_args()
    TrajLogFilePath=args.TrajLogFilePath
    timestep=args.timestep

    # Define string managers
    r1 = re.compile('[_]')
    r2 = re.compile(r'(([-+]?\d+\.\d+)|([-+]?\d+))')


    filepattern = subprocess.check_output(
            " echo {:s} | sed 's/.*\/TrajLog_//' | sed 's/.xyz//'".format(TrajLogFilePath), shell=True)
    filepattern = filepattern.decode("utf-8")[:-1]
    Test5Path = subprocess.check_output(
        " echo {:s} | sed 's/Results_.*//' | sed 's/.xyz//'".format(TrajLogFilePath), shell=True)
    Test5Path = Test5Path.decode("utf-8")[:-2]
    ClustersListFilePath = '{:s}/Analysis/ClustersFiles/ClustersList_{:s}_ts{:d}.dat'.format(Test5Path, filepattern, timestep)
    ClustersListDf = pd.read_csv(ClustersListFilePath, sep=' ')
    AvgRg = ClustersListDf['Rg'].mean()

    for s in r1.split(filepattern):
        if s.startswith('nA'):
            nA = int(r2.split(s)[1])
        if s.startswith('dens'):
            dens = float(r2.split(s)[1])
        if s[0].isdigit():
            real = int(s)
    TwoL=np.sqrt(nA/dens)



        self.num_A = num_A

        self.q_A = q_A

        self.sigma = sigma
        self.PatchRadialDistance = PatchRadialDistance

        self.lattice_factor = lattice_factor
        self.side = side
     
        self.lattice_constant = lattice_factor*self.sigma
        self.Lx = self.side*self.lattice_constant
        self.Ly = self.side*self.lattice_constant
        self.Lz = 0.1

        self.coords = ["\nAtoms \n \n"]
        self.bonds = ["\nBonds \n \n"]
        self.angles= ["\nAngles \n \n"]

        # m is number of atoms and k is number of 3 atom-molecules
        self.numAll = 0
        self.numMol = 0

        self.numBonds = 0

        self.numAngles = 0

        self.numTypes = 0
        self.type_mass_list =[]


        def init_lattice_2d(lat_con):
            num_x = int(self.side)
            num_y = int(self.side)
            lattice_out = np.zeros(shape=(num_x*num_y, 3))
            counter = 0
            for i in range(num_x):
                for j in range(num_y):
                    lattice_out[counter, 0] = lat_con * (i+0.5) - lat_con*num_x * 0.5
                    lattice_out[counter, 1] = lat_con * (j+0.5) - lat_con*num_y * 0.5
                    lattice_out[counter, 2] = 0.0
                    counter += 1
            return lattice_out


        self.lattice_sites = init_lattice_2d(self.lattice_constant)
        assert self.lattice_factor>1.1, "ERROR: lattice_factor<1.1, too many particles."
        assert self.numAll <= self.lattice_sites.shape[0], "ERROR: not enough lattice sites, increase box size."

        lattice_inds = np.arange(self.lattice_sites.shape[0])

        np.random.seed(100+real)
        self.start_inds_A = np.random.choice(lattice_inds, size = self.num_A, replace = False)


    def make_A(self):

        particles_list = [1,2]+[3 for i in range(3,self.q_A+2)] # range(1,self.q_A+2)
        particle_types = len(particles_list)
        self.numTypes += len(np.unique(particles_list))

        MassAlpha = 0.1 * (self.sigma/self.PatchRadialDistance)**2  # fraction of mass due to patches such that moment of inertia equals that of a sphere
        MassCentralAtom = 1-MassAlpha
        MassPatch = MassAlpha/self.q_A
        for i in range(particle_types):
            if i+1==1:
                to_add = [particles_list[i], MassCentralAtom]
            else:
                to_add = [particles_list[i], MassPatch]
            if to_add not in self.type_mass_list:
                self.type_mass_list.append(to_add)

        for i in range(len(self.start_inds_A)):
            indBuf = self.start_inds_A[i]

            self.numMol += 1
            x = self.lattice_sites[indBuf, 0]
            y = self.lattice_sites[indBuf, 1]
            z = self.lattice_sites[indBuf, 2]

            RandomPhase = np.random.rand() * 2 * np.pi
            for n in range(self.q_A+1):
                self.numAll += 1
                if (n == 0):
                    self.coords.append(
                        "\t " + str(self.numAll) + " " + str(self.numMol) + " "+str(particles_list[n])+" 0 " + str(x) + " " + str(y) + " " + str(
                            z) + " 0 0 0 \n")
                position_on_circle = self.PatchRadialDistance*np.array([[ np.cos(RandomPhase + i*2*np.pi/self.q_A), np.sin(RandomPhase + i*2*np.pi/self.q_A) ] for i in range(0,self.q_A)])
                if n > 0:
                    self.coords.append(
                        "\t " + str(self.numAll) + " " + str(self.numMol) + " " +  str(particles_list[n]) + " 0 " + str(x - position_on_circle[n-1,0]) + " " + str(y - position_on_circle[n-1,1]) + " " + str(
                            z) + " 0 0 0 \n")
                    # if n==1:
                    #     self.numBonds += 1
                    #     self.bonds.append(
                    #         "\t {:d} {:d} {:d} {:d} \n".format(self.numBonds, 1, self.numAll-1, self.numAll)
                    #     )
