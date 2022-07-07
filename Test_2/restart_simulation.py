# From code for fluid clusters as described in:
# I. Palaia and A. Saric, Controlling cluster size in 2D phase-separating binary mixtures with specific interactions, J. Chem. Phys. (2022)
#
import argparse
from make_initconfig import make_particles
import os
import glob
import numpy as np



def write_in_script(sigma, numParticleTypes, PatchRange, PatchStrength, IsotropicAttrRange, IsotropicAttrStrength, real,
                    RunSteps, dumpevery, filePattern, RestartFilePath, ResultsFolder, extForce):
    seed = 100 + real
    filename = "Input/Scripts/Input_{:s}.in".format(filePattern)
    Temperature = 1.0
    LangevinDamping = 0.08
    TimestepToTau0 = 0.008

    f = open(filename, "w")

    #################################
    # INITIALISATION AND INTERACTIONS
    #################################

    f.write("log                    {:s}/Log_{:s}.dat \n\n".format(ResultsFolder, filePattern))
    f.write("read_restart           {:s}\n".format(RestartFilePath))
    #f.write("units                  lj \n")
    #f.write("dimension              2 \n")
    #f.write("atom_style             full \n")
    #f.write("boundary               p p p \n")
    #f.write("read_data              {:s} \n\n".format(configName))

    f.write("group                  Central type 1\n")
    f.write("group                  Patch type 2\n\n")
    f.write("group                  CentralAndPatch type 1 2\n\n")

    f.write("neighbor               0.5 bin\n")
    f.write("neigh_modify           every 1 delay 1\n")
    f.write("neigh_modify           exclude molecule/intra all \n\n")

    GlobalCutoff = sigma + IsotropicAttrRange
    f.write("pair_style             cosine/squared {:.4f} \n".format(GlobalCutoff))
    f.write("pair_coeff             * * 0.0 1.0 1.0 wca \n")

    f.write("pair_coeff             1 1  {:.1f} {:.3f} {:.3f} wca \n".format(IsotropicAttrStrength, sigma,
                                                                             sigma + IsotropicAttrRange))

    for id1 in range(2, numParticleTypes + 1):
        for id2 in range(id1, numParticleTypes + 1):
            f.write("pair_coeff             {:d} {:d}  {:.1f} {:.3f} {:.3f} \n".format(id1, id2, PatchStrength, 0.0,
                                                                                       PatchRange))  # mirrored cosine squared if PatchStrength<0 (purely repulsive), or additional attraction if PatchStrength>0


    ######################
    # COMPUTE FORCE TO ADD
    ######################

    # compute per-atom positions, if necessary shifted to account for pbc
    f.write("compute                cMol all chunk/atom molecule \n")
    # extract positions of central atom and 1sdt patch for each molecule
    f.write("compute                cCMCentral Central com/chunk cMol \n")
    f.write("compute                cCMPatch Patch com/chunk cMol \n")
    # create two atom variables that return the position of the central atom of the molecule and the position of the 1st patch of the molecule
    f.write("compute                cCentral all chunk/spread/atom cMol c_cCMCentral[*] \n")
    f.write("compute                cPatch all chunk/spread/atom cMol c_cCMPatch[*] \n")
    # compute unwrapped coordinates (because centre of mass is unwrapped in LAMMPS)
    f.write("compute                cxu all property/atom xu \n")
    f.write("compute                cyu all property/atom yu \n")

    if np.isclose(extForce, 0, atol=1e-6, rtol=0)==False:
        # compute forces
        f.write("# extForce = {:f} \n".format(extForce))
        # here x and y are positions of the patch
        f.write("variable               fxPatch   atom -{:f}*(c_cyu-c_cCentral[2]) \n".format(extForce / PatchRadialDistance))
        f.write("variable               fyPatch   atom  {:f}*(c_cxu-c_cCentral[1]) \n".format(extForce / PatchRadialDistance))
        # now x and y are positions of the central atom
        f.write("variable               fxCentral atom  {:f}*(c_cPatch[2]-c_cyu) \n".format(extForce / PatchRadialDistance))
        f.write("variable               fyCentral atom -{:f}*(c_cPatch[1]-c_cxu) \n\n".format(extForce / PatchRadialDistance))

    # # Move Central and Patch stored positions such that if patch and central are across boundaries, computed force is still correct
    # # Remember that vCentralx and vCentraly are only used for patch atoms and viceversa. Note that ((x>0)-(x<0)) = sign(x)
    # f.write("variable               vCentralx atom \"(abs(c_cCentral[1]-x)<0.5*lx) * c_cCentral[1] + (abs(c_cCentral[1]-x)>=0.5*lx) * (c_cCentral[1]+((x>0)-(x<0))*lx)\"  \n")  #.format(0.5*Lx, 0.5*Lx, Lx))
    # f.write("variable               vCentraly atom \"(abs(c_cCentral[2]-y)<0.5*ly) * c_cCentral[2] + (abs(c_cCentral[2]-y)>=0.5*ly) * (c_cCentral[2]+((y>0)-(y<0))*ly)\"  \n")
    # f.write("variable               vPatchx   atom \"(abs(c_cPatch[1]  -x)<0.5*lx) * c_cPatch[1] +   (abs(c_cPatch[1]  -x)>=0.5*lx) * (c_cPatch[1]  +((x>0)-(x<0))*lx)\"  \n")
    # f.write("variable               vPatchy   atom \"(abs(c_cPatch[2]  -y)<0.5*ly) * c_cPatch[2] +   (abs(c_cPatch[2]  -y)>=0.5*ly) * (c_cPatch[2]  +((y>0)-(y<0))*ly)\"  \n")
    # # compute forces
    # f.write("# extForce = {:f} \n".format(extForce))
    # # here x and y are positions of the patch
    # f.write("variable               fxPatch   atom  {:f}*(y-v_vCentraly) \n".format(extForce / PatchRadialDistance))
    # f.write("variable               fyPatch   atom -{:f}*(x-v_vCentralx) \n".format(extForce / PatchRadialDistance))
    # # now x and y are positions of the central atom
    # f.write("variable               fxCentral atom -{:f}*(v_vPatchy-y) \n".format(extForce/PatchRadialDistance))
    # f.write("variable               fyCentral atom  {:f}*(v_vPatchx-x) \n\n".format(extForce/PatchRadialDistance))

    ##################
    # TIME INTEGRATION
    ##################

    f.write("fix                    fLANG all langevin {:.3f} {:.3f} {:.3f} {:d} \n".format(Temperature, Temperature,
                                                                                            LangevinDamping, seed))
    f.write("fix                    fRigidNVE all rigid/nve molecule\n")
    f.write("fix                    fEnforce2d all enforce2d\n")

    if np.isclose(extForce, 0, atol=1e-6, rtol=0)==False:
        f.write("fix                    fTorqueCentral Central addforce v_fxCentral v_fyCentral 0 \n")
        f.write("fix                    fTorquePatch Patch addforce v_fxPatch v_fyPatch 0 \n\n")

    ########################
    # VARIABLES FOR ANALYSIS
    ########################

    if np.isclose(extForce, 0, atol=1e-6, rtol=0)==False:
        # Check modulus of forces (debug)
        f.write("variable               fCentralModulus atom sqrt(v_fxCentral^2+v_fyCentral^2)\n")
        f.write("variable               fPatchModulus atom sqrt(v_fxPatch^2+v_fyPatch^2)\n\n")

    # compute per-molecule angular momentum, torque, velocity and angular velocity
    f.write("compute                cAngMomMol all angmom/chunk cMol \n")
    f.write("compute                cTorqueMol all torque/chunk cMol \n")
    f.write("compute                cOmegaMol all omega/chunk cMol \n")
    f.write("compute                cVelCM Central vcm/chunk cMol \n\n")

    # compute total angular momentum and torque, wrt centre of box
    f.write("variable               vAngMomTot equal angmom(all,z) \n")
    f.write("fix                    fAngMomTotAvgt all ave/time 10 100 {:d} v_vAngMomTot \n".format(dumpevery))
    f.write("variable               vTorqueTot equal torque(all,z) \n")
    f.write("fix                    fTorqueTotAvgt all ave/time 10 100 {:d} v_vTorqueTot \n".format(dumpevery))
    # compute average angular momentum, torque and angular velocity, about centre of mass of each molecule
    f.write("compute                cAngMomMolSpread all chunk/spread/atom cMol c_cAngMomMol[3] \n")
    f.write("compute                cAngMomAvgm Central reduce ave c_cAngMomMolSpread \n")
    f.write("fix                    fAngMomAvgmAvgt all ave/time 10 100 {:d} c_cAngMomAvgm \n".format(dumpevery))
    f.write("compute                cTorqueMolSpread all chunk/spread/atom cMol c_cTorqueMol[3] \n")
    f.write("compute                cTorqueAvgm Central reduce ave c_cTorqueMolSpread \n")
    f.write("fix                    fTorqueAvgmAvgt all ave/time 10 100 {:d} c_cTorqueAvgm \n".format(dumpevery))
    f.write("compute                cOmegaMolSpread all chunk/spread/atom cMol c_cOmegaMol[3] \n")
    f.write("compute                cOmegaAvgm Central reduce ave c_cOmegaMolSpread \n")
    f.write("fix                    fOmegaAvgmAvgt all ave/time 10 100 {:d} c_cOmegaAvgm \n".format(dumpevery))
    f.write("variable               vTimestepsPerTurn equal 2*PI/{:f}/(f_fOmegaAvgmAvgt+1e-99) \n\n".format(TimestepToTau0))

    if np.isclose(extForce, 0, atol=1e-6, rtol=0)==False:
        # compute quantity that must be 0 if addforce is correct
        f.write("variable               vAddforceCheck0 equal f_fTorqueCentral+f_fTorquePatch \n")
        f.write("fix                    fAddforceCheck0 all ave/time 10 100 {:d} v_vAddforceCheck0 \n\n".format(dumpevery))

    # compute clusters and angular momentum per cluster
    # f.write("compute                cClusters Central cluster/atom {:f} \n".format(sigma+IsotropicAttrRange))
    # f.write("compute                     Central chunk/atom cMol \n")
    # f.write("compute                cClustChunk Central chunk/atom c_cClusters nchunk every compress yes")
    # f.write("compute                cAngMomClust all angmom/chunk cClusterMol  ")

    ########
    # THERMO
    ########

    f.write("thermo                 {:d} \n".format(int(dumpevery)))
    if np.isclose(extForce, 0, atol=1e-6, rtol=0)==False:
        f.write("thermo_style           custom step temp press etotal epair f_fAddforceCheck0 f_fTorqueTotAvgt f_fAngMomTotAvgt f_fTorqueAvgmAvgt f_fAngMomAvgmAvgt f_fOmegaAvgmAvgt v_vTimestepsPerTurn \n")
    else:
        f.write("thermo_style           custom step temp press etotal epair f_fTorqueTotAvgt f_fAngMomTotAvgt f_fTorqueAvgmAvgt f_fAngMomAvgmAvgt f_fOmegaAvgmAvgt v_vTimestepsPerTurn \n")
    f.write("thermo_modify          flush yes \n")
    f.write("timestep               {:f} \n\n".format(TimestepToTau0))

    ########
    # DUMP
    ########

    # f.write("fix                    prova all momentum 1 linear 1 1 0\n")
    f.write("variable        	dumpts equal logfreq(1000,9,10)\n")
    f.write("dump                   1 CentralAndPatch custom {:d} {:s}/Traj_{:s}.xyz id type mol x y z vx vy vz \n".format(dumpevery, ResultsFolder, filePattern))
    f.write("dump_modify            1  sort id  flush yes  first yes\n")
    f.write("dump                   1Log CentralAndPatch custom {:d} {:s}/TrajLog_{:s}.xyz id type mol x y z vx vy vz \n".format(dumpevery, ResultsFolder, filePattern))
    f.write("dump_modify            1Log  sort id  every v_dumpts  first yes\n")

    f.write("fix                    2 all ave/time 1 1 {:d} c_cCMCentral[1] c_cCMCentral[2] c_cVelCM[1] c_cVelCM[2] c_cTorqueMol[3] c_cAngMomMol[3] c_cOmegaMol[3] mode vector file {:s}/RotationStatsMolecule_{:s}.dat \n".format(
            dumpevery, ResultsFolder, filePattern))

    if np.isclose(extForce, 0, atol=1e-6, rtol=0)==False:
        f.write("dump                   3 CentralAndPatch custom {:d} {:s}/CheckTorque_{:s}.dat id type mol x y c_cCentral[1] c_cCentral[2] c_cPatch[1] c_cPatch[2] c_cxu c_cyu v_fxCentral v_fyCentral v_fxPatch v_fyPatch v_fCentralModulus v_fPatchModulus \n".format(
                dumpevery, ResultsFolder, filePattern))
        f.write("dump_modify            3 sort id \n\n")


    ########
    # RUN
    ########

    f.write("restart                1000000 {:s}/RestartA_{:s} {:s}/RestartB_{:s}\n".format(ResultsFolder,filePattern,ResultsFolder,filePattern))
    f.write("run                    {:d} upto\n".format(RunSteps))

    f.close()

    return filename




def writeRunScript(filePattern, Resultsfolder, inputfilename, MaxRunTime, MPInum):
    runfilename = "runscript_{}.sh".format(filePattern)
    if MPInum == 1:
        JobName = "SerJob_{}".format(filePattern)
    else:
        assert MPInum <= 32, "ERROR: Too many cores requested."
        JobName = "ParJob_{}".format(filePattern)
    OutFileName = "{:s}/Out_{}.dat".format(Resultsfolder, filePattern)

    f = open(runfilename, "w")
    f.write("""#!/bin/bash
#
#SBATCH --job-name={}
#SBATCH --output={}
""".format(JobName, OutFileName))
    if MPInum > 1:
        # f.write("#SBATCH -c {:d}\n".format(MPInum))
        f.write("#SBATCH --ntasks-per-node={:d} \n".format(MPInum))
        f.write("#SBATCH --nodes=1 \n")
        f.write("#SBATCH --ntasks={:d} \n".format(MPInum))
        f.write("#SBATCH --mem-per-cpu=800M \n")
    if MPInum == 1:
        f.write("#SBATCH --mem=800M \n")
    #f.write("#SBATCH --exclude beta233,leonid63 \n")
    f.write("#SBATCH --time={:d}:00:00 \n".format(MaxRunTime))
    f.write("#SBATCH --mail-user=ivan.palaia@ist.ac.at \n")
    f.write("#SBATCH --mail-type=NONE \n")
    # f.write("#SBATCH --exclude leonid63,leonid64 \n")
    f.write("#SBATCH --no-requeue \n")
    # if MPInum>1:
    #	if.write("#SBATCH --partition=beta \n")
    f.write("#SBATCH --export=NONE \n")
    f.write("unset SLURM_EXPORT_ENV \n")
    if MPInum>1:
        f.write("module load openmpi/3.1.3\n")
    f.write("export OMP_NUM_THREADS=1 \n")
    if MPInum == 1:
        f.write(
            "srun --cpu_bind=verbose  ~/my_src/lammps/lammps-29Sep2021/src/lmp_serial -in {} ".format(inputfilename))
    if MPInum > 1:
        f.write("srun --cpu_bind=verbose  ~/my_src/lammps/lammps-29Sep2021/src/lmp_mpi -in {} ".format(inputfilename))

    f.close()

    return runfilename




if __name__ == "__main__":

    dumpevery = 100000
    RunSteps = 2000000
    MaxRunTime = 240
    MPInum = 16

    # Simulation parameters
    num_A = 20
    density_A = 0.2  # number of colloids per sigma^-2. A density of 0.25 corresponds to a packing fraction of 0.2
    q_A = 5
    sigma = 1.0

    IsotropicAttrRange = 0.15
    IsotropicAttrStrength = 20.0

    PatchRange = 0.15
    PatchRadialDistance = 0.5
    PatchStrength = -IsotropicAttrStrength

    Temperature = 1

    extTorque = 0.1

    real = 0

    parser = argparse.ArgumentParser(description="Script for spinning, aggregating colloids.")
    # parser.add_argument('--eps_pp', '-eps_pp', dest='eps_pp', action='store', type=float, default=eps_pp, help='Energy depth of repulsion between proteins')
    # parser.add_argument('--eps_ll', '-eps_ll', dest='eps_ll', action='store', type=float, default=eps_ll, help='Energy depth of attraction/repulsion between cross patchs')
    parser.add_argument('--numA', '-nA', dest='num_A', action='store', type=int, default=num_A,
                        help='Number of aggregating colloids.')
    parser.add_argument('--density_A', '-density_A', '-dens', dest='density_A', action='store', type=float,
                        default=density_A,
                        help='Surface density of colloids in units of sigma^-2 (multiply by 0.25*pi*sigma^2 to obtain packing fraction).')
    parser.add_argument('--q_A', '-qA', dest='q_A', action='store', type=int, default=q_A,
                        help='Number of patches on 1 colloid.')
    # parser.add_argument('--patch_sigma', '-ps', dest='patch_sigma', action='store', type=float, default=patch_sigma, help='Diameter of patch')
    parser.add_argument('--PatchRange', '-pr', '-rp', dest='PatchRange', action='store', type=float, default=PatchRange,
                        help='Radius of patches.')
    parser.add_argument('--PatchRadialDistance', '-pd', '-dp', dest='PatchRadialDistance', action='store', type=float,
                        default=PatchRadialDistance, help='Radial distance of patches from centre of colloid.')
    parser.add_argument('--PatchStrength', '-ps', '-ep', dest='PatchStrength', action='store', type=float,
                        default=PatchStrength, help='Radial distance of patches from centre of colloid.')
    parser.add_argument('--IsotropicAttrRange', '-ra', dest='IsotropicAttrRange', action='store', type=float,
                        default=IsotropicAttrRange,
                        help='Attraction range of the cosine-squared isotropic attraction between molecules.')
    parser.add_argument('--IsotropicAttrStrength', '-ea', dest='IsotropicAttrStrength', action='store', type=float,
                        default=IsotropicAttrStrength,
                        help='Strength of the cosine-squared isotropic attraction between molecules.')
    parser.add_argument('--real', '-real', dest='real', action='store', type=int, default=real,
                        help='Number of realisation (statistics), sets the seed.')
    parser.add_argument('--Temperature', '-T', dest='Temperature', action='store', type=float, default=Temperature,
                        help='Temperature.')
    parser.add_argument('-runsteps', '--RunSteps', dest='RunSteps', action='store', type=int, default=RunSteps,
                        help='Number of time steps.')
    parser.add_argument('-dumpevery', '--dumpevery', dest='dumpevery', action='store', type=int, default=dumpevery,
                        help='Number of skipped time steps in dump file.')
    parser.add_argument('-extTorque', '--extTorque', '-eT', dest='extTorque', action='store', type=float,
                        default=extTorque, help='Number of skipped time steps in dump file.')
    parser.add_argument('-MRT', '--MaxRunTime', dest='MaxRunTime', action='store', type=int, default=MaxRunTime,
                        help="Maximum run time on cluster (integer, hours). Default {}.".format(MaxRunTime))
    parser.add_argument('-MPI', '--MPInum', dest='MPInum', action='store', type=int, default=MPInum,
                        help="Number of cores used for simulation. Default {}.".format(MPInum))
    args = parser.parse_args()

    num_A = args.num_A
    density_A = args.density_A
    q_A = args.q_A
    PatchRange = args.PatchRange
    PatchRadialDistance = args.PatchRadialDistance
    PatchStrength = args.PatchStrength
    IsotropicAttrRange = args.IsotropicAttrRange
    IsotropicAttrStrength = args.IsotropicAttrStrength
    real = args.real
    RunSteps = args.RunSteps
    dumpevery = args.dumpevery
    extTorque = args.extTorque
    MaxRunTime = args.MaxRunTime
    MPInum = args.MPInum

    extForce = extTorque / PatchRadialDistance

    # Initial configuration

    side = np.ceil(np.sqrt(num_A))
    lattice_factor = np.sqrt(
        1.0 * num_A / (side * side) / density_A)  # Ensures that A surface density is exactly 0.0075 (protein_sigma)^-2
    densitycheck = num_A / (side * lattice_factor) ** 2  # In units of (protein_sigma)^-2
    assert np.isclose(density_A, densitycheck), "ERROR: Densities do not correspond."

    system = make_particles(sigma, PatchRadialDistance, num_A, q_A, side, lattice_factor, real)
    system.make_A()

    OriginalFilePattern = "qA{:d}_dp{:.2f}_dens{:.2f}_eT{:.2f}_nA{:d}_rp{:.2f}_ra{:.2f}_ep{:.1f}_ea{:.1f}_T{:.1f}_{:d}".format(
        q_A, PatchRadialDistance, density_A, extTorque, num_A, PatchRange, IsotropicAttrRange, PatchStrength, IsotropicAttrStrength, Temperature, real)
    ResultsFolder = "Results_qA{:d}_dp{:.2f}_dens{:.3f}_eT{:.2f}".format(q_A, PatchRadialDistance, density_A, extTorque)

    if not os.path.exists(ResultsFolder):
        os.makedirs(ResultsFolder)

    # Check that restart file exists and determine new filePattern
    filelist = glob.glob("{:s}/Restart*_{:s}".format(ResultsFolder,OriginalFilePattern))
    assert len(filelist)!=0, "ERROR: No restart file for {:s}".format(OriginalFilePattern)
    nexti=2
    for i in np.arange(20,1,-1):
        listi = [x for x in filelist if any(s in x for s in ['RestartA_{:d}_'.format(i), 'RestartB_{:d}_'.format(i)]) ]
        if len(listi)!=0:
            filelist=listi
            nexti = i+1
            break
    OldestRestartFilePath = min(filelist, key=os.path.getctime)
    filePattern = "{:d}_{:s}".format(nexti,OriginalFilePattern)

    '''
    configName = "Input/Configurations/Config_{:s}.dat".format(filePattern, real)

    header = ["LAMMPS Description \n \n",
              "\t " + str(system.numAll) + " atoms \n \t " + str(system.numBonds) +
              " bonds \n \t " + str(system.numAngles) + " angles \n \t 0 dihedrals \n \t 0 impropers \n",
              "\n \t " + str(
                  system.numTypes) + " atom types \n \t 2 bond types \n \t 0 angle types \n \t 0 dihedral types \n \t 0 improper types \n",
              "\n \t " + str(-system.Lx * 0.5) + " " + str(system.Lx * 0.5) + " xlo xhi\n \t",
              str(-system.Ly * 0.5) + " " + str(system.Ly * 0.5) + " ylo yhi \n \t",
              str(-system.Lz * 0.5) + " " + str(system.Lz * 0.5) + " zlo zhi\n"]
    header.append("\nMasses \n \n")
    for i in range(len(system.type_mass_list)):
        header.append("\t {:d} {:.5f} \n".format(system.type_mass_list[i][0], system.type_mass_list[i][1]))

    f = open(configName, "w")
    for item in header:
        f.write("{:s} ".format(item))
    for item in system.coords:
        f.write("{:s} ".format(item))
    # for itm in system.bonds:
    #     f.write("{:s} ".format(item))
    f.close()
    '''

    # Write input script
    inputscriptfile = write_in_script(sigma, system.numTypes, PatchRange, PatchStrength, IsotropicAttrRange,
                                      IsotropicAttrStrength, real, RunSteps, dumpevery, filePattern, OldestRestartFilePath,
                                      ResultsFolder, extForce)

    # Create run file for cluster
    runfilename = writeRunScript(filePattern, ResultsFolder, inputscriptfile, MaxRunTime, MPInum)

    print(inputscriptfile)
    # clusterscriptfile = write_cluster_script(MPInum)
