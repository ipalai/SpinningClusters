# From code for fluid clusters as described in:
# I. Palaia and A. Saric, Controlling cluster size in 2D phase-separating binary mixtures with specific interactions, J. Chem. Phys. (2022)
#
import argparse
import os
import numpy as np
import subprocess


def write_in_script(sigma, numParticleTypes, PatchRange, PatchStrength, IsotropicAttrRange, IsotropicAttrStrength, real,
                    RunSteps, dumpevery, filePattern, configName, ResultsFolder, extForce):
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

    f.write("units                  lj \n")
    f.write("dimension              2 \n")
    f.write("atom_style             full \n")
    f.write("boundary               p p p \n")
    f.write("read_data              {:s} \n\n".format(configName))

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

    f.write("\nvelocity               all create {:.3f} {:d} \n\n".format(Temperature, seed))
    # f.write("\nvelocity               all create {:.3f} {:d} \n".format(0, seed))

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

    f.write("run                    {:d}\n".format(RunSteps))

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
    RunSteps = 1000000
    MaxRunTime = 48
    MPInum = 1

    # Simulation parameters
    q_A = 5
    sigma = 1.0

    IsotropicAttrRange = 0.15
    IsotropicAttrStrength = 20.0

    PatchRange = 0.15
    PatchRadialDistance = 0.5
    PatchStrength = -IsotropicAttrStrength

    Temperature = 1

    extTorque = 10.0

    real = 0

    parser = argparse.ArgumentParser(description="Script for spinning, aggregating colloids.")

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
    parser.add_argument('--configfile', dest='configfile', action='store', type=str, default='',
                        help="Initial configuration file, containing cluster extracted from previous simulation.")

    args = parser.parse_args()

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
    configfile = args.configfile

    extForce = extTorque / PatchRadialDistance

    # Initial configuration
    filePattern = subprocess.check_output(
        " echo {:s} | sed 's/.*\/Config_.*\/Config_//' | sed 's/.dat//'".format(configfile), shell=True)
    filePattern = filePattern.decode("utf-8")[:-1]
    folderPattern = subprocess.check_output(" echo {:s} | sed 's/.*ConfigurationsToRun\/Config_//' | sed  's/\/Config.*//' ".format(configfile), shell=True)
    folderPattern = folderPattern.decode("utf-8")[:-1]
    ResultsFolder = "Results_{:s}".format(folderPattern)

    if not os.path.exists(ResultsFolder):
        os.makedirs(ResultsFolder)


    # Write input script
    inputscriptfile = write_in_script(sigma, 3, PatchRange, PatchStrength, IsotropicAttrRange,
                                      IsotropicAttrStrength, real, RunSteps, dumpevery, filePattern, configfile,
                                      ResultsFolder, extForce)

    # Create run file for cluster
    runfilename = writeRunScript(filePattern, ResultsFolder, inputscriptfile, MaxRunTime, MPInum)

    print(inputscriptfile)
    # clusterscriptfile = write_cluster_script(MPInum)
