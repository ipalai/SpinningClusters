# From code for fluid clusters as described in:
# I. Palaia and A. Saric, Controlling cluster size in 2D phase-separating binary mixtures with specific interactions, J. Chem. Phys. (2022)
#
import argparse
from make_initconfig import make_particles
import os
import numpy as np


def write_in_script(q_A, sigma, PatchRadius, PatchStrength, IsotropicAttrRange, IsotropicAttrStrength, real, RunSteps, dumpevery, filePattern, configName, ResultsFolder):

    seed = 100 + real
    filename = "Input/Scripts/Input_{:s}_{:d}.in".format(filePattern, real)
    Temperature = 1.0
    LangevinDamping = 0.1


    f = open(filename, "w")

    f.write("log                    {:s}/Log_{:s}_{:d}.dat \n\n".format(ResultsFolder,filePattern,real))

    f.write("units                  lj \n")
    f.write("dimension              2 \n")
    f.write("atom_style             full \n")
    f.write("boundary               p p p \n")
    f.write("read_data              {:s} \n\n".format(configName))

    f.write("group                  Central type 1")
    f.write("group                  Patch type 2")

    f.write("neighbor               0.5 bin\n")
    f.write("neigh_modify           every 1 delay 1\n")
    f.write("neigh_modify           exclude molecule/intra all \n\n")

    GlobalCutoff = sigma + IsotropicAttrRange
    f.write("pair_style             cosine/squared {:.4f} \n".format(GlobalCutoff))
    f.write("pair_coeff             * * cosine/squared 0 1.0 1.0 wca \n")

    f.write("pair_coeff             1 1  {:.1f} {:.3f} {:.3f} wca \n".format(IsotropicAttrStrength, sigma, sigma + IsotropicAttrRange))

    for id1 in range(2,q_A+2):
        for id2 in range(id1,q_A+2):
            f.write("pair_coeff             {:d} {:d}  {:.1f} {:.3f} {:.3f} \n".format(id1, id2, PatchStrength, 0.0, PatchRadius))  # mirrored cosine squared if PatchStrength<0 (purely repulsive)



    f.write("\nvelocity               all create {:.3f} {:d} \n".format(Temperature,seed))

    f.write("thermo                 {:d}\n".format(dumpevery))
    f.write("thermo_style           custom step temp press etotal epair\n")
    f.write("thermo_modify          flush yes\n")
    f.write("timestep               0.01\n\n")

    f.write("compute                cmol all chunk/atom molecule \n")
    f.write("compute                cCMmol all com/chunk cmol \n")
    f.write("variable               vCMx atom c_cCMmol[mol][1] \n")
    f.write("variable               vCMy atom c_cCMmol[mol][2] \n")
    f.write("variable               fxPatch   atom  {:f}*(y-v_vCMy)/{:f} \n".format(extForce,PatchRadialDistance))
    f.write("variable               fyPatch   atom -{:f}*(x-v_vCMx)/{:f} \n".format(extForce,PatchRadialDistance))
    f.write("variable               fxCentral atom -{:f}*(y-v_vCMy)/{:f} \n".format(extForce,PatchRadialDistance))
    f.write("variable               fyCentral atom  {:f}*(x-v_vCMx)/{:f} \n\n".format(extForce,PatchRadialDistance))

    f.write("fix                    fLANG all langevin {:.3f} {:.3f} {:.3f} {:d} \n".format(Temperature, Temperature, LangevinDamping, seed))
    # f.write("fix                    fNVE all nve\n")
    f.write("fix                    fRigidNVE all rigid/nve molecule\n")
    f.write("fix                    fEnforce2d all enforce2d\n")
    f.write("fix                    fTorqueCentral Central v_fxCentral v_fyCentral 0 \n")
    f.write("fix                    fTorquePatch Patch v_fxPatch v_fyPatch 0 \n\n")

    f.write("dump                   1 all custom {:d} {:s}/Traj_{:s}_{:d}.xyz id type mol x y z \n\n".format(dumpevery, ResultsFolder, filePattern, real))

    f.write("run                    {:d}\n".format(RunSteps))

    f.close()

    return filename







if __name__ == "__main__":

    dumpevery = 5e5
    EqSteps = 1e5
    RunSteps = 5e7

    ## Simulation parameters

    Energy_Patches = 10.0
    Energy_A = 10.0
    
    num_A = 200
    num_tot = num_A
    density_A = 0.25 # number of colloids per sigma^-2. A density of 0.25 corresponds to a packing fraction of 0.2

    qA = 6

    sigma = 1.0
    IsotropicAttrRange = 0.2
    IsotropicAttrStrength = 10.0

    PatchRadius = 0.1
    PatchRadialDistance = 0.5
    PatchStrength = -IsotropicAttrStrength

    Temperature = 1

    real = 0


    parser = argparse.ArgumentParser(description="Script for spinning, aggregating colloids.")
    #parser.add_argument('--eps_pp', '-eps_pp', dest='eps_pp', action='store', type=float, default=eps_pp, help='Energy depth of repulsion between proteins')
    #parser.add_argument('--eps_ll', '-eps_ll', dest='eps_ll', action='store', type=float, default=eps_ll, help='Energy depth of attraction/repulsion between cross patchs')
    parser.add_argument('--numA', '-nA', dest='num_A', action='store', type=int, default=num_A, help='Number of aggregating colloids.')
    parser.add_argument('--density_A', '-density_A', '-dens', dest='density_A', action='store', type=float, default=density_A, help='Surface density of colloids in units of sigma^-2 (multiply by 0.25*pi*sigma^2 to obtain packing fraction).')
    parser.add_argument('--q_A', '-qA', dest='qA', action='store', type=int, default=q_A, help='Number of patches on 1 colloid.')
    #parser.add_argument('--patch_sigma', '-ps', dest='patch_sigma', action='store', type=float, default=patch_sigma, help='Diameter of patch')
    parser.add_argument('--PatchRadius', '-pr', '-rp', dest='PatchRadius', action='store', type=float, default=PatchRadius, help='Radius of patches.')
    parser.add_argument('--PatchRadialDistance', '-pd', '-dp' dest='PatchRadialDistance', action='store', type=float, default=PatchRadialDistance, help='Radial distance of patches from centre of colloid.')
    parser.add_argument('--PatchStrength', '-ps', '-ep', dest='PatchStrength', action='store', type=float, default=PatchStrength, help='Radial distance of patches from centre of colloid.')
    parser.add_argument('--IsotropicAttrRange', '-ra', dest='IsotropicAttrRange', action='store', type=float, default=IsotropicAttrRange, help='Attraction range of the cosine-squared isotropic attraction between molecules.')
    parser.add_argument('--IsotropicAttrStrength', '-ea', dest='IsotropicAttrStrength', action='store', type=float, default=IsotropicAttrStrength, help='Strength of the cosine-squared isotropic attraction between molecules.')
    parser.add_argument('--real','-real', dest='real', action='store', type=int, default=0, help='Number of realisation (statistics), sets the seed.')
    parser.add_argument('--Temperature','-T', dest='Temperature', action='store', type=float, default=0, help='Temperature.')
    parser.add_argument('-runsteps','--RunSteps', dest='RunSteps', action='store', type=int, default=RunSteps, help='Number of time steps.')
    parser.add_argument('-dumpevery','--dumpevery', dest='dumpevery', action='store', type=int, default=dumpevery, help='Number of skipped time steps in dump file.')

    args = parser.parse_args()

    num_A = args.num_A
    density_A = args.density_A
    q_A = args.q_A
    PatchRadius = args.PatchRadius
    PatchRadialDistance = args.PatchRadialDistance
    PatchStrength = args.PatchStrength
    IsotropicAttrRange = args.IsotropicAttrRange
    IsotropicAttrStrength = args.IsotropicAttrStrength
    real = args.real
    RunSteps = args.RunSteps
    dumpevery = args.dumpevery

    GlobalCutoff = sigma + IsotropicAttrRange


    ## Initial configuration

    side = np.ceil(np.sqrt(num_tot))
    lattice_factor = np.sqrt(1.0*num_A/(side*side)/density_A)       ## Ensures that A surface density is exactly 0.0075 (protein_sigma)^-2
    densitycheck = num_A/(side*lattice_factor)**2                         ## In units of (protein_sigma)^-2
    assert np.isclose(density_A,densitycheck), "ERROR: Densities do not correspond."

    system = make_particles(sigma, PatchRadius, PatchRadialDistance, num_A, q_A, side, lattice_factor, real)
    system.make_A()


    filePattern = "qA{:d}_dp{:.2f}_dens{:.3f}_nA{:d}_rp{:.2f}_ra{:.2f}_ep{:.1f}_ea{:.1f}_T{:.1f}".format(q_A,PatchRadialDistance,density_A,num_A,PatchRadius,IsotropicAttrRange,PatchStrength,IsotropicAttrStrength,Temperature)
    ResultsFolder = "Results_qA{:d}_dp{:.2f}_dens{:.3f}".format(q_A,PatchRadialDistance,density_A)

    if not os.path.exists(ResultsFolder):
        os.makedirs(ResultsFolder)

    configName = "Config_{:s}_{:d}.dat".format(filePattern,real)

    header = ["LAMMPS Description \n \n",
              "\t " + str(system.numAll) + " atoms \n \t " + str(system.numBonds) +
              " bonds \n \t " + str(system.numAngles) + " angles \n \t 0 dihedrals \n \t 0 impropers \n",
              "\n \t "+str(system.numTypes)+" atom types \n \t 2 bond types \n \t 0 angle types \n \t 0 dihedral types \n \t 0 improper types \n",
              "\n \t " + str(-system.Lx*0.5) + " " + str(system.Lx*0.5) + " xlo xhi\n \t", str(-system.Ly*0.5) + " " + str(system.Ly*0.5) + " ylo yhi \n \t",
              str(-system.Lz*0.5) + " " + str(system.Lz*0.5) + " zlo zhi\n"]

    header.append("\nMasses \n \n")
    for i in range(len(system.type_mass_list)):
        header.append("\t {:d} {:.5f} \n".format(system.type_mass_list[i][0],system.type_mass_list[i][1]))

    f = open(configName, "w")

    for item in header:
        f.write("{:s} ".format(item))

    for item in system.coords:
        f.write("{:s} ".format(item))


    f.close()

    ## Write input script
    inputscriptfile = write_in_script(q_A, sigma, PatchRadius, PatchStrength, IsotropicAttrRange, IsotropicAttrStrength, real, RunSteps, dumpevery, filePattern, configName, ResultsFolder)
    clusterscriptfile = write_cluster_script(MPInum)