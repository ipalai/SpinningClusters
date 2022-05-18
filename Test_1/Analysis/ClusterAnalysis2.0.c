// 2.0 Adapted from TheoryClustersFluidljcossq's ClusterAnalysis1.2.c as a clustering code for SpinningClusters

// 1.1 Branch from ClusterAnalysis1.1.c, to account for fluid clusters.
//     The adjacency matrix (and therefore clusters) are now based on physical distance only (not bonds); I keep counting bonds for statistics though
// 1.2 NMAX changed from 10000 to 12000 for large-nB simulations (8Jun2021)
//     dist 2d-array is now a static int with length NMOLECULESMAX, to avoid Segmentation Fault due to stack overflow.

// PRE-EXISTING LOG FROM ClusterAnalysis1.1.c
// 0.2 Verifies preliminarily if atoms interact according to InteractionsMatrix, before computing distances in adjacencyfunction.
//     This avoids considering clustered two particles that are just close to each other, but not interacting.
// 0.3 Computes AvgPLCConnectivity and CoalescenceLikelihood1
// 0.4 Classifies both existing and theoretically possible bond, according to their type (as per the InteractionsMatrix)
// 0.5 Computes moment of inertia
// 0.6 Computes graph-theoretical moment of inertia and small-world network compactness.
//     Corrected ColescenceLikelihood (last molecule was not considered)
// 1.0 Adapted from 0.6 for TheoryClusters project
// 1.1 Accepts output file name as a 3rd argument from command line


// USAGE: ./$executable $xyzInputFile $Timestep $OutputFile

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#define NMAX 12000
#define NMOLECULESMAX 2000
#define NATOMTYPESMAX 14
#define MOLECULETYPES 1
#define VALENCEMAX 6
#define NINTERACTIONSMAX 36
#define strlength 1500
#define SQR(x) ((x)*(x))
#define PI 3.14159265359


float Position[NMAX][2], Velocity[NMAX][2], CutoffDistanceBond, CutoffDistanceBondSq, CutoffDistanceDist, CutoffDistanceDistSq, BoxLength;
int NumberOfInteractions;
int NumberOfMolecules;
int AtomType[NMAX];
int ClusterID[NMOLECULESMAX];
static int AdjacencyMatrixDist[NMOLECULESMAX][NMOLECULESMAX];
static int AdjacencyMatrixBond[NMOLECULESMAX][NMOLECULESMAX];
int ClusterSize[NMOLECULESMAX];
int OrderedClustersIndex[NMOLECULESMAX];
int ClustersCounter;
long long int gtinertia;
int gtdiameter;
float gtMooreBound;
int gtNodesByDegree[VALENCEMAX+1]; //First column (0) and first row (#edges=0) not used!
int qA;


float fminfunction(float a, float b)
{
	if(a<b)
		return a;
	else
		return b;
}



int adjacencydistfunction(int i, int j)
{
    float dx, dy;
    dx = fabsf(Position[i][0]-Position[j][0]);
    dx = fminfunction(dx, BoxLength-dx);
    if(dx>CutoffDistanceDist)
        return 0;
    dy = fabsf(Position[i][1]-Position[j][1]);
    dy = fminfunction(dy, BoxLength-dy);
    if(dy>CutoffDistanceDist)
        return 0;
    if(SQR(dx)+SQR(dy)<CutoffDistanceDistSq)
        return 1;
    else
        return 0;
}




/*
float unwrappeddistance2function(int i, int j, int WindingNumberjx, int WindingNumberjy)
{
    dx = (Position[j][0] + BoxLength*WindingNumberjx) - Position[i][0];
    dy = (Position[j][1] + BoxLength*WindingNumberjy) - Position[i][1];
    return (SQR(dx)+SQR(dy));
}
*/


void sortclustersfunction()
{
    int i, j, k;
    int SkipFlag=0;
    int x=0;
    for(k=0; k<ClustersCounter; k++)
        OrderedClustersIndex[k]=-1;
    for(i=0; i<ClustersCounter; i++)
    {
        x=0;
        for(j=0; j<ClustersCounter; j++)
        {
            // check if index j has already been selected previously. If so, skip j
            SkipFlag=0;
            for(k=0; k<i; k++)
            {
                if(j==OrderedClustersIndex[k])
                    SkipFlag=1;
            }
            // Check if ClusterSize[j][0]>x
            if(SkipFlag==0 && ClusterSize[j]>x)
            {
                OrderedClustersIndex[i]=j;
                x=ClusterSize[j];
            }
        }
    }
}






////////////////////////////////
////////////////////////////////

int main(int argc, char *argv[])
{
	float x, y, vx, vy, AvgClusterSizeAll, AvgClusterSizeA, AvgClusterSizeB, DistrClusterSize[NMOLECULESMAX], AConnectivity, BConnectivity;
	float AvgClusterSizeAllBigger, AvgClusterSizeABigger, AvgClusterSizeBBigger;
	int ClustersCounterBigger;
	int i, j, k, ki, kj, timestep, FlagRead, FirstCharToCopy;
	int NumberOfAClusters, NumberOfBClusters, NumberOfAMolecules, NumberOfBMolecules, CheckCounter;
	int NumberOfAtoms, AtomID, MoleculeIDDummy, MoleculeID[NMAX], AtomTypeDummy;
	int iToDo[NMOLECULESMAX], iToDoCounter, iToDoSegno, Sum, iDone[NMOLECULESMAX];
	int CentralAtom[NMOLECULESMAX];
	int WindingNumber[NMOLECULESMAX][2];
	int ci, cj, mi;
	long long int BondCounter[NINTERACTIONSMAX];
	long long int CoalescenceLikelihood1;
    double CenterOfMass[NMOLECULESMAX][2];
    double Inertia[NMOLECULESMAX][6];   // Ixx Iyy Ixy Izz lambda1 lambda2
    double AngularMomentumz[NMOLECULESMAX];  // Angular momentum of clusters along z
    double Omegaz[NMOLECULESMAX];  // Angular velocity of clusters along z
    float XVyMinusYVx, mass;
    float SigmaLJ, ra, PatchRadialDistance;
    float MassAlpha, MassCentralAtom, MassPatch;
    float BoxEdgeForUnwinding;
    double UnwoundPosition[2];
    double trace, determinant, Sphericity;
    int EdgesInCluster[NMOLECULESMAX];
    int NumberOfTerminalNodes[NMOLECULESMAX];
    int NumberOfUnboundNodes[NMOLECULESMAX];
    int NumberOfEdgesForThisMolecule;
	char string[strlength], FileNameOutput[strlength], FileNameTraj[strlength], timestepString[strlength];
	char *p;
	char *pch;
	FILE *FileOutput;
	FILE *FileTraj;

	printf("Initialization...\n");

	// Get arguments
    strcpy(FileNameTraj,argv[1]);
    sprintf(timestepString, "%s\n", argv[2]);
    sscanf(timestepString, "%d", &timestep);
    timestep = strtol(argv[2], &p, 10);
    strcpy(FileNameOutput, argv[3]);
    // Get qA
    pch = strpbrk (FileNameTraj, "q");
    while (pch != NULL)
    {
        if(sscanf(pch, "qA%d_%*s", &qA)==1)
            break;
        pch = strpbrk (pch+1,"q");
    }
    //Get CutoffDistanceDist = range of isotropic attraction
    pch = strpbrk (FileNameTraj, "_");
    while (pch != NULL)
    {
        if(sscanf(pch, "_ra%f_%*s", &ra)==1)
            break;
        pch = strpbrk (pch+1,"_");
    }
    SigmaLJ=1.0;
    CutoffDistanceDist=ra+SigmaLJ;
    CutoffDistanceDistSq=SQR(CutoffDistanceDist);
    //Get PatchRadialDistance to compute masses
    pch = strpbrk (FileNameTraj, "_");
    while (pch != NULL)
    {
        if(sscanf(pch, "_dp%f_%*s", &PatchRadialDistance)==1)
            break;
        pch = strpbrk (pch+1,"_");
    }
    MassAlpha = 0.1*SQR(SigmaLJ/PatchRadialDistance);  // fraction of mass due to patches such that moment of inertia equals that of a sphere
    MassCentralAtom = 1.0-MassAlpha;
    MassPatch = MassAlpha/qA;
    // Summary
    printf("qA %d, ra %f, CutoffDistanceDist %f, PatchRadialDistance %f. Masses %f %f\n", qA, ra, CutoffDistanceDist, PatchRadialDistance, MassCentralAtom, MassPatch);



    // Initialization
	for(i=0; i<NMOLECULESMAX; i++)
	{
        for(j=0;j<2;j++)
        {
            WindingNumber[i][j]=0;
            CenterOfMass[i][j]=0;
        }
        for(j=0; j<6; j++)
            Inertia[i][j]=0;
        AngularMomentumz[i]=0;
        Omegaz[i]=0;
        EdgesInCluster[i]=0;
        NumberOfTerminalNodes[i]=0;
	}
    AConnectivity=0;


    // Read file
    FileTraj=fopen(FileNameTraj,"r");
  	printf("Reading file %s...\n", FileNameTraj);

    // Find desired TIMESTEP in Traj file
    FlagRead=0;
	while(FlagRead==0)
    {
        fgets(string, strlength, FileTraj);
        if(strcmp(string, "ITEM: TIMESTEP\n")==0)
        {
            fgets(string, strlength, FileTraj);
            if(strcmp(string,timestepString)==0)
                FlagRead=1;
            //printf("%s", string);
        }
        if(feof(FileTraj)!=0)
        {
            printf("*** ERROR: FileTraj does not contain TIMESTEP %s", timestepString);
	    return 1;
        }
    }
    fgets(string, strlength, FileTraj);
    fgets(string, strlength, FileTraj);
    NumberOfAtoms = strtol(string, &p, 10);

    fgets(string, strlength, FileTraj);
    fgets(string, strlength, FileTraj);
    sscanf(string, "%f %f", &x, &BoxLength);
    BoxLength*=2;
    //printf("BoxLength = %f\n",BoxLength);

    // Read atoms' positions
    FlagRead=0;
	while(FlagRead==0)
    {
        fgets(string, strlength, FileTraj);
        if(strncmp(string, "ITEM: ATOMS", strlen("ITEM: ATOMS"))==0)
                FlagRead=1;
    }

    NumberOfMolecules=0;
    NumberOfAMolecules=0;
    for(i=0; i<NumberOfAtoms; i++)
    {
        fgets(string, strlength, FileTraj);
        sscanf(string, "%d %d %d %f %f %*f %f %f %*f", &AtomID, &AtomTypeDummy, &MoleculeIDDummy, &x, &y, &vx, &vy );
        AtomType[AtomID]=AtomTypeDummy;
        MoleculeID[AtomID]=MoleculeIDDummy;
        Position[AtomID][0]=x;
        Position[AtomID][1]=y;
        Velocity[AtomID][0]=vx;
        Velocity[AtomID][1]=vy;
        if(AtomTypeDummy==1)
        {
            NumberOfAMolecules++;
            CentralAtom[MoleculeIDDummy]=AtomID;
        }
        if(MoleculeIDDummy>NumberOfMolecules)
            NumberOfMolecules=MoleculeIDDummy;

    }


    fclose(FileTraj);



	// *****************************************************************************************
	// ************ Building adjacency matrix (based on distance between central atoms) ********

	printf("Building adjacency matrix...");

    for(i=1; i<=NumberOfMolecules; i++)
    {
        //printf("\n%d: ", i);
        AdjacencyMatrixDist[i][i]=0;
        for(j=i+1; j<=NumberOfMolecules; j++)
        {
            Sum=adjacencydistfunction(CentralAtom[i], CentralAtom[j]);
            AdjacencyMatrixDist[i][j]=Sum;
            AdjacencyMatrixDist[j][i]=Sum;
            //if(Sum!=0)
                //printf("%d %d, ", i,j);
        }
    }


	// **********************************************************************
    // *************** Compute clusters *************************************

	printf("\nAnalyzing adjacency matrices...\n");

    ClustersCounter=0;
    iToDoCounter=0;
    iToDo[0]=-1;
    iToDoSegno=0;
    for(i=0; i<NMOLECULESMAX; i++)
        iDone[i]=-1; // Tells whether molecule i has been analysed or not:  -1 = not done nor in the to-do list; 0 = in the to-do list; 1 = done
    i=1;
    CheckCounter=0;
    BoxEdgeForUnwinding=0.5*BoxLength-SigmaLJ-CutoffDistanceDist;
    //printf("Box %f %f\n", BoxLength, BoxEdgeForUnwinding);

	while(i<=NumberOfMolecules)
    {
        CheckCounter++;
        ClusterID[i]=ClustersCounter;
        iDone[i]=1;
        NumberOfEdgesForThisMolecule=0;
        for(j=1; j<=NumberOfMolecules; j++)
        {
            if(AdjacencyMatrixDist[i][j]>0)
            {
                if(iDone[j]==-1)
                {
                    iToDo[iToDoCounter]=j;
                    iToDoCounter++;
                    iToDo[iToDoCounter]=-1;
                    iDone[j]=0;

                    ci=CentralAtom[i];              // WindingNumber for moment of inertia
                    cj=CentralAtom[j];
                    for(k=0; k<2; k++)
                    {
                        if( Position[ci][k]>BoxEdgeForUnwinding && Position[cj][k]<-BoxEdgeForUnwinding )
                            WindingNumber[j][k]=WindingNumber[i][k]+1;
                        else if( Position[ci][k]<-BoxEdgeForUnwinding && Position[cj][k]>BoxEdgeForUnwinding )
                            WindingNumber[j][k]=WindingNumber[i][k]-1;
                        else
                            WindingNumber[j][k]=WindingNumber[i][k];
                    }
                }
                NumberOfEdgesForThisMolecule++;
                AConnectivity+=AdjacencyMatrixBond[i][j];
            }
        }

        EdgesInCluster[ClustersCounter]+=NumberOfEdgesForThisMolecule;  // later divided by 2, as I'm counting twice every edge
        if(NumberOfEdgesForThisMolecule==1)
            NumberOfTerminalNodes[ClustersCounter]++;

        if(iToDo[iToDoSegno]==-1)   // (Cluster is complete)
        {
            if(EdgesInCluster[ClustersCounter]%2==0)        // Divide by 2 EdgesInClusters, because I counted twice every edge
                EdgesInCluster[ClustersCounter]/=2;
            else
            {
                printf("*** ERROR: EdgesInCluster[%d]=%d is not even!", ClustersCounter, EdgesInCluster[ClustersCounter]);
                //printf("Molecules in error cluster: ");
                //for(k=0;k<=NumberOfMolecules;k++)
                //    if(ClusterID[k]==ClustersCounter)
                //        printf("%d ", k);
                return 1;
            }

            j=1;
            while(iDone[j]==1)      // find first not-done molecule (number j): the next cluster will start from here
                j++;
            i=j;
            ClustersCounter++;
        }
        else
        {
            i=iToDo[iToDoSegno];
            iToDoSegno++;
        }
    }

    AConnectivity/=NumberOfAMolecules;

    if(CheckCounter!=NumberOfMolecules)
    {
        printf("*** ERROR: CheckCounter!=NumberOfMolecules\n");
        return 1;
    }


    // ************************************************************************************************
    // **************** Compute Cluster sizes, size distribution, center of mass, inertia tensor ******

    // Cluster Size distribution is currently Molecules distribution according to size of the cluster they are in

    for(k=0; k<NMOLECULESMAX; k++)
    {
        ClusterSize[k]=0;
        DistrClusterSize[k]=0;
    }
    NumberOfAClusters=0;

    for(i=1; i<=NumberOfMolecules; i++)
    {
        ClusterSize[ClusterID[i]]++;

        // Inertia
        ci=CentralAtom[i];
        for(k=0; k<2; k++)
        {
            UnwoundPosition[k]= Position[ci][k] + BoxLength*WindingNumber[i][k];
            CenterOfMass[ClusterID[i]][k] += UnwoundPosition[k];
        }
        Inertia[ClusterID[i]][0] += SQR(UnwoundPosition[1]);  // I_xx
        Inertia[ClusterID[i]][1] += SQR(UnwoundPosition[0]);  // I_yy
        Inertia[ClusterID[i]][2] -= UnwoundPosition[0]*UnwoundPosition[1];  // I_xy
    
    }

    AvgClusterSizeAll=0;
    AvgClusterSizeAllBigger=0;
    ClustersCounterBigger=0;

    for(k=0; k<ClustersCounter; k++)
    {
        AvgClusterSizeAll+=ClusterSize[k];
        if(ClusterSize[k]>=2) // Threshold for considering a clusters a BiggerCluster
        {
            AvgClusterSizeAllBigger+=ClusterSize[k];
            ClustersCounterBigger+=1;
        }
        DistrClusterSize[ClusterSize[k]]+=ClusterSize[k]*1.0/NumberOfMolecules;

        // Approximate moment of gyration and center of mass
        for(j=0; j<2; j++)
            CenterOfMass[k][j]/=ClusterSize[k];
        Inertia[k][0] -= ClusterSize[k]*SQR(CenterOfMass[k][1]); // I_xx
        Inertia[k][1] -= ClusterSize[k]*SQR(CenterOfMass[k][0]); // I_yy
        Inertia[k][2] += ClusterSize[k]*CenterOfMass[k][0]*CenterOfMass[k][1]; // I_xy
        Inertia[k][3] = Inertia[k][0]+Inertia[k][1];    // I_zz
        //printf("%d %g %g %g %g %g %g %g\n", ClusterSize[k][0], CenterOfMass[k][0], CenterOfMass[k][1], Inertia[k][0], Inertia[k][1], Inertia[k][2], Inertia[k][3], sqrt(3.0)/4/PI*SQR(SigmaLJ)*SQR(ClusterSize[k][0])/Inertia[k][3] );
        trace = Inertia[k][0]+Inertia[k][1];
        determinant= Inertia[k][0]*Inertia[k][1] - SQR(Inertia[k][2]);
        Inertia[k][4] = 0.5*(trace-sqrt(SQR(trace)-4*determinant));  // first principal component of the moment of gyration
        Inertia[k][5] = 0.5*(trace+sqrt(SQR(trace)-4*determinant));  // second principal component of the moment of gyration

        // Angular momentum and angular velocity (Omega)
        AngularMomentumz[k]=0;
        Omegaz[k]=0;
        for(AtomID=1; AtomID<=NumberOfAtoms; i++)
        {
            mi=MoleculeID[AtomID];
            if(ClusterID[mi]==k)
            {
                if(CentralAtom[mi]==AtomID)
                {
                    mass=MassCentralAtom;
                    if(ClusterSize[ClusterID[mi]]==1)
                        continue;
                }
                else
                    mass=MassPatch;
                for(k=0; k<2; k++)
                    UnwoundPosition[k]= Position[AtomID][k] + BoxLength*WindingNumber[mi][k];
                XVyMinusYVx= (UnwoundPosition[0]-CenterOfMass[k][0])*Velocity[AtomID][1] - (UnwoundPosition[1]-CenterOfMass[k][1])*Velocity[AtomID][0];
                AngularMomentumz[k]+= mass * XVyMinusYVx;
                Omegaz[k]+= XVyMinusYVx/(SQR(UnwoundPosition[0])+SQR(UnwoundPosition[1]));  // r_i cross v_i / r_i^2 
            }
        }
        

    }
    AvgClusterSizeAll/=ClustersCounter;

    if(ClustersCounterBigger!=0)
        AvgClusterSizeAllBigger/=ClustersCounterBigger;
    sortclustersfunction();



  	// **********************************************************************
    // ************** Print output ******************************************

    printf("Printing output file...\n");
    printf("%s\n", FileNameOutput);


    FileOutput=fopen(FileNameOutput,"wb+");


    fprintf(FileOutput, "InputFile\n%s\nData analysis program\nClusterAnalysisFluidljcossq1.2.exe\nts\n%d\n\n", FileNameTraj,timestep);
    fprintf(FileOutput, "NumberOfClusters (number of clusters, irrespective of the type of molecule they are composed of)\n%d\n", ClustersCounter);
    fprintf(FileOutput, "NumberOfClustersBigger (number of clusters containing at least 2 A and 2 B molecules)\n%d\n", ClustersCounterBigger);
    fprintf(FileOutput, "\nAvgClusterSize (total number of molecules / NumberOfClusters)\n%f\n", AvgClusterSizeAll);
    fprintf(FileOutput, "\nAvgClusterSizeBigger (total number of molecules in clusters with at least 2 A and 2 B / number of such clusters)\n%f\n", AvgClusterSizeAllBigger);
    fprintf(FileOutput, "\nAvgAConnectivity (average number of molecules bound to one A molecule)\n%f\n", AConnectivity);

    fprintf(FileOutput, "\nClusters detail\n");
    fprintf(FileOutput, "ClusterSize[], EdgesInCluster, NumberOfTerminalNodes, CenterOfMass_x, CenterOfMass_y, Inertia_lambda1, Inertia_lambda2, Inertia_zz, Compactness, Roundness\n"); // Roundeness = Sphericity
 
    for(k=0; k<ClustersCounter; k++)
    {
        i=OrderedClustersIndex[k];
        Sphericity=Inertia[i][4]/Inertia[i][5];
        fprintf(FileOutput, "%d %d %d %d %d %g %g %g %g %g %g %g\n", ClusterSize[i], EdgesInCluster[i], NumberOfTerminalNodes[i], CenterOfMass[i][0], CenterOfMass[i][1], Inertia[i][4], Inertia[i][5], Inertia[i][3],  sqrt(3.0)/8/PI*SQR(SigmaLJ)*SQR(ClusterSize[i])*(1+Sphericity)/sqrt(Sphericity)/Inertia[i][3],  Sphericity );
    }

    fprintf(FileOutput, "\nDistrClusterSizeAll (fraction of particles in clusters, by size of the clusters)\n");
    for(k=0; k<=NumberOfMolecules; k++)
        fprintf(FileOutput, "%d %f %f %f\n", k, DistrClusterSize[k]);

    fclose(FileOutput);


	return 0;
}



