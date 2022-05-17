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
#define MOLECULETYPES 2
#define VALENCEMAX 6
#define NINTERACTIONSMAX 36
#define strlength 1500
#define SQR(x) ((x)*(x))
#define PI 3.14159265359


float Position[NMAX][2], CutoffDistanceBond, CutoffDistanceBondSq, CutoffDistanceDist, CutoffDistanceDistSq, BoxLength;
int NumberOfInteractions;
int NumberOfMolecules;
int AtomType[NMAX];
int MoleculeType[NMOLECULESMAX];
int ClusterID[NMOLECULESMAX];
int FreeLinkersNumberInCluster[NMOLECULESMAX][NATOMTYPESMAX+1];
static int AdjacencyMatrixDist[NMOLECULESMAX][NMOLECULESMAX];
static int AdjacencyMatrixBond[NMOLECULESMAX][NMOLECULESMAX];
int ClusterSize[NMOLECULESMAX][MOLECULETYPES+1];
int OrderedClustersIndex[NMOLECULESMAX];
long long int CoalescenceBondCounter[NINTERACTIONSMAX];
int ClustersCounter;
long long int gtinertia;
int gtdiameter;
float gtMooreBound;
int gtNodesByDegree[MOLECULETYPES+1][VALENCEMAX+1]; //First column (moleculeType=0) and first row (#edges=0) not used!
int qA, qB;
int InteractionsMatrix[VALENCEMAX*VALENCEMAX][2];


float fminfunction(float a, float b)
{
	if(a<b)
		return a;
	else
		return b;
}


int interactionflagfunction(int a, int b)
{
    int k, InteractionFlag;
    if(a>b)
    {
        k=a;
        a=b;
        b=k;
    }
    InteractionFlag=-1;
    for(k=0; k<NumberOfInteractions; k++)
    {
        if(a==InteractionsMatrix[k][0])
        {
            if(b==InteractionsMatrix[k][1])
                InteractionFlag=k;
        }
        if(a<InteractionsMatrix[k][0])
            k=NumberOfInteractions;
    }
    return InteractionFlag;
}


int adjacencybondfunction(int i, int j)
{
    float dx, dy;
    if(interactionflagfunction(AtomType[i],AtomType[j])==-1)
        return 0;
    dx = fabsf(Position[i][0]-Position[j][0]);
    dx = fminfunction(dx, BoxLength-dx);
    if(dx>CutoffDistanceBond)
        return 0;
    dy = fabsf(Position[i][1]-Position[j][1]);
    dy = fminfunction(dy, BoxLength-dy);
    if(dy>CutoffDistanceBond)
        return 0;
    if(SQR(dx)+SQR(dy)<CutoffDistanceBondSq)
        return 1;
    else
        return 0;
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

int moleculetypefunction(int AtomTypeDummy)
{
    if(AtomTypeDummy>=1 && AtomTypeDummy<=qA+1)
        return 1; // A
    if(AtomTypeDummy>=qA+2 && AtomTypeDummy<=qA+qB+2)
        return 2; // B
    else
        return -1; // error
}


long long int coalescence1function(int iCluster, int jCluster, int FlagCoalescenceBondCounter)
{
    int k, sum, dummy1, dummy2;
    sum=0;
    for(k=0; k<NumberOfInteractions; k++)
    {
        dummy1 = FreeLinkersNumberInCluster[iCluster][InteractionsMatrix[k][0]] * FreeLinkersNumberInCluster[jCluster][InteractionsMatrix[k][1]];
        dummy2 = FreeLinkersNumberInCluster[iCluster][InteractionsMatrix[k][1]] * FreeLinkersNumberInCluster[jCluster][InteractionsMatrix[k][0]];
        sum += dummy1+dummy2;
        if(FlagCoalescenceBondCounter==1)
            CoalescenceBondCounter[k] += dummy1+dummy2;
    }
    return sum;
}


float coalescence2function(int iCluster, int jCluster)
{
    return 1.0*coalescence1function(iCluster,jCluster,0)/4.0/sqrt(ClusterSize[iCluster][0]*ClusterSize[jCluster][0]);
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
            if(SkipFlag==0 && ClusterSize[j][0]>x)
            {
                OrderedClustersIndex[i]=j;
                x=ClusterSize[j][0];
            }
        }
    }
}



// Solves the all-pairs shortest path problem using Floyd Warshall algorithm
int graphtheoreticalinertiafunction (int ClusterIndex)
{
    /* dist[][] will be the output matrix that will finally have the shortest
distances between every pair of vertices */
    static int dist[NMOLECULESMAX][NMOLECULESMAX];
    int i, j, k, csize, ci, cj, imin, diam, NumberOfEdgesFromMoleculei;
    long long int Sum, MinSum;
    csize=ClusterSize[ClusterIndex][0];
    ci=-1;

    for(i=0;i<=MOLECULETYPES;i++)
    {
        for(j=0;j<=VALENCEMAX;j++)
            gtNodesByDegree[i][j]=0;
    }

    /* Restrict AdjacencyMatrix to the nodes belonging to cluster #ClusterIndex */
    for (i=1; i<=NumberOfMolecules; i++)
    {
        if(ClusterID[i]==ClusterIndex)
        {
            ci++;
            cj=-1;
            NumberOfEdgesFromMoleculei=0;
            for (j=1; j<=NumberOfMolecules; j++)
            {
                if(ClusterID[j]==ClusterIndex)
                {
                    cj++;
                    if(AdjacencyMatrixBond[i][j]>0)
                    {
                        dist[ci][cj]=1;
                        NumberOfEdgesFromMoleculei++;
                    }
                    else
                        dist[ci][cj]=csize;
                    if(i==j)
                    {
                        if(ci!=cj)
                        {
                            printf("*** ERROR: ci!=cj ");
                            return 200;
                        }
                        dist[ci][cj]=0;
                    }
                }
            }
            gtNodesByDegree[MoleculeType[i]][NumberOfEdgesFromMoleculei]++;
        }
    }
    if(ci!=csize-1 || cj!=csize-1)
    {
        printf("*** ERROR: ci!=csize-1 || cj!=csize-1   %d %d %d\n ", ci, cj, csize);
        return 201;
    }

    /* Floyd Warshall algorithm to determine matrix of shorter distances*/
    for (k=0; k<csize; k++)
    {
        // Pick all vertices as source one by one
        for (i=0; i<csize; i++)
        {
            // Pick all vertices as destination for the above picked source
            for (j=0; j<csize; j++)
            {
                // If vertex k is on the shortest path from i to j, then update the value of dist[i][j]
                if (dist[i][k] + dist[k][j] < dist[i][j])
                    dist[i][j] = dist[i][k] + dist[k][j];
            }
        }
    }

    // Computer minimum graph-theoretical moment of inertia
    MinSum=csize*csize*csize+1;
    imin=-1;
    diam=0;
    for(i=0; i<csize; i++)
    {
        Sum=0;
        for(j=0; j<csize; j++)
        {
            Sum+=SQR(dist[i][j]);
            if(dist[i][j]>diam)
                diam=dist[i][j];
        }
        if(Sum<MinSum)
        {
            MinSum=Sum;
            imin=i;
        }
    }

    // find graph-theoretical center-of-mass node (i)
    k=-1;
    for(i=1; i<=NumberOfMolecules; i++)
    {
        if(ClusterID[i]==ClusterIndex)
            k++;
        if(k==imin)
            break;
    }

    gtMooreBound=log(1+0.5*(csize-1))/log(3.0);
    gtinertia=MinSum;
    gtdiameter=diam;
    return 0;
}



int main(int argc, char *argv[])
{
	float x, y, AvgClusterSizeAll, AvgClusterSizeA, AvgClusterSizeB, DistrClusterSize[NMOLECULESMAX][3], AConnectivity, BConnectivity;
	float AvgClusterSizeAllBigger, AvgClusterSizeABigger, AvgClusterSizeBBigger;
	int ClustersCounterBigger;
	int i, j, k, ki, kj, timestep, FlagRead, FirstCharToCopy;
	int NumberOfAClusters, NumberOfBClusters, NumberOfAMolecules, NumberOfBMolecules, CheckCounter;
	int NumberOfAtoms, AtomID, MoleculeIDDummy, MoleculeID[NMAX], AtomTypeDummy;
	int Linkers[NMOLECULESMAX][VALENCEMAX], linker;
	int iToDo[NMOLECULESMAX], iToDoCounter, iToDoSegno, Sum, iDone[NMOLECULESMAX], FreeLinkersFlag[NMAX];
	int CentralAtom[NMOLECULESMAX];
	int WindingNumber[NMOLECULESMAX][2];
	int ci, cj;
	long long int BondCounter[NINTERACTIONSMAX];
	long long int CoalescenceLikelihood1;
    float CoalescenceLikelihood2;
    float FractionOfAInClustersOfSizeUpTo1, FractionOfAInClustersOfSizeUpTo5, FractionOfAInClustersOfSizeUpTo10;
    double CenterOfMass[NMOLECULESMAX][2];
    double Inertia[NMOLECULESMAX][6];   // Ixx Iyy Ixy Izz lambda1 lambda2
    float FractionOfFreeB;
    float SigmaLJ, ra;
    float BoxEdgeForUnwinding;
    double UnwoundPosition[2];
    double trace, determinant, Sphericity;
    int EdgesInCluster[NMOLECULESMAX];
    int NumberOfTerminalNodes[NMOLECULESMAX];
    int NumberOfUnboundNodes[NMOLECULESMAX];
    int NumberOfEdgesForThisMolecule;
	char string[strlength], FileNameOutput[strlength], FileNameMovie[strlength], timestepString[strlength];
	char *p;
	char *pch;
	FILE *FileOutput;
	FILE *FileMovie;

	printf("Initialization...\n");

	// Get arguments
    strcpy(FileNameMovie,argv[1]);
    sprintf(timestepString, "%s\n", argv[2]);
    sscanf(timestepString, "%d", &timestep);
    timestep = strtol(argv[2], &p, 10);
    strcpy(FileNameOutput, argv[3]);
    //printf("%d", timestep);
        // Get qA and qB
    pch = strpbrk (FileNameMovie, "q");
    while (pch != NULL)
    {
        if(sscanf(pch, "qA%d_nB%*d_qB%d_%*s", &qA, &qB)==2)
            break;
        pch = strpbrk (pch+1,"q");
    }
        //Get CutoffDistanceDist = range of isotropic attraction
    pch = strpbrk (FileNameMovie, "_");
    while (pch != NULL)
    {
        if(sscanf(pch, "_ra%f_%*s", &ra)==1)
            break;
        pch = strpbrk (pch+1,"_");
    }

    SigmaLJ=2.0;
    CutoffDistanceDist=ra+SigmaLJ;
    CutoffDistanceDistSq=SQR(CutoffDistanceDist);
    printf("qA %d, qB %d, ra %f, CutoffDistanceDist %f\n", qA, qB, ra, CutoffDistanceDist);

    CutoffDistanceBond=0.3;
   	CutoffDistanceBondSq=SQR(CutoffDistanceBond);


	// Compute InteractionsMatrix
	NumberOfInteractions=qA*qB;
	k=0;
    for(i=2; i<2+qA; i++)
    {
        for(j=3+qA; j<3+qA+qB; j++)
        {
            InteractionsMatrix[k][0]=i;
            InteractionsMatrix[k][1]=j;
            BondCounter[k]=0;
            CoalescenceBondCounter[k]=0;
            k++;
        }
    }
    if(k!=NumberOfInteractions)
    {
        printf("*** ERROR: Computation of InteractionsMatrix");
        return 1;
    }

    // Initialization
	for(i=0; i<NMOLECULESMAX; i++)
	{
        for(j=0; j<VALENCEMAX; j++)
            Linkers[i][j]=-1;       // AtomID of jth linker of molecule i
        for(j=0; j<=NATOMTYPESMAX; j++)
            FreeLinkersNumberInCluster[i][j]=0;
        for(j=0;j<2;j++)
        {
            WindingNumber[i][j]=0;
            CenterOfMass[i][j]=0;
        }
        for(j=0; j<6; j++)
            Inertia[i][j]=0;
        EdgesInCluster[i]=0;
        NumberOfTerminalNodes[i]=0;
        NumberOfUnboundNodes[i]=0;
	}
    AConnectivity=0;
    BConnectivity=0;

//    if (*p != '\0' || timestep > INT_MAX)
//    {
//        return 1;
//    }

    // Read file
    FileMovie=fopen(FileNameMovie,"r");
  	printf("Reading file %s...\n", FileNameMovie);

    // Find desired TIMESTEP in Movie file
    FlagRead=0;
	while(FlagRead==0)
    {
        fgets(string, strlength, FileMovie);
        if(strcmp(string, "ITEM: TIMESTEP\n")==0)
        {
            fgets(string, strlength, FileMovie);
            if(strcmp(string,timestepString)==0)
                FlagRead=1;
            //printf("%s", string);
        }
        if(feof(FileMovie)!=0)
        {
            printf("*** ERROR: FileMovie does not contain TIMESTEP %s", timestepString);
	    return 1;
        }
    }
    fgets(string, strlength, FileMovie);
    fgets(string, strlength, FileMovie);
    NumberOfAtoms = strtol(string, &p, 10);

    fgets(string, strlength, FileMovie);
    fgets(string, strlength, FileMovie);
    sscanf(string, "%f %f", &x, &BoxLength);
    BoxLength*=2;
    //printf("BoxLength = %f\n",BoxLength);

    // Read atoms' positions
    FlagRead=0;
	while(FlagRead==0)
    {
        fgets(string, strlength, FileMovie);
        if(strncmp(string, "ITEM: ATOMS", strlen("ITEM: ATOMS"))==0)
                FlagRead=1;
    }

    NumberOfMolecules=0;
    NumberOfAMolecules=0;
    NumberOfBMolecules=0;
    for(i=0; i<NumberOfAtoms; i++)
    {
        fgets(string, strlength, FileMovie);
        sscanf(string, "%d %d %d %f %f %*f", &AtomID, &AtomTypeDummy, &MoleculeIDDummy, &x, &y);
        AtomType[AtomID]=AtomTypeDummy;
        MoleculeID[AtomID]=MoleculeIDDummy;
        Position[AtomID][0]=x;
        Position[AtomID][1]=y;
        MoleculeType[MoleculeIDDummy]=moleculetypefunction(AtomTypeDummy);
        if(AtomTypeDummy==1)
        {
            NumberOfAMolecules++;
            CentralAtom[MoleculeIDDummy]=AtomID;
        }
        if(AtomTypeDummy==qA+2)
        {
            NumberOfBMolecules++;
            CentralAtom[MoleculeIDDummy]=AtomID;
        }
        if(AtomTypeDummy!=1 && AtomTypeDummy!=qA+2)
        {
            j=0;
            while(Linkers[MoleculeIDDummy][j]!=-1)
                j++;
            Linkers[MoleculeIDDummy][j]=AtomID;
        }
        if(MoleculeIDDummy>NumberOfMolecules)
            NumberOfMolecules=MoleculeIDDummy;

        FreeLinkersFlag[i+1]=0; // Initialization (0 means free, 1 means bound)
    }

    FreeLinkersFlag[0]=-1; // Initialization

    fclose(FileMovie);


	// *************************************************************
	// ************ Building adjacency matrix (bonds) **************

	printf("Building adjacency matrices...");

    for(i=1; i<=NumberOfMolecules; i++)
    {
        //printf("\n%d: ", i);
        AdjacencyMatrixBond[i][i]=0;
        for(j=i+1; j<=NumberOfMolecules; j++)
        {
            Sum=0;
            for(ki=0; Linkers[i][ki]!=-1; ki++)
            {
                for(kj=0; Linkers[j][kj]!=-1; kj++)
                {
                    k=adjacencybondfunction(Linkers[i][ki], Linkers[j][kj]);
                    if(k==1)
                    {
                        Sum+=k;
                        FreeLinkersFlag[Linkers[i][ki]]+=1;
                        FreeLinkersFlag[Linkers[j][kj]]+=1;
                        BondCounter[interactionflagfunction(AtomType[Linkers[i][ki]], AtomType[Linkers[j][kj]])]+=1;
                        // ERROR if interactionflagfunction returns -1
                    }
                }
            }
            AdjacencyMatrixBond[i][j]=Sum;
            AdjacencyMatrixBond[j][i]=Sum;
            //if(Sum!=0)
                //printf("%d %d, ", i,j);
        }
    }

//    for(i=1;i<=NumberOfMolecules;i++)
//    {
//        printf("\n");
//        for(j=1;j<=NumberOfMolecules;j++)
//            printf("%d ", AdjacencyMatrix[i][j]);
//    }


	// *****************************************************************************************
	// ************ Building adjacency matrix (based on distance between central atoms) ********

	//printf("Building adjacency matrix (Dist)...");

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
        iDone[i]=-1; // -1 = not done nor in the to-do list, 0 = in the to-do list, 1 = done
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
            }

            if(AdjacencyMatrixBond[i][j]>0)
            {
                NumberOfEdgesForThisMolecule++;
                if(MoleculeType[i]==1)
                    AConnectivity+=AdjacencyMatrixBond[i][j];
                if(MoleculeType[i]==2)
                    BConnectivity+=AdjacencyMatrixBond[i][j];
            }
        }

        EdgesInCluster[ClustersCounter]+=NumberOfEdgesForThisMolecule;  // later divided by 2, as I'm counting twice every edge
        if(NumberOfEdgesForThisMolecule==1)
            NumberOfTerminalNodes[ClustersCounter]++;
        if(NumberOfEdgesForThisMolecule==0)
            NumberOfUnboundNodes[ClustersCounter]++;

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
    BConnectivity/=NumberOfBMolecules;

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
        ClusterSize[k][0]=0;
        ClusterSize[k][1]=0;
        ClusterSize[k][2]=0;
        DistrClusterSize[k][0]=0;
        DistrClusterSize[k][1]=0;
        DistrClusterSize[k][2]=0;
    }
    NumberOfAClusters=0;
    NumberOfBClusters=0;

    for(i=1; i<=NumberOfMolecules; i++)
    {
        ClusterSize[ClusterID[i]][0]++;
        if(MoleculeType[i]==1)
        {
            ClusterSize[ClusterID[i]][1]++;
            if(ClusterSize[ClusterID[i]][1]==1)
                NumberOfAClusters++;
        }
        if(MoleculeType[i]==2)
        {
            ClusterSize[ClusterID[i]][2]++;
            if(ClusterSize[ClusterID[i]][2]==1)
                NumberOfBClusters++;
        }

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
    AvgClusterSizeA=0;
    AvgClusterSizeB=0;
    AvgClusterSizeAllBigger=0;
    AvgClusterSizeABigger=0;
    AvgClusterSizeBBigger=0;
    ClustersCounterBigger=0;

    for(k=0; k<ClustersCounter; k++)
    {
        AvgClusterSizeAll+=ClusterSize[k][0];
        if(ClusterSize[k][1]>=2 && ClusterSize[k][2]>=2)
        {
            AvgClusterSizeAllBigger+=ClusterSize[k][0];
            AvgClusterSizeABigger+=ClusterSize[k][1];
            AvgClusterSizeBBigger+=ClusterSize[k][2];
            ClustersCounterBigger+=1;
        }
        DistrClusterSize[ClusterSize[k][0]][0]+=ClusterSize[k][0]*1.0/NumberOfMolecules;
        if(ClusterSize[k][1]>0)
        {
            AvgClusterSizeA+=ClusterSize[k][1];
            DistrClusterSize[ClusterSize[k][1]][1]+=ClusterSize[k][1]*1.0/NumberOfAMolecules;
        }
        if(ClusterSize[k][2]>0)
        {
            AvgClusterSizeB+=ClusterSize[k][2];
            DistrClusterSize[ClusterSize[k][2]][2]+=ClusterSize[k][2]*1.0/NumberOfBMolecules;
        }

        // Inertia
        for(j=0; j<2; j++)
            CenterOfMass[k][j]/=ClusterSize[k][0];
        Inertia[k][0] -= ClusterSize[k][0]*SQR(CenterOfMass[k][1]); // I_xx
        Inertia[k][1] -= ClusterSize[k][0]*SQR(CenterOfMass[k][0]); // I_yy
        Inertia[k][2] += ClusterSize[k][0]*CenterOfMass[k][0]*CenterOfMass[k][1]; // I_xy
        Inertia[k][3] = Inertia[k][0]+Inertia[k][1];    // I_zz
        //printf("%d %g %g %g %g %g %g %g\n", ClusterSize[k][0], CenterOfMass[k][0], CenterOfMass[k][1], Inertia[k][0], Inertia[k][1], Inertia[k][2], Inertia[k][3], sqrt(3.0)/4/PI*SQR(SigmaLJ)*SQR(ClusterSize[k][0])/Inertia[k][3] );
        trace = Inertia[k][0]+Inertia[k][1];
        determinant= Inertia[k][0]*Inertia[k][1] - SQR(Inertia[k][2]);
        Inertia[k][4] = 0.5*(trace-sqrt(SQR(trace)-4*determinant));
        Inertia[k][5] = 0.5*(trace+sqrt(SQR(trace)-4*determinant));

    }
    AvgClusterSizeAll/=ClustersCounter;
    AvgClusterSizeA/=NumberOfAClusters;
    AvgClusterSizeB/=NumberOfBClusters;

    if(ClustersCounterBigger!=0)
    {
        AvgClusterSizeAllBigger/=ClustersCounterBigger;
        AvgClusterSizeABigger/=ClustersCounterBigger;
        AvgClusterSizeBBigger/=ClustersCounterBigger;
    }
    sortclustersfunction();


    // ************************************************************
    // **************** Compute fraction of small clusters  *******

    // Same thing as the current way of calculating the Cluster Size distribution (which is actually molecules distribution according to size of the cluster they are in)

    FractionOfAInClustersOfSizeUpTo1=DistrClusterSize[1][1];
    FractionOfAInClustersOfSizeUpTo5=0;
    for(k=1; k<=5; k++)
        FractionOfAInClustersOfSizeUpTo5+=DistrClusterSize[k][1];
    FractionOfAInClustersOfSizeUpTo10=0;
    for(k=1; k<=10; k++)
        FractionOfAInClustersOfSizeUpTo10+=DistrClusterSize[k][1];

    FractionOfFreeB=0;
    for(i=1; i<=NumberOfMolecules; i++)
    {
        if(MoleculeType[i]==2)
        {
            if(ClusterSize[ClusterID[i]][0]==1)
                FractionOfFreeB+=1;
        }
    }
    FractionOfFreeB/=NumberOfBMolecules;


    // ******************************************************
    // **************** Compute CoalescenceLikelihood *******

    for(i=1; i<=NumberOfMolecules; i++)
    {
        for(k=0; k<5; k++)
        {
            linker=Linkers[i][k];
            if(linker!=-1)
            {
                if(FreeLinkersFlag[linker]==0)
                    FreeLinkersNumberInCluster[ClusterID[i]][AtomType[linker]]++;
            }
        }
    }

    CoalescenceLikelihood1=0;
    CoalescenceLikelihood2=0;

    for(i=0; i<ClustersCounter; i++)
    {
        for(j=i+1; j<ClustersCounter; j++)
        {
            CoalescenceLikelihood1+=coalescence1function(i,j,1);
            CoalescenceLikelihood2+=coalescence2function(i,j);
        }
    }


  	// **********************************************************************
    // ************** Print output ******************************************

    printf("Printing output file...\n");
    /*FirstCharToCopy=16;
    for(i=FirstCharToCopy; i<strlen(FileNameMovie)-4; i++)
        string[i-FirstCharToCopy]=FileNameMovie[i];
    string[i-FirstCharToCopy]='\0';
    sprintf(FileNameOutput,"ClusterAnalysis_");
    strcat(FileNameOutput,string);
    strcat(FileNameOutput,"_ts");
    strcat(FileNameOutput,argv[2]);
    strcat(FileNameOutput,".dat");
    */
    printf("%s\n", FileNameOutput);


    FileOutput=fopen(FileNameOutput,"wb+");


    fprintf(FileOutput, "InputFile\n%s\nData analysis program\nClusterAnalysisFluidljcossq1.2.exe\nts\n%d\n\n", FileNameMovie,timestep);
    fprintf(FileOutput, "NumberOfClusters (number of clusters, irrespective of the type of molecule they are composed of)\n%d\n", ClustersCounter);
    fprintf(FileOutput, "NumberOfClustersBigger (number of clusters containing at least 2 A and 2 B molecules)\n%d\n", ClustersCounterBigger);
    fprintf(FileOutput, "NumberOfAClusters (number of clusters containing at least one A molecule)\n%d\n", NumberOfAClusters);
    fprintf(FileOutput, "\nAvgClusterSize (total number of molecules / NumberOfClusters)\n%f\n", AvgClusterSizeAll);
    fprintf(FileOutput, "AvgAClusterSize (number of A molecules / NumberOfAClusters)\n%f\n", AvgClusterSizeA);
    fprintf(FileOutput, "AvgBClusterSize (number of B molecules / NumberOfBClusters)\n%f\n", AvgClusterSizeB);
    fprintf(FileOutput, "\nAvgClusterSizeBigger (total number of molecules in clusters with at least 2 A and 2 B / number of such clusters)\n%f\n", AvgClusterSizeAllBigger);
    fprintf(FileOutput, "AvgAClusterSizeBigger (number of A molecules in clusters with at least 2 A and 2 B / number of such clusters)\n%f\n", AvgClusterSizeABigger);
    fprintf(FileOutput, "AvgBClusterSizeBigger (number of B molecules in clusters with at least 2 A and 2 B / number of such clusters)\n%f\n", AvgClusterSizeBBigger);
    fprintf(FileOutput, "\nAvgAConnectivity (average number of molecules bound to one A molecule)\n%f\n", AConnectivity);
    if(NumberOfBMolecules!=0)
        fprintf(FileOutput, "AvgBConnectivity (average number of molecules bound to one B molecule)\n%f\n", BConnectivity);
    else
        fprintf(FileOutput, "AvgBConnectivity (average number of molecules bound to one B molecule)\n0\n");
    fprintf(FileOutput, "\nCoalescenceLikelihood1 (number of combinations between: 1- free linkers of each cluster and 2- free linkers of each other cluster that the former linkers can bind to)\n%lld\n", CoalescenceLikelihood1);
    fprintf(FileOutput, "CoalescenceLikelihood2 (as before, but clusters are supposed circular and compact, so that only a number of free linkers proportional to sqrt(sizeofthecluster) is considered)\n%f\n", CoalescenceLikelihood2);
    fprintf(FileOutput, "\nFractionOfAInClustersOfSizeUpTo1 (fraction of A molecules in clusters with only 1 A molecule)\n%f\n", FractionOfAInClustersOfSizeUpTo1);
    fprintf(FileOutput, "FractionOfAInClustersOfSizeUpTo5 (fraction of A molecules in clusters with no more than 5 A molecules)\n%f\n", FractionOfAInClustersOfSizeUpTo5);
    fprintf(FileOutput, "FractionOfAInClustersOfSizeUpTo10 (fraction of A molecules in clusters with no more than 10 A molecules)\n%f\n", FractionOfAInClustersOfSizeUpTo10);
    fprintf(FileOutput, "\nFractionOfFreeB (fraction of B molecules in a cluster of unit size)\n%f\n", FractionOfFreeB);

    fprintf(FileOutput, "\nBondCounter (number of active bonds between binding sites of a given kind), CoalescenceBondCounter (number of possible bonds between free binding sites of a given kind, belonging to two different clusters)\n");
    for(k=0; k<NumberOfInteractions; k++)
        fprintf(FileOutput, "%d-%d %lld %lld\n", InteractionsMatrix[k][0], InteractionsMatrix[k][1], BondCounter[k], CoalescenceBondCounter[k]);

    fprintf(FileOutput, "\nClusters detail\n");
    fprintf(FileOutput, "ClusterSize[][0All], ClusterSize[][1A], ClusterSize[][2B], EdgesInCluster, NumberOfTerminalNodes, CenterOfMass_x, CenterOfMass_y, Inertia_lambda1, Inertia_lambda2, Inertia_zz, Compactness, Roundness\n"); // Roundeness = Sphericity
    fprintf(FileOutput, " > NumberOfUnboundNodes\n"); //gtinertia, gtinertia*/ClusterSize[][0]^2, gtdiameter, gtdiameter/gtMooreBound\n");
    fprintf(FileOutput, " >> ");
    for (i=0; i<=qA; i++)
        fprintf(FileOutput, "Adeg%d, ", i);
    for (i=0; i<qB; i++)
        fprintf(FileOutput, "Bdeg%d, ", i);
    fprintf(FileOutput, "Bdeg%d\n", qB);
    for(k=0; k<ClustersCounter; k++)
    {
        i=OrderedClustersIndex[k];
        Sphericity=Inertia[i][4]/Inertia[i][5];
        fprintf(FileOutput, "%d %d %d %d %d %g %g %g %g %g %g %g\n", ClusterSize[i][0], ClusterSize[i][1], ClusterSize[i][2], EdgesInCluster[i], NumberOfTerminalNodes[i], CenterOfMass[i][0], CenterOfMass[i][1], Inertia[i][4], Inertia[i][5], Inertia[i][3], /*sqrt(3.0)/4/PI*SQR(SigmaLJ)*SQR(ClusterSize[i][0])/Inertia[i][3],*/ sqrt(3.0)/8/PI*SQR(SigmaLJ)*SQR(ClusterSize[i][0])*(1+Sphericity)/sqrt(Sphericity)/Inertia[i][3], Sphericity );
        graphtheoreticalinertiafunction(i); // Only to compute gtNodesByDegree
        if(NumberOfUnboundNodes[i]!=gtNodesByDegree[1][0]+gtNodesByDegree[2][0])
        {
            printf("*** ERROR: NumberOfUnboundNodes!=Adeg0+Bdeg0 for cluster %d\n", k);
            return 501;
        }
        fprintf(FileOutput, " > %d\n", NumberOfUnboundNodes[i]); //%lli %f %d %f\n", gtinertia, gtinertia*1.0/SQR(ClusterSize[i][0]), gtdiameter, gtdiameter/gtMooreBound);
        fprintf(FileOutput, " >> ");
        for (i=0; i<=qA; i++)
            fprintf(FileOutput, "%d ", gtNodesByDegree[1][i]);
        for (i=0; i<=qB; i++)
            fprintf(FileOutput, "%d ", gtNodesByDegree[2][i]);
        fprintf(FileOutput, "\n");
    }

    fprintf(FileOutput, "\nDistrClusterSizeAll (fraction of particles (A+B, A, B) in clusters, by size of the clusters)\n");
    for(k=0; k<=NumberOfMolecules; k++)
        fprintf(FileOutput, "%d %f %f %f\n", k, DistrClusterSize[k][0], DistrClusterSize[k][1], DistrClusterSize[k][2]);

    fclose(FileOutput);


	return 0;
}












