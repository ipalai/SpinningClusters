// // 2.0 Adapted from TheoryClustersFluidljcossq's ClusterInternalDiffusivity1.2.c for SpinningClusters

// From ClusterAnalysisFluidljcossq1.2.c

// USAGE: ./$executable $xyzInputFile $OutputFile $Timestep1 $Timestep2 $Timestep3 ...

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#define NMAX 12000
#define NMOLECULESMAX 2000
#define NATOMTYPESMAX 14
#define MOLECULETYPES 2
#define VALENCEMAX 6
#define NINTERACTIONSMAX 36
#define NTIMESTEPSMAX 15
#define NPAIRSMAX 101
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
int CentralAtom[NMOLECULESMAX];
float CutoffDistanceDiffu, CutoffDistanceDiffuSq;


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


float pairdist2function(int i, int j)
{
    float dx, dy;
    dx = fabsf(Position[i][0]-Position[j][0]);
    dx = fminfunction(dx, BoxLength-dx);
    dy = fabsf(Position[i][1]-Position[j][1]);
    dy = fminfunction(dy, BoxLength-dy);
    return (SQR(dx)+SQR(dy));
}


int numberdiffneighboursfunction(int m0)  // Computes #neighbouring molecules of molecule m0, with cutoff distance CutoffDistanceDiffu (in general different from CutoffDistanceDist)
{
    int Sum, m;
    Sum=0;
    for(m=0;m<NumberOfMolecules;m++)
    {
        if(AdjacencyMatrixDist[m0][m]>0)
        {
            if(pairdist2function(CentralAtom[m0],CentralAtom[m])<CutoffDistanceDiffuSq)
                Sum++;
        }
    }
    return Sum;
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
	int Timesteps[NTIMESTEPSMAX];
	int MoleculesInBiggestCluster[NMOLECULESMAX];
	int NumberOfAClusters, NumberOfBClusters, NumberOfAMolecules, NumberOfBMolecules, CheckCounter;
	int NumberOfAtoms, AtomID, MoleculeIDDummy, MoleculeID[NMAX], AtomTypeDummy;
	int Linkers[NMOLECULESMAX][VALENCEMAX], linker;
	int iToDo[NMOLECULESMAX], iToDoCounter, iToDoSegno, Sum, iDone[NMOLECULESMAX], FreeLinkersFlag[NMAX];
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
	char string[strlength], FileNameOutput[strlength], FileNameMovie[strlength], timestepString[strlength], timestepStrings[NTIMESTEPSMAX][strlength];
	char *p;
	char *pch;
	FILE *FileOutput;
	FILE *FileMovie;

    int NTimesteps, TimestepCounter;
	int PairsToFollow[NPAIRSMAX][2];
	int MoleculesToFollow[2*NPAIRSMAX];
	int MoleculesToFollowCounter;
	float r1, r2;
    int dummyint;
    int BiggestClusterSize0, NPairsToFollow;
    int NewPairFlag;
    float PairMSD[NTIMESTEPSMAX][NPAIRSMAX];
    float PairDistance[NTIMESTEPSMAX][NPAIRSMAX]; //PairInitialDistanceVector[NPAIRSMAX][2], PairInitialDistance2[NPAIRSMAX];
    int PairInClusterFlag[NPAIRSMAX];
    int FollowFlag0, FollowFlag1;
    int MoleculesInThisCluster;
    int m0, m1, m;
    float dummyDisplacementx, dummyDisplacementy, dummyNormProduct;
    float CosTheta, SinTheta;
    float PositionsCMRFTime0[NMOLECULESMAX][2]; // Positions in the center of mass ref frame at first timestep read
    float RelativePositionx, RelativePositiony, DeltaDeltaDiffusionx, DeltaDeltaDiffusiony;
    int NumberOfNeighboursOfm0, NumberOfNeighboursOfm1;
    int ThresholdNeighbours;
    int PairAttemptsCounter=0;
    //clock_t clockt;

	printf("\nInitialization...\n");
    srand(1234);

	// Get arguments
    strcpy(FileNameMovie,argv[1]);
    strcpy(FileNameOutput, argv[2]);
    NTimesteps=argc-3;
    if(NTimesteps>NTIMESTEPSMAX)
    {
        printf("*** ERROR: Too many timesteps (max=%d)", NTIMESTEPSMAX);
        return 300;
    }
    printf("Timesteps: ");
    for(i=0;i<NTimesteps;i++)
    {
        sprintf(timestepStrings[i], "%s\n", argv[i+3]);
        sscanf(timestepStrings[i], "%d", &Timesteps[i]);
        Timesteps[i] = strtol(argv[i+3], &p, 10);
        printf("%d ", Timesteps[i]);
    }

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
    printf("\nqA %d, qB %d, ra %f, CutoffDistanceDist %f\n", qA, qB, ra, CutoffDistanceDist);

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


    NPairsToFollow=100;
    ThresholdNeighbours=3;
    CutoffDistanceDiffu=1.3*SigmaLJ;
    CutoffDistanceDiffuSq=SQR(CutoffDistanceDiffu);
    printf("NPairsToFollow %d, ThresholdNeighbours %d, CutoffDistanceDiffu % f\n",NPairsToFollow, ThresholdNeighbours, CutoffDistanceDiffu);
    if(NPairsToFollow>NPAIRSMAX)
    {
        printf("*** ERROR: NPairsToFollow>NPAIRSMAX (%d>%d)", NPairsToFollow, NPAIRSMAX);
        return 305;
    }



    // Iterate over Timesteps

    FileMovie=fopen(FileNameMovie,"r");   // Read file
    printf("\nReading file %s...\n", FileNameMovie);

    for(TimestepCounter=0;TimestepCounter<NTimesteps;TimestepCounter++)
    {
        timestep=Timesteps[TimestepCounter];
        strcpy(timestepString, timestepStrings[TimestepCounter]);
        printf("\ntimestep %d=%s", timestep, timestepString);

        // Initialization of clustering algorithm
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



        // Find desired TIMESTEP in Movie file
        FlagRead=0;
        while(FlagRead==0)
        {
            fgets(string, strlength, FileMovie);
            //printf("%s",string);
            if(strcmp(string, "ITEM: TIMESTEP\n")==0)
            {
                fgets(string, strlength, FileMovie);
                //printf("%s",string);
                if(strcmp(string,timestepStrings[TimestepCounter])==0)
                    FlagRead=1;
            }
            if(feof(FileMovie)!=0)
            {
                printf("*** ERROR: FileMovie does not contain TIMESTEP %s", timestepString);
                return 1;
            }
        }
        //printf("Checkpoint7");
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


        // *************************************************************
        // ************ Building adjacency matrix (bonds) **************
        //clockt=clock();

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

        //clockt = clock()-clockt;

        // ************************************************************
        // **************** Compute fraction of small clusters  *******

        // Same thing as the current way of calculating the Cluster Size distribution (which is actually molecules distribution according to size of the cluster they are in)

        /*
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
        */

        // ******************************************************
        // **************** Compute CoalescenceLikelihood *******

        /*
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
        */


        // ************************************************************************************************
        // **************** Compute Internal MSD for a bunch of particle pairs ******

        // If this is the first analysed timestep, choose which particle pairs to track
        if(TimestepCounter==0)
        {
            printf("Choosing molecules to follow (only at initial timestep)\n");
            k=OrderedClustersIndex[0];
            BiggestClusterSize0=ClusterSize[k][0];
            printf("BiggestClusterSize0 %d\n", BiggestClusterSize0);
            // Check that NPairsToFollow is not too large for the size of the cluster. If so, set it to half of its max theoretical value (which is BiggestClusterSize0 choose 2)
            if(NPairsToFollow > BiggestClusterSize0) //*(BiggestClusterSize0-1)/2/2)
            {
                NPairsToFollow = BiggestClusterSize0; //*(BiggestClusterSize0-1)/2/2;
            }
            // Build list of particles in biggest cluster (only at 1st timestep) and compute positions in the CM reference frame
            j=0;
            for(i=0; i<NumberOfMolecules; i++)
            {
                if(ClusterID[i]==k)
                {
                    MoleculesInBiggestCluster[j]=i;
                    j++;
                }
                PositionsCMRFTime0[i][0]= Position[CentralAtom[i]][0] + BoxLength*WindingNumber[i][0] - CenterOfMass[ClusterID[i]][0];
                PositionsCMRFTime0[i][1]= Position[CentralAtom[i]][1] + BoxLength*WindingNumber[i][1] - CenterOfMass[ClusterID[i]][1];
            }
            j=0;
            //printf("Built list of molecules in biggest cluster. NPairsToFollow %d\n",NPairsToFollow);


            // Select NPairsToFollow pairs at random from the biggest cluster: the two molecules in the pair have to be nearest neighbors at the initial timestep
            MoleculesToFollowCounter=0;
            PairAttemptsCounter=0;
            for (i=0;i<NPairsToFollow;i++)
            {
                NewPairFlag=1;
                PairInClusterFlag[i]=1;
                while(NewPairFlag==1)
                {
                    // Extract the pair at random
                    r1 = rand();      // Returns a pseudo-random integer between 0 and RAND_MAX.
                    PairsToFollow[i][0]=MoleculesInBiggestCluster[(int) (r1/RAND_MAX * BiggestClusterSize0)];
                    r2 = rand();
                    while( ((int)(r2/RAND_MAX * BiggestClusterSize0))==((int)(r1/RAND_MAX * BiggestClusterSize0)) )
                        r2=rand();
                    PairsToFollow[i][1]=MoleculesInBiggestCluster[(int) (r2/RAND_MAX * BiggestClusterSize0)];
                    // Order the pair
                    if(PairsToFollow[i][0]>PairsToFollow[i][1])
                    {
                        dummyint=PairsToFollow[i][0];
                        PairsToFollow[i][0]=PairsToFollow[i][1];
                        PairsToFollow[i][1]=dummyint;
                    }

                    // Check that the two chosen molecules are neighbors
                    NewPairFlag=0;
                    m0=PairsToFollow[i][0];
                    m1=PairsToFollow[i][1];
                    RelativePositionx = (Position[CentralAtom[m1]][0] + BoxLength*WindingNumber[m1][0]) - (Position[CentralAtom[m0]][0] + BoxLength*WindingNumber[m0][0]) ;  // x1(t)-x0(t) component x
                    RelativePositiony = (Position[CentralAtom[m1]][1] + BoxLength*WindingNumber[m1][1]) - (Position[CentralAtom[m0]][1] + BoxLength*WindingNumber[m0][1]) ;  // component y
                    if(SQR(RelativePositionx)+SQR(RelativePositiony)>CutoffDistanceDiffuSq )
                    {
                        NewPairFlag=1;
                        continue;
                    }

                    // Check that the extracted pair was not already chosen
                    for(j=0;j<i;j++)
                    {
                        if((PairsToFollow[j][0]==PairsToFollow[i][0]) && (PairsToFollow[j][0]==PairsToFollow[i][0]))
                        {
                            NewPairFlag=1;
                            break;
                        }
                    }

                    // Check that both chosen molecules have at least 3 neighbours
                    if(numberdiffneighboursfunction(PairsToFollow[i][0])<ThresholdNeighbours || numberdiffneighboursfunction(PairsToFollow[i][1])<ThresholdNeighbours)
                        NewPairFlag=1;

                    // If too many attempts, end program
                    PairAttemptsCounter++;
                    if(PairAttemptsCounter>NPairsToFollow*1000)
                    {
                        printf("\n----- Impossible to build initial PairsToFollow list (stuck at %d th pair). Probably cluster is 1D.\n", i);
                        FileOutput=fopen(FileNameOutput,"wb+");
                        fprintf(FileOutput, "InputFile\n%s\nData analysis program\nClusterInternalDiffusivityFluidljcossq1.0.exe\n", FileNameMovie);
                        fprintf(FileOutput, "\nNPairsToFollow\n%d\nThresholdNeighbours\n%d\nCutoffDistanceDiffu\n%f\n",NPairsToFollow, ThresholdNeighbours,CutoffDistanceDiffu);
                        fprintf(FileOutput, "\nTimesteps\n");
                        for(j=0;j<NTimesteps;j++)
                            fprintf(FileOutput, "%d ", Timesteps[j]);
                        fprintf(FileOutput, "\n\nImpossible to build initial PairsToFollow list (stuck at %d th pair). Probably cluster is 1D.\n", i);
                        return 133;
                    }
                }
                //printf("Chosen pair %d %d\n", PairsToFollow[i][0],PairsToFollow[i][1]);
                // Store initial distance vector
                m0=PairsToFollow[i][0];
                m1=PairsToFollow[i][1];
                //PairInitialDistanceVector[i][0]= Position[CentralAtom[m1]][0] + BoxLength*WindingNumber[m1][0] - ( Position[CentralAtom[m0]][0] + BoxLength*WindingNumber[m0][0] );
                //PairInitialDistanceVector[i][1]= Position[CentralAtom[m1]][1] + BoxLength*WindingNumber[m1][1] - ( Position[CentralAtom[m0]][1] + BoxLength*WindingNumber[m0][1] );
                //PairInitialDistance2[i]= SQR(PairInitialDistanceVector[i][0]) + SQR(PairInitialDistanceVector[i][1]);


                // Create list of molecules to follow
                if(MoleculesToFollowCounter==0)
                {
                    MoleculesToFollow[0]=PairsToFollow[i][0];
                    MoleculesToFollow[1]=PairsToFollow[i][1];
                    MoleculesToFollowCounter=2;
                    //printf("MoleculesToFollow %d %d\n", PairsToFollow[i][0], PairsToFollow[i][1]);
                }
                else
                {
                    FollowFlag0=0;
                    FollowFlag1=0;
                    for(j=0;j<=MoleculesToFollowCounter;j++)
                    {
                        if(MoleculesToFollow[j]==PairsToFollow[i][0])
                            FollowFlag0=1;
                        if(MoleculesToFollow[j]==PairsToFollow[i][1])
                            FollowFlag1=1;
                    }
                    if(FollowFlag0==0)
                    {
                        MoleculesToFollow[MoleculesToFollowCounter]=PairsToFollow[i][0];
                        MoleculesToFollowCounter++;
                        //printf("MoleculesToFollow %d\n", MoleculesToFollow[MoleculesToFollowCounter-1]);
                    }
                    if(FollowFlag1==0)
                    {
                        MoleculesToFollow[MoleculesToFollowCounter]=PairsToFollow[i][1];
                        MoleculesToFollowCounter++;
                        //printf("MoleculesToFollow %d\n", MoleculesToFollow[MoleculesToFollowCounter-1]);
                    }
                }
            }
        }

        // At any timestep, compute rotation and diffusive MSD
        printf("Computing rotation and MSD...");
        // Identify what cluster k contains most of particles, meanwhile compute rotation
        for(j=0; j<ClustersCounter; j++)
        {
            CosTheta=0;
            SinTheta=0;
            k=OrderedClustersIndex[j];
            MoleculesInThisCluster=0;
            for(i=0;i<MoleculesToFollowCounter;i++)
            {
                m0=MoleculesToFollow[i];
                if(ClusterID[m0]==k)
                {
                    MoleculesInThisCluster++;
                    RelativePositionx = (Position[CentralAtom[m0]][0] + BoxLength*WindingNumber[m0][0]) - CenterOfMass[k][0] ; // Relative to CM
                    RelativePositiony = (Position[CentralAtom[m0]][1] + BoxLength*WindingNumber[m0][1]) - CenterOfMass[k][1] ;
                    dummyNormProduct = sqrt( SQR(RelativePositionx) + SQR(RelativePositiony) ) * sqrt( SQR(PositionsCMRFTime0[m0][0]) + SQR(PositionsCMRFTime0[m0][1]) ) ;
                    CosTheta+= (RelativePositionx*PositionsCMRFTime0[m0][0] + RelativePositiony*PositionsCMRFTime0[m0][1] ) / dummyNormProduct ; // from dot product x'(t).x'(0)
                    SinTheta+= (PositionsCMRFTime0[m0][0]*RelativePositiony - PositionsCMRFTime0[m0][1]*RelativePositionx ) / dummyNormProduct ; // from cross product x'(0) cross x'(t)
                }
            }
            CosTheta/=MoleculesInThisCluster;
            SinTheta/=MoleculesInThisCluster;
            if(MoleculesInThisCluster>=MoleculesToFollowCounter/2)
                break;
        }
        printf(" Cluster index k=%d, current MoleculesInThisCluster=%d, initial MoleculesToFollowCounter=%d, CosTheta=%f, SinTheta=%f\n", k, MoleculesInThisCluster,MoleculesToFollowCounter,CosTheta,SinTheta);

        // If cluster still exists (i.e. >50% molecules are still in it), compute distances

        if(MoleculesInThisCluster>MoleculesToFollowCounter/2)
        {
            for(i=0; i<NPairsToFollow; i++)
            {
                // Check that a pair stayed in the same cluster and that it has always belonged to that cluster. If so compute MSD and Distance, if not set them to -1 and remove pair from PairInClusterFlag.
                m0=PairsToFollow[i][0];
                m1=PairsToFollow[i][1];
                if((ClusterID[PairsToFollow[i][0]]==k) && (ClusterID[PairsToFollow[i][1]]==k) && (numberdiffneighboursfunction(m0)>=ThresholdNeighbours) && (numberdiffneighboursfunction(m1)>=ThresholdNeighbours) && PairInClusterFlag[i]==1 )
                {
                    // Invert rotation and compute DeltaDeltaDiffusion, i.e. the relative displacement vector due to diffusion, using as rotation state the one of t=0 (i.e. no rotation)
                    RelativePositionx = (Position[CentralAtom[m1]][0] + BoxLength*WindingNumber[m1][0]) - (Position[CentralAtom[m0]][0] + BoxLength*WindingNumber[m0][0]) ;  // x1(t)-x0(t) component x
                    RelativePositiony = (Position[CentralAtom[m1]][1] + BoxLength*WindingNumber[m1][1]) - (Position[CentralAtom[m0]][1] + BoxLength*WindingNumber[m0][1]) ;  // component y
                    DeltaDeltaDiffusionx = CosTheta*RelativePositionx + SinTheta*RelativePositiony; // rotate relative position at time t back in the rotation state of t=0
                    DeltaDeltaDiffusiony = -SinTheta*RelativePositionx + CosTheta*RelativePositiony;
                    DeltaDeltaDiffusionx -= PositionsCMRFTime0[m1][0] - PositionsCMRFTime0[m0][0] ; // Subtract initial relative distance (PositionsCMRFTime0 is x'(0), i.e. in the CM ref frame, but it doesn't matter)
                    DeltaDeltaDiffusiony -= PositionsCMRFTime0[m1][1] - PositionsCMRFTime0[m0][1] ; //
                    PairMSD[TimestepCounter][i] = SQR(DeltaDeltaDiffusionx) + SQR(DeltaDeltaDiffusiony); // this is the diffusive MSD
                    //    pairdist2function(CentralAtom[PairsToFollow[i][0]],CentralAtom[PairsToFollow[i][j]]);
                    //printf("MSD %d %d %f\n", m0, m1, PairMSD[TimestepCounter][i]);
                    PairDistance[TimestepCounter][i] = sqrt(SQR(RelativePositionx)+SQR(RelativePositiony));
                }
                else
                {
                    PairMSD[TimestepCounter][i]=-1;
                    PairDistance[TimestepCounter][i]=-1;
                    PairInClusterFlag[i]=0;
                }
            }
        }
        // If cluster has broken, set PairMSD and PairDistance to -1
        else
        {
            printf("Previous cluster does not exist anymore at timestep %d (half of its molecules have escaped).\n", timestep);
            for(i=0;i<NPairsToFollow;i++)
            {
                PairMSD[TimestepCounter][i]=-1;
                PairDistance[TimestepCounter][i]=-1;
            }
        }

        printf("End of MSD calculation.\n");
        //printf("fun() took %f seconds to execute \n", ((double)clockt)/CLOCKS_PER_SEC);



    // ***************************
    // End of cycle over timesteps
    }


    fclose(FileMovie);


    // **********************************************************************
    // ************** Print output ******************************************

    printf("\nPrinting output file...\n");
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


    fprintf(FileOutput, "InputFile\n%s\nData analysis program\nClusterInternalDiffusivityFluidljcossq1.0.exe\n", FileNameMovie);
    fprintf(FileOutput, "\nNPairsToFollow\n%d\nThresholdNeighbours\n%d\nCutoffDistanceDiffu\n%f\n",NPairsToFollow, ThresholdNeighbours,CutoffDistanceDiffu);
    fprintf(FileOutput, "\nTimesteps\n");
    for(j=0;j<NTimesteps;j++)
        fprintf(FileOutput, "%d ", Timesteps[j]);

    fprintf(FileOutput,"\n\nMSDAtTimesteps\nMolecule_i Molecule_j Timesteps");
    for(i=0;i<NPairsToFollow;i++)
    {
        fprintf(FileOutput, "\n%d %d ", PairsToFollow[i][0], PairsToFollow[i][1]);
        for(j=0;j<NTimesteps;j++)
            fprintf(FileOutput, "%f ", PairMSD[j][i]);
    }

    fprintf(FileOutput,"\n\nDistanceAtTimesteps\nMolecule_i Molecule_j Timesteps");
    for(i=0;i<NPairsToFollow;i++)
    {
        fprintf(FileOutput, "\n%d %d ", PairsToFollow[i][0], PairsToFollow[i][1]);
        for(j=0;j<NTimesteps;j++)
            fprintf(FileOutput, "%f ", PairDistance[j][i]);
    }

    /*
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
    fprintf(FileOutput, " > NumberOfUnboundNodes\n"); //gtinertia, gtinertia_/ClusterSize[][0]^2, gtdiameter, gtdiameter/gtMooreBound\n");
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
        fprintf(FileOutput, "%d %d %d %d %d %g %g %g %g %g %g %g\n", ClusterSize[i][0], ClusterSize[i][1], ClusterSize[i][2], EdgesInCluster[i], NumberOfTerminalNodes[i], CenterOfMass[i][0], CenterOfMass[i][1], Inertia[i][4], Inertia[i][5], Inertia[i][3], //sqrt(3.0)/4/PI*SQR(SigmaLJ)*SQR(ClusterSize[i][0])/Inertia[i][3], // sqrt(3.0)/8/PI*SQR(SigmaLJ)*SQR(ClusterSize[i][0])*(1+Sphericity)/sqrt(Sphericity)/Inertia[i][3], Sphericity );
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

    */
    fclose(FileOutput);


	return 0;
}












