// 2.0 Adapted from TheoryClustersFluidljcossq's ClusterAnalysisStatistics1.1.c as a clustering code for SpinningClusters

// ClusterAnalysisStatisticsFluidljcossq1.1.c Branch from ClusterAnalysisStatistics1.0.c, to account for fluid clusters.
// 	 Corrected a few things, then introduced in original ClusterAnalysisStatistics1.1.c as version 1.1

// 0.1 Computes standard deviations for every quantity
// 0.3 Reads CoalescenceLikelyhood's, when present
// 0.4 Reads BondCounter and CoalescenceBondCounter
// 0.5 Reads moments of inertia
// 0.6 Reads graph-theoretical data (diameter, moment of inertia, NodesByDegree)

// 1.0 adjusted for TheoryClusters project
//      updated ERROR sscanf 1a, 1b, 1c to avoid error when ClusterSize is too small and Compactness=inf

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#define NMAX 5000
#define NINTERACTIONSMAX 36
#define strlength 1500
#define SQR(x) ((x)*(x))
#define PI 3.14159265359
#define MAXFILES 100


int main(int argc, char *argv[])
{
    FILE *InputFiles[argc-2];
    FILE *OutputFile;
    int i, j, k;
    int NumberOfInputFiles;
    char string[argc-2][strlength];
    char NumberOfClustersString[strlength], NumberOfAClustersString[strlength], AvgClusterSizeString[strlength], AvgAClusterSizeString[strlength], AvgBClusterSizeString[strlength];
    char NumberOfClustersBiggerString[strlength], AvgClusterSizeBiggerString[strlength], AvgAClusterSizeBiggerString[strlength], AvgBClusterSizeBiggerString[strlength];
    char AvgAConnectivityString[strlength]="\0";
    char AvgBConnectivityString[strlength]="\0";
    char CoalescenceLikelihood1String[strlength]="\0";
    char CoalescenceLikelihood2String[strlength]="\0";
    char FractionOfAInClustersOfSizeUpTo1String[strlength]="\0";
    char FractionOfAInClustersOfSizeUpTo5String[strlength]="\0";
    char FractionOfAInClustersOfSizeUpTo10String[strlength]="\0";
    char FractionOfFreeBString[strlength]="\0";
    char BondCountersString[strlength]="\0";
    char ClustersInertiaString[strlength]="\0";
    char DistrClusterSizeString[strlength];
    char *token;
    char *p;
    int NumberOfClusters;
    float NumberOfClustersSum=0;
    float NumberOfClustersVariance=0;
    int NumberOfClustersCounter=0;
    int NumberOfClustersBigger;
    float NumberOfClustersBiggerSum=0;
    float NumberOfClustersBiggerVariance=0;
    int NumberOfClustersBiggerCounter=0;
    int NumberOfAClusters;
    float NumberOfAClustersSum=0;
    float NumberOfAClustersVariance=0;
    int NumberOfAClustersCounter=0;
    float AvgClusterSize;
    float AvgClusterSizeSum=0;
    float AvgClusterSizeVariance=0;
    int AvgClusterSizeCounter=0;
    float AvgAClusterSize;
    float AvgAClusterSizeSum=0;
    float AvgAClusterSizeVariance=0;
    int AvgAClusterSizeCounter=0;
    float AvgBClusterSize;
    float AvgBClusterSizeSum=0;
    float AvgBClusterSizeVariance=0;
    int AvgBClusterSizeCounter=0;
    float AvgClusterSizeBigger;
    float AvgClusterSizeBiggerSum=0;
    float AvgClusterSizeBiggerVariance=0;
    int AvgClusterSizeBiggerCounter=0;
    float AvgAClusterSizeBigger;
    float AvgAClusterSizeBiggerSum=0;
    float AvgAClusterSizeBiggerVariance=0;
    int AvgAClusterSizeBiggerCounter=0;
    float AvgBClusterSizeBigger;
    float AvgBClusterSizeBiggerSum=0;
    float AvgBClusterSizeBiggerVariance=0;
    int AvgBClusterSizeBiggerCounter=0;
    float AvgAConnectivity;
    float AvgAConnectivitySum=0;
    float AvgAConnectivityVariance=0;
    int AvgAConnectivityCounter=0;
    float AvgBConnectivity;
    float AvgBConnectivitySum=0;
    float AvgBConnectivityVariance=0;
    int AvgBConnectivityCounter=0;
    float CoalescenceLikelihood1;
    float CoalescenceLikelihood1Sum=0;
    float CoalescenceLikelihood1Variance=0;
    int CoalescenceLikelihood1Counter=0;
    float CoalescenceLikelihood2;
    float CoalescenceLikelihood2Sum=0;
    float CoalescenceLikelihood2Variance=0;
    int CoalescenceLikelihood2Counter=0;
    float FractionOfAInClustersOfSizeUpTo1;
    float FractionOfAInClustersOfSizeUpTo1Sum=0;
    float FractionOfAInClustersOfSizeUpTo1Variance=0;
    int FractionOfAInClustersOfSizeUpTo1Counter=0;
    float FractionOfAInClustersOfSizeUpTo5;
    float FractionOfAInClustersOfSizeUpTo5Sum=0;
    float FractionOfAInClustersOfSizeUpTo5Variance=0;
    int FractionOfAInClustersOfSizeUpTo5Counter=0;
    float FractionOfAInClustersOfSizeUpTo10;
    float FractionOfAInClustersOfSizeUpTo10Sum=0;
    float FractionOfAInClustersOfSizeUpTo10Variance=0;
    int FractionOfAInClustersOfSizeUpTo10Counter=0;
    float FractionOfFreeB;
    float FractionOfFreeBSum=0;
    float FractionOfFreeBVariance=0;
    int FractionOfFreeBCounter=0;
    float BondCounters[2];
    float BondCountersSum[NINTERACTIONSMAX][2];
    float BondCountersVariance[NINTERACTIONSMAX][2];
    int BondCountersCounter=0;
    char BondTypeString[NINTERACTIONSMAX][strlength];
    char DummyString[strlength];
    int NumberOfInteractions;

    int ClustersInertiaCounter=0;
    //double ClustersInertia[4];
    //double ClustersInertiaSum[4];   // NormalizedIzz, asphericity
    //double ClustersInertiaVariance[4];
    int ClusterSize, EdgesInCluster, NumberOfTerminalNodes, NumberOfUnboundNodes;
    int ClusterSizeThreasholdForInertia=30;
    float Compactness, Roundness;
    double EdgesInClusterSum=0, NumberOfTerminalNodesSum=0;
    double CompactnessSum=0;
    double RoundnessSum=0;
    double NumberOfUnboundNodesSum=0;
    double EdgesInClusterVariance=0;
    double NumberOfTerminalNodesVariance=0;
    double CompactnessVariance=0;
    double RoundnessVariance=0;
    double NumberOfUnboundNodesVariance=0;
    float gtDiameterNorm, gtInertiaNorm;
    double gtDiameterNormSum=0;
    double gtDiameterNormVariance=0;
    double gtInertiaNormSum=0;
    double gtInertiaNormVariance=0;
    int gtNodesByDegree[NINTERACTIONSMAX];
    double gtNodesByDegreeSum[NINTERACTIONSMAX], gtNodesByDegreeVariance[NINTERACTIONSMAX];
    int gtNodesByDegreeLength;

    float DistrClusterSize[3];
    float DistrClusterSizeSum[NMAX][3];
    float DistrClusterSizeVariance[NMAX][3];
    int DistrClusterSizeCounter=0;

    int x[NMAX], x0, kMax;
    double dummydouble;
    int dummyint;



    NumberOfInputFiles=argc-2;

    for(j=0; j<3; j++)
    {
        DistrClusterSize[j]=0;
        for(i=0; i<NMAX; i++)
        {
            DistrClusterSizeSum[i][j]=0;
            DistrClusterSizeVariance[i][j]=0;
        }
    }
    for(j=0;j<2;j++)
    {
        for(i=0; i<NINTERACTIONSMAX; i++)
        {
            BondCountersSum[i][j]=0;
            BondCountersVariance[i][j]=0;
        }
    }
    for(i=0; i<NINTERACTIONSMAX; i++)
    {
        gtNodesByDegreeSum[i]=0;
        gtNodesByDegreeVariance[i]=0;
    }


    // Read input

    OutputFile=fopen(argv[argc-1],"wb+");

    printf("\nInputFiles: %d\n", NumberOfInputFiles);
    fprintf(OutputFile, "Data analysis program\nClusterAnalysisStatistics1.1.exe\n\n");
    fprintf(OutputFile, "InputFiles: %d\n", NumberOfInputFiles);

    for(i=0; i<NumberOfInputFiles; i++)
    {
        InputFiles[i] = fopen(argv[i+1], "r");
        printf("%s\n",argv[i+1]);

        while(fgets(string[i],strlength,InputFiles[i])!=NULL)
        {

            if(strcmp(string[i], "InputFile\n")==0)
            {
                fgets(string[i], strlength, InputFiles[i]);
                fputs(argv[i+1], OutputFile);
            }

            if(strncmp(string[i], "NumberOfClusters ", strlen("NumberOfClusters "))==0)
            {
                strcpy(NumberOfClustersString, string[i]);
                fgets(string[i], strlength, InputFiles[i]);
                sscanf(string[i], "%d\n", &NumberOfClusters);
                NumberOfClustersSum+=NumberOfClusters;
                NumberOfClustersVariance+=SQR(NumberOfClusters);
                NumberOfClustersCounter++;
            }

            if(strncmp(string[i], "NumberOfClustersBigger ", strlen("NumberOfClustersBigger "))==0)
            {
                strcpy(NumberOfClustersBiggerString, string[i]);
                fgets(string[i], strlength, InputFiles[i]);
                sscanf(string[i], "%d\n", &NumberOfClustersBigger);
                NumberOfClustersBiggerSum+=NumberOfClustersBigger;
                NumberOfClustersBiggerVariance+=SQR(NumberOfClustersBigger);
                NumberOfClustersBiggerCounter++;
            }

            if(strncmp(string[i], "NumberOfAClusters ", strlen("NumberOfAClusters "))==0)
            {
                strcpy(NumberOfAClustersString, string[i]);
                fgets(string[i], strlength, InputFiles[i]);
                sscanf(string[i], "%d\n", &NumberOfAClusters);
                NumberOfAClustersSum+=NumberOfAClusters;
                NumberOfAClustersVariance+=SQR(NumberOfAClusters);
                NumberOfAClustersCounter++;
            }


            if(strncmp(string[i], "AvgClusterSize ", strlen("AvgClusterSize "))==0)
            {
                strcpy(AvgClusterSizeString, string[i]);
                fgets(string[i], strlength, InputFiles[i]);
                sscanf(string[i], "%f\n", &AvgClusterSize);
                AvgClusterSizeSum+=AvgClusterSize;
                AvgClusterSizeVariance+=SQR(AvgClusterSize);
                AvgClusterSizeCounter++;
            }

            if(strncmp(string[i], "AvgAClusterSize ", strlen("AvgAClusterSize "))==0)
            {
                strcpy(AvgAClusterSizeString, string[i]);
                fgets(string[i], strlength, InputFiles[i]);
                sscanf(string[i], "%f\n", &AvgAClusterSize);
                AvgAClusterSizeSum+=AvgAClusterSize;
                AvgAClusterSizeVariance+=SQR(AvgAClusterSize);
                AvgAClusterSizeCounter++;
            }

            if(strncmp(string[i], "AvgBClusterSize ", strlen("AvgBClusterSize "))==0)
            {
                strcpy(AvgBClusterSizeString, string[i]);
                fgets(string[i], strlength, InputFiles[i]);
                sscanf(string[i], "%f\n", &AvgBClusterSize);
                AvgBClusterSizeSum+=AvgBClusterSize;
                AvgBClusterSizeVariance+=SQR(AvgBClusterSize);
                AvgBClusterSizeCounter++;
            }


            if(strncmp(string[i], "AvgClusterSizeBigger ", strlen("AvgClusterSizeBigger "))==0)
            {
                strcpy(AvgClusterSizeBiggerString, string[i]);
                fgets(string[i], strlength, InputFiles[i]);
                sscanf(string[i], "%f\n", &AvgClusterSizeBigger);
                AvgClusterSizeBiggerSum+=AvgClusterSizeBigger;
                AvgClusterSizeBiggerVariance+=SQR(AvgClusterSizeBigger);
                AvgClusterSizeBiggerCounter++;
            }

            if(strncmp(string[i], "AvgAClusterSizeBigger ", strlen("AvgAClusterSizeBigger "))==0)
            {
                strcpy(AvgAClusterSizeBiggerString, string[i]);
                fgets(string[i], strlength, InputFiles[i]);
                sscanf(string[i], "%f\n", &AvgAClusterSizeBigger);
                AvgAClusterSizeBiggerSum+=AvgAClusterSizeBigger;
                AvgAClusterSizeBiggerVariance+=SQR(AvgAClusterSizeBigger);
                AvgAClusterSizeBiggerCounter++;
            }

            if(strncmp(string[i], "AvgBClusterSizeBigger ", strlen("AvgBClusterSizeBigger "))==0)
            {
                strcpy(AvgBClusterSizeBiggerString, string[i]);
                fgets(string[i], strlength, InputFiles[i]);
                sscanf(string[i], "%f\n", &AvgBClusterSizeBigger);
                AvgBClusterSizeBiggerSum+=AvgBClusterSizeBigger;
                AvgBClusterSizeBiggerVariance+=SQR(AvgBClusterSizeBigger);
                AvgBClusterSizeBiggerCounter++;
            }


            if(strncmp(string[i], "AvgAConnectivity ", strlen("AvgAConnectivity "))==0)
            {
                strcpy(AvgAConnectivityString, string[i]);
                fgets(string[i], strlength, InputFiles[i]);
                sscanf(string[i], "%f\n", &AvgAConnectivity);
                AvgAConnectivitySum+=AvgAConnectivity;
                AvgAConnectivityVariance+=SQR(AvgAConnectivity);
                AvgAConnectivityCounter++;
            }

            if(strncmp(string[i], "AvgBConnectivity ", strlen("AvgBConnectivity "))==0)
            {
                strcpy(AvgBConnectivityString, string[i]);
                fgets(string[i], strlength, InputFiles[i]);
                sscanf(string[i], "%f\n", &AvgBConnectivity);
                AvgBConnectivitySum+=AvgBConnectivity;
                AvgBConnectivityVariance+=SQR(AvgBConnectivity);
                AvgBConnectivityCounter++;
            }


            if(strncmp(string[i], "CoalescenceLikelihood1 ", strlen("CoalescenceLikelihood1 "))==0)
            {
                strcpy(CoalescenceLikelihood1String, string[i]);
                fgets(string[i], strlength, InputFiles[i]);
                sscanf(string[i], "%f\n", &CoalescenceLikelihood1);
                CoalescenceLikelihood1Sum+=CoalescenceLikelihood1;
                CoalescenceLikelihood1Variance+=SQR(CoalescenceLikelihood1);
                CoalescenceLikelihood1Counter++;
            }

            if(strncmp(string[i], "CoalescenceLikelihood2 ", strlen("CoalescenceLikelihood2 "))==0)
            {
                strcpy(CoalescenceLikelihood2String, string[i]);
                fgets(string[i], strlength, InputFiles[i]);
                sscanf(string[i], "%f\n", &CoalescenceLikelihood2);
                CoalescenceLikelihood2Sum+=CoalescenceLikelihood2;
                CoalescenceLikelihood2Variance+=SQR(CoalescenceLikelihood2);
                CoalescenceLikelihood2Counter++;
            }


            if(strncmp(string[i], "FractionOfAInClustersOfSizeUpTo1 ", strlen("FractionOfAInClustersOfSizeUpTo1 "))==0)
            {
                strcpy(FractionOfAInClustersOfSizeUpTo1String, string[i]);
                fgets(string[i], strlength, InputFiles[i]);
                sscanf(string[i], "%f\n", &FractionOfAInClustersOfSizeUpTo1);
                FractionOfAInClustersOfSizeUpTo1Sum+=FractionOfAInClustersOfSizeUpTo1;
                FractionOfAInClustersOfSizeUpTo1Variance+=SQR(FractionOfAInClustersOfSizeUpTo1);
                FractionOfAInClustersOfSizeUpTo1Counter++;
            }

            if(strncmp(string[i], "FractionOfAInClustersOfSizeUpTo5 ", strlen("FractionOfAInClustersOfSizeUpTo5 "))==0)
            {
                strcpy(FractionOfAInClustersOfSizeUpTo5String, string[i]);
                fgets(string[i], strlength, InputFiles[i]);
                sscanf(string[i], "%f\n", &FractionOfAInClustersOfSizeUpTo5);
                FractionOfAInClustersOfSizeUpTo5Sum+=FractionOfAInClustersOfSizeUpTo5;
                FractionOfAInClustersOfSizeUpTo5Variance+=SQR(FractionOfAInClustersOfSizeUpTo5);
                FractionOfAInClustersOfSizeUpTo5Counter++;
            }

            if(strncmp(string[i], "FractionOfAInClustersOfSizeUpTo10 ", strlen("FractionOfAInClustersOfSizeUpTo10 "))==0)
            {
                strcpy(FractionOfAInClustersOfSizeUpTo10String, string[i]);
                fgets(string[i], strlength, InputFiles[i]);
                sscanf(string[i], "%f\n", &FractionOfAInClustersOfSizeUpTo10);
                FractionOfAInClustersOfSizeUpTo10Sum+=FractionOfAInClustersOfSizeUpTo10;
                FractionOfAInClustersOfSizeUpTo10Variance+=SQR(FractionOfAInClustersOfSizeUpTo10);
                FractionOfAInClustersOfSizeUpTo10Counter++;
            }

            if(strncmp(string[i], "FractionOfFreeB ", strlen("FractionOfFreeB "))==0)
            {
                strcpy(FractionOfFreeBString, string[i]);
                fgets(string[i], strlength, InputFiles[i]);
                sscanf(string[i], "%f\n", &FractionOfFreeB);
                FractionOfFreeBSum+=FractionOfFreeB;
                FractionOfFreeBVariance+=SQR(FractionOfFreeB);
                FractionOfFreeBCounter++;
            }

            if(strncmp(string[i], "BondCounter ", strlen("BondCounter "))==0)
            {
                strcpy(BondCountersString, string[i]);
                BondCountersCounter++;
                k=0;
                while( fgets(string[i],strlength,InputFiles[i])!=NULL && strcmp(string[i],"\n")!=0 )
                {
                    if(BondCountersCounter==1)
                        sscanf(string[i], "%s %f %f\n", BondTypeString[k], &BondCounters[0], &BondCounters[1]);
                    else
                    {
                        sscanf(string[i], "%s %f %f\n", DummyString, &BondCounters[0], &BondCounters[1]);
                        if( strcmp(DummyString, BondTypeString[k])!=0 )
                        {
                            printf("!!!!!!!! ERROR: Incompatible BondTypeString\n%s %s", DummyString, BondTypeString[k]);
                            return 3;
                        }
                    }
                    BondCountersSum[k][0]+=BondCounters[0];
                    BondCountersSum[k][1]+=BondCounters[1];
                    BondCountersVariance[k][0]+=SQR(BondCounters[0]);
                    BondCountersVariance[k][1]+=SQR(BondCounters[1]);
                    k++;
                }
                NumberOfInteractions=k;
            }


            if(strncmp(string[i], "Clusters detail", strlen("Clusters detail"))==0)
            {
                strcpy(ClustersInertiaString, string[i]);
                fgets(string[i],strlength,InputFiles[i]);
                strcat(ClustersInertiaString, string[i]);
                fgets(string[i],strlength,InputFiles[i]);
                strcat(ClustersInertiaString, string[i]);
                fgets(string[i],strlength,InputFiles[i]);
                strcat(ClustersInertiaString, string[i]);
                while( fgets(string[i],strlength,InputFiles[i])!=NULL && strcmp(string[i],"\n")!=0 )
                {
                    if(sscanf(string[i], "%d %*s\n", &ClusterSize)!=1)
                        printf("*** ERROR sscanf 1a - Clustersize = %d\n", ClusterSize);
                    if(ClusterSize>=ClusterSizeThreasholdForInertia)
                    {
                        if(sscanf(string[i], "%d %*d %*d %d %d %*g %*g %*g %*g %*g %g %g\n", &dummyint, &EdgesInCluster, &NumberOfTerminalNodes, &Compactness, &Roundness)!=5)
                            printf("*** ERROR sscanf 1b - '%d %d %d %g %g'\n", ClusterSize, EdgesInCluster, NumberOfTerminalNodes, Compactness, Roundness);
                        if(dummyint!=ClusterSize)
                            printf("*** ERROR sscanf 1c - '%d %d %d %g %g'\n", ClusterSize, EdgesInCluster, NumberOfTerminalNodes, Compactness, Roundness);
                        ClustersInertiaCounter++;
                        EdgesInClusterSum+=EdgesInCluster*1.0/ClusterSize;
                        EdgesInClusterVariance+=SQR(EdgesInCluster*1.0/ClusterSize);
                        NumberOfTerminalNodesSum+=NumberOfTerminalNodes*1.0/ClusterSize;
                        NumberOfTerminalNodesVariance+=SQR(NumberOfTerminalNodes*1.0/ClusterSize);
                        CompactnessSum+=Compactness;
                        CompactnessVariance+=SQR(Compactness);
                        RoundnessSum+=Roundness;
                        RoundnessVariance+=SQR(Roundness);
                    }
                    fgets(string[i],strlength,InputFiles[i]);
                    if(ClusterSize>=ClusterSizeThreasholdForInertia)
                    {
                        //if(sscanf(string[i], " > %*d %f %*d %f\n", &gtInertiaNorm, &gtDiameterNorm)!=2)
                        if(sscanf(string[i], " > %d\n", &NumberOfUnboundNodes)!=1)
                            printf("*** ERROR sscanf 2\n");
                        NumberOfUnboundNodesSum+=NumberOfUnboundNodes*1.0/ClusterSize;
                        NumberOfUnboundNodesVariance+=SQR(NumberOfUnboundNodes*1.0/ClusterSize);
                    }
                    fgets(string[i],strlength,InputFiles[i]);
                    if(ClusterSize>=ClusterSizeThreasholdForInertia)
                    {
                        if(strncmp(" >> ", string[i], 4)==0)
                            sscanf(string[i], " >> %[^\t]", DummyString);
                            //*DummyString=*string[i]+4;
                        else
                            printf("*** ERROR sscanf 3\n");
                        k=0;
                        token = strtok(DummyString, " ");
                        //printf("%s, %s, %s\n", string[i], DummyString, token);
                        while( token != NULL )
                        {
                            //printf("%s %s %s\n", string[i], DummyString, token);
                            gtNodesByDegree[k] = (int) strtol(token, &p, 10);
                            dummydouble=gtNodesByDegree[k]*1.0/ClusterSize;
                            gtNodesByDegreeSum[k]+=dummydouble;
                            gtNodesByDegreeVariance[k]+=SQR(dummydouble);
                            token = strtok(NULL, " ");
                            k++;
                        }
                        gtNodesByDegreeLength=k;
                    }
                }
            }


            if(strncmp(string[i], "DistrClusterSizeAll", strlen("DistrClusterSizeAll"))==0)
            {
                strcpy(DistrClusterSizeString, string[i]);
                fgets(string[i],strlength,InputFiles[i]);
                DistrClusterSizeCounter++;
                k=0;
                while( fgets(string[i],strlength,InputFiles[i])!=NULL && strcmp(string[i],"\n")!=0 )
                {
                    sscanf(string[i], "%d %f %f %f\n", &x0, &DistrClusterSize[0], &DistrClusterSize[1], &DistrClusterSize[2]);
                    if(DistrClusterSizeCounter==1)
                        x[k]=x0;
                    else
                    {
                        if(x0!=x[k])
                            return 2;
                    }
                    DistrClusterSizeSum[k][0]+=DistrClusterSize[0];
                    DistrClusterSizeSum[k][1]+=DistrClusterSize[1];
                    DistrClusterSizeSum[k][2]+=DistrClusterSize[2];
                    DistrClusterSizeVariance[k][0]+=SQR(DistrClusterSize[0]);
                    DistrClusterSizeVariance[k][1]+=SQR(DistrClusterSize[1]);
                    DistrClusterSizeVariance[k][2]+=SQR(DistrClusterSize[2]);
                    k++;
                }
                kMax=k;
            }
        }
        fclose(InputFiles[i]);
    }



    // Print output

    printf("OutputFile:\n%s\n", argv[argc-1]);

    fprintf(OutputFile,"\n\n");
    NumberOfClustersSum/=NumberOfClustersCounter;
    NumberOfClustersVariance=(NumberOfClustersVariance-NumberOfClustersCounter*SQR(NumberOfClustersSum))/(NumberOfClustersCounter-1.0);
    fputs(NumberOfClustersString,OutputFile);
    fprintf(OutputFile, "%f %+f\n", NumberOfClustersSum, sqrt(NumberOfClustersVariance));

    NumberOfClustersBiggerSum/=NumberOfClustersBiggerCounter;
    NumberOfClustersBiggerVariance=(NumberOfClustersBiggerVariance-NumberOfClustersBiggerCounter*SQR(NumberOfClustersBiggerSum))/(NumberOfClustersBiggerCounter-1.0);
    fputs(NumberOfClustersBiggerString,OutputFile);
    fprintf(OutputFile, "%f %+f\n", NumberOfClustersBiggerSum, sqrt(NumberOfClustersBiggerVariance));

    NumberOfAClustersSum/=NumberOfAClustersCounter;
    NumberOfAClustersVariance=(NumberOfAClustersVariance-NumberOfAClustersCounter*SQR(NumberOfAClustersSum))/(NumberOfAClustersCounter-1.0);
    fputs(NumberOfAClustersString,OutputFile);
    fprintf(OutputFile, "%f %+f\n\n", NumberOfAClustersSum, sqrt(NumberOfAClustersVariance));

    AvgClusterSizeSum/=AvgClusterSizeCounter;
    AvgClusterSizeVariance=(AvgClusterSizeVariance-AvgClusterSizeCounter*SQR(AvgClusterSizeSum))/(AvgClusterSizeCounter-1.0);
    fputs(AvgClusterSizeString,OutputFile);
    fprintf(OutputFile, "%f %+f\n", AvgClusterSizeSum, sqrt(AvgClusterSizeVariance));

    AvgAClusterSizeSum/=AvgAClusterSizeCounter;
    AvgAClusterSizeVariance=(AvgAClusterSizeVariance-AvgAClusterSizeCounter*SQR(AvgAClusterSizeSum))/(AvgAClusterSizeCounter-1.0);
    fputs(AvgAClusterSizeString,OutputFile);
    fprintf(OutputFile, "%f %+f\n", AvgAClusterSizeSum, sqrt(AvgAClusterSizeVariance));

    AvgBClusterSizeSum/=AvgBClusterSizeCounter;
    AvgBClusterSizeVariance=(AvgBClusterSizeVariance-AvgBClusterSizeCounter*SQR(AvgBClusterSizeSum))/(AvgBClusterSizeCounter-1.0);
    fputs(AvgBClusterSizeString,OutputFile);
    fprintf(OutputFile, "%f %+f\n\n", AvgBClusterSizeSum, sqrt(AvgBClusterSizeVariance));

    AvgClusterSizeBiggerSum/=AvgClusterSizeBiggerCounter;
    AvgClusterSizeBiggerVariance=(AvgClusterSizeBiggerVariance-AvgClusterSizeBiggerCounter*SQR(AvgClusterSizeBiggerSum))/(AvgClusterSizeBiggerCounter-1.0);
    fputs(AvgClusterSizeBiggerString,OutputFile);
    fprintf(OutputFile, "%f %+f\n", AvgClusterSizeBiggerSum, sqrt(AvgClusterSizeBiggerVariance));

    AvgAClusterSizeBiggerSum/=AvgAClusterSizeBiggerCounter;
    AvgAClusterSizeBiggerVariance=(AvgAClusterSizeBiggerVariance-AvgAClusterSizeBiggerCounter*SQR(AvgAClusterSizeBiggerSum))/(AvgAClusterSizeBiggerCounter-1.0);
    fputs(AvgAClusterSizeBiggerString,OutputFile);
    fprintf(OutputFile, "%f %+f\n", AvgAClusterSizeBiggerSum, sqrt(AvgAClusterSizeBiggerVariance));

    AvgBClusterSizeBiggerSum/=AvgBClusterSizeBiggerCounter;
    AvgBClusterSizeBiggerVariance=(AvgBClusterSizeBiggerVariance-AvgBClusterSizeBiggerCounter*SQR(AvgBClusterSizeBiggerSum))/(AvgBClusterSizeBiggerCounter-1.0);
    fputs(AvgBClusterSizeBiggerString,OutputFile);
    fprintf(OutputFile, "%f %+f\n\n", AvgBClusterSizeBiggerSum, sqrt(AvgBClusterSizeBiggerVariance));

    AvgAConnectivitySum/=AvgAConnectivityCounter;
    AvgAConnectivityVariance=(AvgAConnectivityVariance-AvgAConnectivityCounter*SQR(AvgAConnectivitySum))/(AvgAConnectivityCounter-1.0);
    fputs(AvgAConnectivityString,OutputFile);
    fprintf(OutputFile, "%f %+f\n", AvgAConnectivitySum, sqrt(AvgAConnectivityVariance));

    AvgBConnectivitySum/=AvgBConnectivityCounter;
    AvgBConnectivityVariance=(AvgBConnectivityVariance-AvgBConnectivityCounter*SQR(AvgBConnectivitySum))/(AvgBConnectivityCounter-1.0);
    if(strcmp(AvgBConnectivityString,"\0")==0)
        fprintf(OutputFile, "AvgBConnectivity not found\n0 +0\n\n");
    else
    {
        fputs(AvgBConnectivityString,OutputFile);
        fprintf(OutputFile, "%f %+f\n\n", AvgBConnectivitySum, sqrt(AvgBConnectivityVariance));
    }

    CoalescenceLikelihood1Sum/=CoalescenceLikelihood1Counter;
    CoalescenceLikelihood1Variance=(CoalescenceLikelihood1Variance-CoalescenceLikelihood1Counter*SQR(CoalescenceLikelihood1Sum))/(CoalescenceLikelihood1Counter-1.0);
    if(strcmp(CoalescenceLikelihood1String,"\0")==0)
        fprintf(OutputFile, "CoalescenceLikelihood1 not found\n0 +0\n");
    else
    {
        fputs(CoalescenceLikelihood1String,OutputFile);
        fprintf(OutputFile, "%f %+f\n", CoalescenceLikelihood1Sum, sqrt(CoalescenceLikelihood1Variance));
    }

    CoalescenceLikelihood2Sum/=CoalescenceLikelihood2Counter;
    CoalescenceLikelihood2Variance=(CoalescenceLikelihood2Variance-CoalescenceLikelihood2Counter*SQR(CoalescenceLikelihood2Sum))/(CoalescenceLikelihood2Counter-1.0);
    if(strcmp(CoalescenceLikelihood2String,"\0")==0)
        fprintf(OutputFile, "CoalescenceLikelihood2 not found\n0 +0\n\n");
    else
    {
        fputs(CoalescenceLikelihood2String,OutputFile);
        fprintf(OutputFile, "%f %+f\n\n", CoalescenceLikelihood2Sum, sqrt(CoalescenceLikelihood2Variance));
    }

    FractionOfAInClustersOfSizeUpTo1Sum/=FractionOfAInClustersOfSizeUpTo1Counter;
    FractionOfAInClustersOfSizeUpTo1Variance=(FractionOfAInClustersOfSizeUpTo1Variance-FractionOfAInClustersOfSizeUpTo1Counter*SQR(FractionOfAInClustersOfSizeUpTo1Sum))/(FractionOfAInClustersOfSizeUpTo1Counter-1.0);
    if(strcmp(FractionOfAInClustersOfSizeUpTo1String,"\0")==0)
        fprintf(OutputFile, "FractionOfAInClustersOfSizeUpTo1 not found\n-1 +1\n");
    else
    {
        fputs(FractionOfAInClustersOfSizeUpTo1String,OutputFile);
        fprintf(OutputFile, "%f %+f\n", FractionOfAInClustersOfSizeUpTo1Sum, sqrt(FractionOfAInClustersOfSizeUpTo1Variance));
    }

    FractionOfAInClustersOfSizeUpTo5Sum/=FractionOfAInClustersOfSizeUpTo5Counter;
    FractionOfAInClustersOfSizeUpTo5Variance=(FractionOfAInClustersOfSizeUpTo5Variance-FractionOfAInClustersOfSizeUpTo5Counter*SQR(FractionOfAInClustersOfSizeUpTo5Sum))/(FractionOfAInClustersOfSizeUpTo5Counter-1.0);
    if(strcmp(FractionOfAInClustersOfSizeUpTo5String,"\0")==0)
        fprintf(OutputFile, "FractionOfAInClustersOfSizeUpTo5 not found\n-1 +1\n");
    else
    {
        fputs(FractionOfAInClustersOfSizeUpTo5String,OutputFile);
        fprintf(OutputFile, "%f %+f\n", FractionOfAInClustersOfSizeUpTo5Sum, sqrt(FractionOfAInClustersOfSizeUpTo5Variance));
    }

    FractionOfAInClustersOfSizeUpTo10Sum/=FractionOfAInClustersOfSizeUpTo10Counter;
    FractionOfAInClustersOfSizeUpTo10Variance=(FractionOfAInClustersOfSizeUpTo10Variance-FractionOfAInClustersOfSizeUpTo10Counter*SQR(FractionOfAInClustersOfSizeUpTo10Sum))/(FractionOfAInClustersOfSizeUpTo10Counter-1.0);
    if(strcmp(FractionOfAInClustersOfSizeUpTo10String,"\0")==0)
        fprintf(OutputFile, "FractionOfAInClustersOfSizeUpTo10 not found\n-1 +1\n\n");
    else
    {
        fputs(FractionOfAInClustersOfSizeUpTo10String,OutputFile);
        fprintf(OutputFile, "%f %+f\n\n", FractionOfAInClustersOfSizeUpTo10Sum, sqrt(FractionOfAInClustersOfSizeUpTo10Variance));
    }

    FractionOfFreeBSum/=FractionOfFreeBCounter;
    FractionOfFreeBVariance=(FractionOfFreeBVariance-FractionOfFreeBCounter*SQR(FractionOfFreeBSum))/(FractionOfFreeBCounter-1.0);
    if(strcmp(FractionOfFreeBString,"\0")==0)
        fprintf(OutputFile, "FractionOfFreeB not found\n-1 +1\n\n");
    else
    {
        fputs(FractionOfFreeBString,OutputFile);
        fprintf(OutputFile, "%f %+f\n\n", FractionOfFreeBSum, sqrt(FractionOfFreeBVariance));
    }

    if(strcmp(BondCountersString,"\0")==0)
        fprintf(OutputFile, "BondCounters not found\n\n");
    else
    {
        fputs(BondCountersString,OutputFile);
        for(k=0; k<NumberOfInteractions; k++)
        {
            for(j=0; j<=1; j++)
            {
                BondCountersSum[k][j]/=BondCountersCounter;
                BondCountersVariance[k][j]=(BondCountersVariance[k][j]-BondCountersCounter*SQR(BondCountersSum[k][j]))/(BondCountersCounter-1.0);
            }
            fprintf(OutputFile, "%s %f %+f %f %+f\n", BondTypeString[k], BondCountersSum[k][0], sqrt(BondCountersVariance[k][0]), BondCountersSum[k][1], sqrt(BondCountersVariance[k][1]));
        }
        fprintf(OutputFile, "\n");
    }


    EdgesInClusterSum/=ClustersInertiaCounter;
    EdgesInClusterVariance=(EdgesInClusterVariance-ClustersInertiaCounter*SQR(EdgesInClusterSum))/(ClustersInertiaCounter-1.0);
    NumberOfTerminalNodesSum/=ClustersInertiaCounter;
    NumberOfTerminalNodesVariance=(NumberOfTerminalNodesVariance-ClustersInertiaCounter*SQR(NumberOfTerminalNodesSum))/(ClustersInertiaCounter-1.0);
    CompactnessSum/=ClustersInertiaCounter;
    CompactnessVariance=(CompactnessVariance-ClustersInertiaCounter*SQR(CompactnessSum))/(ClustersInertiaCounter-1.0);
    RoundnessSum/=ClustersInertiaCounter;
    RoundnessVariance=(RoundnessVariance-ClustersInertiaCounter*SQR(RoundnessSum))/(ClustersInertiaCounter-1.0);
    NumberOfUnboundNodesSum/=ClustersInertiaCounter;
    NumberOfUnboundNodesVariance=(NumberOfUnboundNodesVariance-ClustersInertiaCounter*SQR(NumberOfUnboundNodesSum))/(ClustersInertiaCounter-1.0);
    //gtInertiaNormSum/=ClustersInertiaCounter;
    //gtInertiaNormVariance=(gtInertiaNormVariance-ClustersInertiaCounter*SQR(gtInertiaNormSum))/(ClustersInertiaCounter-1.0);
    //gtDiameterNormSum/=ClustersInertiaCounter;
    //gtDiameterNormVariance=(gtDiameterNormVariance-ClustersInertiaCounter*SQR(gtDiameterNormSum))/(ClustersInertiaCounter-1.0);
    for(k=0; k<gtNodesByDegreeLength; k++)
    {
        gtNodesByDegreeSum[k]/=ClustersInertiaCounter;
        gtNodesByDegreeVariance[k]=(gtNodesByDegreeVariance[k]-ClustersInertiaCounter*SQR(gtNodesByDegreeSum[k]))/(ClustersInertiaCounter-1.0);
    }

    if(strcmp(ClustersInertiaString,"\0")==0)
        fprintf(OutputFile, "Clusters details not found\n-1 +1\n\n");
    else
    {
        fprintf(OutputFile, "Clusters details (based on %d clusters of size >%d)\nEdgesInCluster per node\n", ClustersInertiaCounter, ClusterSizeThreasholdForInertia);
        fprintf(OutputFile, "%g %+g\n", EdgesInClusterSum, sqrt(EdgesInClusterVariance));
        fputs("NumberOfTerminalNodes [%]\n", OutputFile);
        fprintf(OutputFile, "%g %+g\n", NumberOfTerminalNodesSum, sqrt(NumberOfTerminalNodesVariance));
        fputs("Compactness\n", OutputFile);
        fprintf(OutputFile, "%g %+g\n", CompactnessSum, sqrt(CompactnessVariance));
        fputs("Roundness\n", OutputFile);
        fprintf(OutputFile, "%g %+g\n", RoundnessSum, sqrt(RoundnessVariance));
        fputs("NumberOfUnboundNodes [%]\n", OutputFile);
        fprintf(OutputFile, "%g %+g\n", NumberOfUnboundNodesSum, sqrt(NumberOfUnboundNodesVariance));
        //fputs("gtInertia/ClusterSize[][0]^2\n", OutputFile);
        //fprintf(OutputFile, "%g %+g\n", gtInertiaNormSum, sqrt(gtInertiaNormVariance));
        //fputs("gtDiameter/gtMooreBound\n", OutputFile);
        //fprintf(OutputFile, "%g %+g\n", gtDiameterNormSum, sqrt(gtDiameterNormVariance));
        fputs("Adeg0, Adeg1, Adeg2, Adeg3, Adeg4, Bdeg0, Bdeg1, Bdeg2, Bdeg3, Bdeg4\n", OutputFile);
        dummydouble=0;
        for(k=0; k<14; k++)
        {
            fprintf(OutputFile, "%g %+g ", gtNodesByDegreeSum[k], sqrt(gtNodesByDegreeVariance[k]));
            dummydouble+=gtNodesByDegreeSum[k];
        }
        fprintf(OutputFile, "   %g\n\n", dummydouble);
    }


    fputs(DistrClusterSizeString,OutputFile);
    for(k=0; k<kMax; k++)
    {
        DistrClusterSizeSum[k][0]/=DistrClusterSizeCounter;
        DistrClusterSizeVariance[k][0]=(DistrClusterSizeVariance[k][0]-DistrClusterSizeCounter*SQR(DistrClusterSizeSum[k][0]))/(DistrClusterSizeCounter-1.0);

        DistrClusterSizeSum[k][1]/=DistrClusterSizeCounter;
        DistrClusterSizeVariance[k][1]=(DistrClusterSizeVariance[k][1]-DistrClusterSizeCounter*SQR(DistrClusterSizeSum[k][1]))/(DistrClusterSizeCounter-1.0);

        DistrClusterSizeSum[k][2]/=DistrClusterSizeCounter;
        DistrClusterSizeVariance[k][2]=(DistrClusterSizeVariance[k][2]-DistrClusterSizeCounter*SQR(DistrClusterSizeSum[k][2]))/(DistrClusterSizeCounter-1.0);

        fprintf(OutputFile, "%d %f %+f %f %+f %f %+f\n", x[k], DistrClusterSizeSum[k][0], sqrt(DistrClusterSizeVariance[k][0]), DistrClusterSizeSum[k][1], sqrt(DistrClusterSizeVariance[k][1]), DistrClusterSizeSum[k][2], sqrt(DistrClusterSizeVariance[k][2]));
    }

    fclose(OutputFile);


    return 0;

}











