#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "initial_curves.h"
#include "curve_similarity.h"
#include "initialization.h"

int myRandom(int size, int k)
{
    int i, n, curvePos;

    static int numNums = 0;
    static int timesCalled = 0;
    static int *numArr = NULL;

    timesCalled++;

    if (size >= 0)
    {
        if (numArr != NULL)
            free (numArr);
        if ((numArr = malloc (sizeof(int) * size)) == NULL)
            return ERR_NO_MEM;
        for (i = 0; i  < size; i++)
            numArr[i] = i;
        numNums = size;
    }

    if (numNums == 0)
       return ERR_NO_NUM;

    n = rand() % numNums;
    i = numArr[n];
    numArr[n] = numArr[numNums-1];
    numNums--;

    if ( timesCalled == k )
    {
        free (numArr);
        numArr = 0;
    }

    return i;
}

int* randomSelection(curves* crvs, int totalNumOfCurves, int k)
{
	int i, count, curvePos;
	int *kcenters;

	kcenters = malloc( k * sizeof(int) );

	curvePos = myRandom(totalNumOfCurves, k);
	count = 0;
	kcenters[count] = curvePos;
	count++;

    while (count < k) 
    {
        curvePos = myRandom(-1, k);
        kcenters[count] = curvePos;
        count++;
    }
    return kcenters;
}

int binarySearch(pSum* partialSumsArray, double randValue, int totalNumOfCurves, int currentNumOfCentroids)
{
    int first, last, half;

    first = 0;
    last = totalNumOfCurves - currentNumOfCentroids - 1;
    half = (first+last)/2;

    while( first <= last )
    {
        if( randValue > partialSumsArray[half].partialSum)
        {
            if(randValue < partialSumsArray[half+1].partialSum)
                return partialSumsArray[half+1].originalPos;
            else
                first = half + 1; 
        }
        else if( randValue < partialSumsArray[half].partialSum)
        {
        	if( half == first )
        		return partialSumsArray[half].originalPos;
        	else
        	{	
	            if(randValue > partialSumsArray[half-1].partialSum)
	                return partialSumsArray[half].originalPos;
	            else
	                last = half - 1; 
	        }
        }
        else
            return partialSumsArray[half].originalPos;

        half = (first+last)/2;
    }
}

int* kmeansplusplus(curves* crvs, int totalNumOfCurves, int k, char* metric)
{
    int firstCentroidPos, currentNumOfCentroids, i, j, centroidPos, newCentroidPos, nextAvailablePosPartial;
    int *kcentroids;
    double minDistance, distance, sum, randValue;
    char *DTW="DTW", *Frechet="Frechet";
    plusplusInfo *plusplusArray, *pointBack;
    pSum *partialSumsArray;
    curveInfo* curveArray;

    curveArray = crvs -> curveArray;
    // arxikopoiiseis
    currentNumOfCentroids = 0;

    plusplusArray = malloc( totalNumOfCurves * sizeof(plusplusInfo) );
    for( i=0; i<totalNumOfCurves; i++ )
    {
        plusplusArray[i].minDistanceToCentroid = INFINITY;
        plusplusArray[i].minDistanceToCentroidSqrt = -1;
        plusplusArray[i].closestCentroid = -1;
        plusplusArray[i].isCentroid = 0;
    }

    kcentroids = malloc( k * sizeof(int) );
    //^

    /* dialegoume to proto centroid stin tixi */
    firstCentroidPos = rand() % totalNumOfCurves;

    plusplusArray[firstCentroidPos].minDistanceToCentroid = -1;
    plusplusArray[firstCentroidPos].minDistanceToCentroidSqrt = -1;
    plusplusArray[firstCentroidPos].closestCentroid = -1;
    plusplusArray[firstCentroidPos].isCentroid = 1;

    kcentroids[currentNumOfCentroids] = firstCentroidPos;

    currentNumOfCentroids++;
    //printf("KMEANSPP\n");
    while( currentNumOfCentroids != k )
    {
        partialSumsArray = malloc( (totalNumOfCurves-currentNumOfCentroids) * sizeof(pSum) );
        nextAvailablePosPartial = 0;
        sum = 0;

        for(i=0; i<totalNumOfCurves; i++)   //gia kathe mia kampili
        {
            if( plusplusArray[i].isCentroid == 0 )  //an i kampili den einai centroid
            {
                /* ipologismos distance kai sigkrisi me minDistance */
                if( strcmp(metric, Frechet) == 0 )
                    distance = discreteFrechetDistance( curveArray[i], curveArray[ kcentroids[currentNumOfCentroids-1] ]);
                else if( strcmp(metric, DTW) == 0 )
                    distance = dynamicTimeWarping( curveArray[i], curveArray[ kcentroids[currentNumOfCentroids-1] ]); 

                if( distance < plusplusArray[i].minDistanceToCentroid)
                {
                    plusplusArray[i].minDistanceToCentroid = distance;
                    plusplusArray[i].minDistanceToCentroidSqrt = distance * distance;
                    plusplusArray[i].closestCentroid = kcentroids[currentNumOfCentroids-1];
                    plusplusArray[i].isCentroid = 0;
                }
    
                sum += plusplusArray[i].minDistanceToCentroidSqrt;

                partialSumsArray[nextAvailablePosPartial].originalPos = i;
                partialSumsArray[nextAvailablePosPartial].partialSum = sum;
                nextAvailablePosPartial++;
            }

        }

        /* diadikasia epilogis kainouriou centroid */

        /* __rand sto diastima [0, sum] */
        randValue = ( (double)rand() / (double)RAND_MAX ) * sum;
        //printf("randValue: %f, sum: %f\n", randValue, sum);

        /* __binary search */
        newCentroidPos = binarySearch(partialSumsArray, randValue, totalNumOfCurves, currentNumOfCentroids); 

        /* __"metatropi" tis kampilis se centroid */
        plusplusArray[newCentroidPos].minDistanceToCentroid = -1;
        plusplusArray[newCentroidPos].minDistanceToCentroidSqrt = -1;
        plusplusArray[newCentroidPos].closestCentroid = -1;
        plusplusArray[newCentroidPos].isCentroid = 1;

        kcentroids[currentNumOfCentroids] = newCentroidPos;

        currentNumOfCentroids++;

        free(partialSumsArray);
        //printf("ananeosa pinaka kai evala kainourio centroid \n");
    }


    free(plusplusArray);
    return kcentroids;
} 