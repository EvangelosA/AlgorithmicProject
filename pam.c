#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "initial_curves.h"
#include "grid.h"
#include "hash.h"
#include "curve_similarity.h"
#include "assignment.h"
#include "pam.h"

double** createAndInitialize2D(int numOfCurves)
{
	int i, j;
	double **distanceStorage;

	distanceStorage = malloc(numOfCurves * sizeof(double*));

	for( i=0; i<numOfCurves; i++ ){
		distanceStorage[i] = malloc(numOfCurves * sizeof(double));

		for( j=0; j<numOfCurves; j++ )
			distanceStorage[i][j] = -1;
	}

	return distanceStorage;
}

void delete2D(double** distanceStorage, int numOfCurves)
{
	int i, j;

	for( i=0; i<numOfCurves; i++ ){
		free(distanceStorage[i]);
	}
	free(distanceStorage);

	return;
}

double calculateObjectiveFunction(curveInfo* curveArray, int* crvs, int numOfCurves, int exceptionID, double** distanceStorage)
{
	int i;
	double distance, objectiveFunctionValue;

	objectiveFunctionValue = 0;
	for( i=0; i<numOfCurves; i++ )
	{
		if( crvs[i] != exceptionID ) 
		{
			if( distanceStorage[exceptionID][crvs[i]] == -1 )
			{
				distance = dynamicTimeWarping( curveArray[ exceptionID ], curveArray[ crvs[i] ]);

				/*** UPDATE distanceStorage ***/
				distanceStorage[ exceptionID ][ crvs[i] ] = distance;
				distanceStorage[ crvs[i] ][ exceptionID ] = distance;
			}
			else
				distance = distanceStorage[exceptionID][crvs[i]];

			/*** UPDATE objectiveFunctionValue ***/
			objectiveFunctionValue += distance;
		}
	}

	return objectiveFunctionValue;
}


int* PartitionAroundMedoids(cluster *cls, int k, minDist *centroidsArray, int totalNumOfCurves, curveInfo* curveArray, double** distanceStorage)
{
	int i, j, l, m, curvesInClusterNum, exceptionID, exceptionIDPos, newCentroid;
	int *updatedClusters;
	double objectiveFunctionValue, objectiveFunctionValue_new, distance;

	updatedClusters = malloc( k * sizeof(int) );

	for(i=0; i<k; i++)
	{
		/*** IPOLOGISMOS OBJECTIVE FUNCTION ***/
		exceptionID = cls[i].centroidId;
		newCentroid = cls[i].centroidId;
		//printf("calculating obj funct 1\n");
		objectiveFunctionValue = calculateObjectiveFunction(curveArray, cls[i].crvs, cls[i].numOfCurves, exceptionID, distanceStorage);
		//printf("done\n" );

		/*** EPILOGI IPOPSIFION CENTROIDS ***/
		curvesInClusterNum = cls[i].numOfCurves;
		for(j=0; j<curvesInClusterNum; j++)
		{
			//printf("j:%d \n", j);
			/*** IPOLOGISMOS OBJECTIVE FUNCTION GIA TO NEO CENTROID STI THESI j ***/
			exceptionID = cls[i].crvs[j];
		 	objectiveFunctionValue_new = calculateObjectiveFunction(curveArray, cls[i].crvs, cls[i].numOfCurves, exceptionID, distanceStorage);

		 	if( objectiveFunctionValue_new < objectiveFunctionValue )
		 	{
		 		newCentroid = exceptionID;
		 		objectiveFunctionValue = objectiveFunctionValue_new;
		 	}
		}
		
		updatedClusters[i] = newCentroid;
	}

	return updatedClusters;
}

