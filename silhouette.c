#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "initial_curves.h"
#include "curve_similarity.h"
#include "grid.h"
#include "hash.h"
#include "assignment.h"
#include "silhouette.h"

void secondBestClusters(cluster* cls, int k, minDist* centroidsArray, curveInfo* curveArray, curves* centroidCurves, int totalNumOfCurves, int metric)
{
	int i, j, centroid, secondBestCentroid, secondBestCentroidPos;
	double secondBestDistance, distance;
	curveInfo* curveArray_mean;

	if( metric == 1 )
	{
		curveArray_mean = centroidCurves->curveArray;

		for( i=0; i<totalNumOfCurves; i++ )
		{	
			centroid = centroidsArray[i].centroid;
			secondBestCentroid = -1;
			secondBestCentroidPos = -1;
			secondBestDistance = INFINITY;

			for( j=0; j<k; j++ )
			{
				if( centroid != j )
				{
					distance = discreteFrechetDistance(curveArray[i], curveArray_mean[j]);

					if(distance < secondBestDistance)
					{
						secondBestCentroid = j;
						secondBestCentroidPos = j;
						secondBestDistance = distance;
					}
				}
			}
			centroidsArray[i].second_centroid = secondBestCentroid;
			centroidsArray[i].second_centroid_Pos = secondBestCentroidPos;
			centroidsArray[i].second_distance = secondBestDistance;
		}


	}
	else if( metric == 2 )
	{
		for(i=0; i<totalNumOfCurves; i++)
		{
			centroid = centroidsArray[i].centroid;
			secondBestCentroid = -1;
			secondBestCentroidPos = -1;
			secondBestDistance = INFINITY;

			for(j=0; j<k; j++)
			{
				if( centroid != cls[j].centroidId )		//den exei noima na ksanaipologisoume to idio
				{
					distance = dynamicTimeWarping( curveArray[i], curveArray[cls[j].centroidId] );

					if(distance < secondBestDistance)
					{
						secondBestCentroid = cls[j].centroidId;
						secondBestCentroidPos = j;
						secondBestDistance = distance;
					}
				}
			}

			centroidsArray[i].second_centroid = secondBestCentroid;
			centroidsArray[i].second_centroid_Pos = secondBestCentroidPos;
			centroidsArray[i].second_distance = secondBestDistance;

		}
	}
}

double meandDistanceOfObject(curveInfo* curveArray, int* crvs, int curvesInClusterNum, int exceptionID, double** distanceStorage, int metric)
{
	int i;
	double distance, sum, MO;

	sum = 0;
	for( i=0; i<curvesInClusterNum; i++ )
	{
		if( crvs[i] != exceptionID )
		{
			if( distanceStorage[exceptionID][crvs[i]] == -1 )
			{
				if( metric == 1)
					distance = discreteFrechetDistance( curveArray[ exceptionID ], curveArray[ crvs[i] ]);
				else
					distance = dynamicTimeWarping( curveArray[ exceptionID ], curveArray[ crvs[i] ]);

				/*** UPDATE distanceStorage ***/
				distanceStorage[ exceptionID ][ crvs[i] ] = distance;
				distanceStorage[ crvs[i] ][ exceptionID ] = distance;
			}
			else
				distance = distanceStorage[exceptionID][crvs[i]];

			/*** UPDATE sum ***/
			sum += distance;
		}
	}
	
	MO = (double)(sum / curvesInClusterNum );	//mean of object i

	return MO;
}

double* silhouetteCalculation(cluster* cls, int k, minDist* centroidsArray, curveInfo* curveArray, curves *centroidCurves, int totalNumOfCurves, int metric, double** distanceStorage_a, double** distanceStorage_b)
{
	int i, j, l, curvesInClusterNum, exceptionID, secondBestCentroidPos;
	double meanDistance_a, meanDistance_b, maxDistance_ab, silhouette_object, silhouette_cluster, sum, overallSum, overallSilhouetteCoefficient;
	double *resultsArray;
	double **distanceStorage;

	/**** PROTA IPOLOGISE TA DEUTERA KALITERA CLUSTERS **/
	//printf("second best cluster calculation...\n");
	secondBestClusters(cls, k, centroidsArray, curveArray, centroidCurves, totalNumOfCurves, metric);
	//printf("done\n");

	resultsArray = malloc( (k+1) * sizeof(double) ); 	//+1 thesi sto telos gia to sinoliko silhouette

	/*** GIA KATHE CLUSTER ***/
	overallSum = 0;
	for(i=0; i<k; i++)
	{
		sum = 0;
		curvesInClusterNum = cls[i].numOfCurves;

		/*** IPOLOGISMOS SILHOUETTE GIA KATHE OBJECT TOU CLUSTER ***/
		for( j=0; j<curvesInClusterNum; j++ )
		{
			/*** IPOLOGISMOS a(i), b(i) ***/
			exceptionID = cls[i].crvs[j];
		 	meanDistance_a = meandDistanceOfObject(curveArray, cls[i].crvs, cls[i].numOfCurves, exceptionID, distanceStorage_a, metric);	//a(i)

		 	secondBestCentroidPos = centroidsArray[exceptionID].second_centroid_Pos;
		 	meanDistance_b = meandDistanceOfObject(curveArray, cls[secondBestCentroidPos].crvs, cls[secondBestCentroidPos].numOfCurves, exceptionID, distanceStorage_b, metric);	//b(i)

		 	/*** s(i) = ( b(i)-a(i) ) / max{ a(i),b(i) } ***/
		 	maxDistance_ab = meanDistance_b;
			if( meanDistance_a > meanDistance_b )
				maxDistance_ab = meanDistance_a;

			silhouette_object = (meanDistance_b - meanDistance_a) / maxDistance_ab;		//s(i)

			/*** UPDATE sum(= ATHROISMA APO TA SILHOUETTES TON OBJECTS) ***/
			sum += silhouette_object;
		}

		/*** IPOLOGISMOS SILHOUETTE TOU CLUSTER ***/
		silhouette_cluster = (double)(sum / curvesInClusterNum);
		resultsArray[i] = silhouette_cluster;

		/*** UPDATE overallSum(= ATHROISMA APO SILHOUETTES OLON TON OBJECTS TOU DATASET) ***/
		overallSum += sum;

	}

	/*** IPOLOGISMOS OVERALL SILHOUETTE COEFFICIENT KAI TOPOTHETISI STIN TELEUTAIA THESI TOU PINAKA GIA i=k ***/
	overallSilhouetteCoefficient = (double)(overallSum / totalNumOfCurves);
	resultsArray[i] = overallSilhouetteCoefficient;

	return resultsArray;
}