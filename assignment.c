#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "initial_curves.h"
#include "curve_similarity.h"
#include "grid.h"
#include "hash.h"
#include "range.h"
#include "assignment.h"

minDist* loydsAssignment(curves* crvs, int totalNumOfCurves, int* centroids, int numOfCentroids, char* metric, curves* centroidCurves)
{
	int i, j, nearestCentroid;
	double minDistance, tempDist;
	minDist* distanceArray;
	if(centroidCurves==NULL)
	{
		distanceArray=(minDist*)malloc(totalNumOfCurves*sizeof(minDist));
		if( strcmp(metric,"Frechet")==0 )
		{
			for( i=0; i<totalNumOfCurves ; i++ )
			{
				minDistance=discreteFrechetDistance(crvs->curveArray[i],crvs->curveArray[centroids[0]]);
				nearestCentroid=centroids[0];
				for( j=1; j<numOfCentroids; j++)
				{
					tempDist=discreteFrechetDistance( crvs->curveArray[i] , crvs->curveArray[centroids[j]] );
					if( tempDist<minDistance )
					{
						minDistance=tempDist;
						nearestCentroid=centroids[j];
					}
				}
				distanceArray[i].distance=minDistance;
				distanceArray[i].centroid=nearestCentroid;
			}
		}
		else if( strcmp(metric,"DTW")==0 )
		{
			for( i=0; i<totalNumOfCurves ; i++ )
			{
				minDistance=dynamicTimeWarping(crvs->curveArray[i],crvs->curveArray[centroids[0]]);
				nearestCentroid=centroids[0];
				for( j=1; j<numOfCentroids; j++)
				{
					tempDist=dynamicTimeWarping( crvs->curveArray[i] , crvs->curveArray[centroids[j]] );
					if( tempDist<minDistance )
					{
						minDistance=tempDist;
						nearestCentroid=centroids[j];
					}
				}
				distanceArray[i].distance=minDistance;
				distanceArray[i].centroid=nearestCentroid;
			}
		}
		return distanceArray;
	}
	else
	{
		distanceArray=(minDist*)malloc(totalNumOfCurves*sizeof(minDist));
		if( strcmp(metric,"Frechet")==0 )
		{
			for( i=0; i<totalNumOfCurves ; i++ )
			{
				minDistance=discreteFrechetDistance(crvs->curveArray[i],centroidCurves->curveArray[0]);
				nearestCentroid=0;
				for( j=1; j<numOfCentroids; j++)
				{
					tempDist=discreteFrechetDistance( crvs->curveArray[i] , centroidCurves->curveArray[j] );
					if( tempDist<minDistance )
					{
						minDistance=tempDist;
						nearestCentroid=j;
					}
				}
				distanceArray[i].distance=minDistance;
				distanceArray[i].centroid=nearestCentroid;
			}
		}
		else if( strcmp(metric,"DTW")==0 )
		{
			for( i=0; i<totalNumOfCurves ; i++ )
			{
				minDistance=dynamicTimeWarping(crvs->curveArray[i],centroidCurves->curveArray[0]);
				nearestCentroid=0;
				for( j=1; j<numOfCentroids; j++)
				{
					tempDist=dynamicTimeWarping( crvs->curveArray[i] , centroidCurves->curveArray[j] );
					if( tempDist<minDistance )
					{
						minDistance=tempDist;
						nearestCentroid=j;
					}
				}
				distanceArray[i].distance=minDistance;
				distanceArray[i].centroid=nearestCentroid;
			}
		}
		return distanceArray;
	}
}

minDist* rangeAssignment(curves* crvs, int totalNumOfCurves, grids* grds, int* centroids, int numOfCentroids, char* metric, ls* hashStructures, int numOfHashStructures, double startR, int k, long long int m, int table_size, curves* centroidCurves, grids* centroidGrids)
{
	minDist* dist;
	int *check;
	int *rangeArray;
	int added, distType, i, j, centroids_added, jmin;
	double r,r2, distance, minDistance;
	if(centroidCurves==NULL)
	{
		check=malloc(sizeof(int)*totalNumOfCurves);
		dist=malloc(totalNumOfCurves*sizeof(minDist));
		for(i=0; i<totalNumOfCurves; i++)
		{
			check[i]=0;
		}
		centroids_added=0;
		r=startR;
		r2=0;
		if(strcmp(metric,"Frechet")==0)
		{
			distType=1;
		}
		else
		{
			distType=2;
		}
		while(4*centroids_added<3*numOfCentroids)
		{
			centroids_added=0;
			for(i=0; i<numOfCentroids; i++)
			{
				added=0;
				rangeArray=hash_structure_range(crvs, grds, crvs->curveArray[centroids[i]], grds->gridArray[centroids[i]], hashStructures, numOfHashStructures, r, r2, distType, k, m, table_size);
				for(j=0; j<totalNumOfCurves; j++)
				{
					if(check[j]==0)
					{
						check[j]=1;
						if(distType==1)
							dist[j].distance=discreteFrechetDistance(crvs->curveArray[centroids[i]],crvs->curveArray[j]);
						else
							dist[j].distance=dynamicTimeWarping(crvs->curveArray[centroids[i]],crvs->curveArray[j]);
						dist[j].centroid=centroids[i];
						added=1;
					}
					else if(check[j]==1)
					{
						if(distType==1)
							distance=discreteFrechetDistance(crvs->curveArray[centroids[i]],crvs->curveArray[j]);
						else
							distance=dynamicTimeWarping(crvs->curveArray[centroids[i]],crvs->curveArray[j]);
						if(distance<dist[j].distance)
						{
							dist[j].distance=distance;
							dist[j].centroid=centroids[i];
							added=1;
						}
					}
				}
				free(rangeArray);
				centroids_added=centroids_added+added;
			}
			for(j=0; j<totalNumOfCurves; j++)
			{
				if(check[j]==1)
					check[j]=2;
			}
			r2=r;
			r=2*r;
		}
		for(i=0; i<totalNumOfCurves; i++)
		{
			if(check[i]==0)
			{
				if(distType==1)
					distance=discreteFrechetDistance(crvs->curveArray[centroids[0]],crvs->curveArray[j]);
				else
					distance=dynamicTimeWarping(crvs->curveArray[centroids[0]],crvs->curveArray[j]);
				jmin=0;
				minDistance=distance;
				for(j=1; j<numOfCentroids; j++)
				{
					if(distType==1)
						distance=discreteFrechetDistance(crvs->curveArray[centroids[i]],crvs->curveArray[j]);
					else
						distance=dynamicTimeWarping(crvs->curveArray[centroids[i]],crvs->curveArray[j]);
					if(minDistance>distance)
					{
						jmin=j;
						minDistance=distance;
					}
				}
			}
		}
		free(check);
		return dist;
	}
	else
	{
		check=malloc(sizeof(int)*totalNumOfCurves);
		dist=malloc(totalNumOfCurves*sizeof(minDist));
		for(i=0; i<totalNumOfCurves; i++)
		{
			check[i]=0;
		}		centroids_added=0;
		r=startR;
		r2=0;
		if(strcmp(metric,"Frechet")==0)
		{
			distType=1;
		}
		else
		{
			distType=2;
		}
		while(4*centroids_added<3*numOfCentroids)
		{
			centroids_added=0;
			for(i=0; i<numOfCentroids; i++)
			{
				added=0;
				rangeArray=hash_structure_range(crvs, grds, centroidCurves->curveArray[i], centroidGrids->gridArray[i], hashStructures, numOfHashStructures, r, r2, distType, k, m, table_size);
				for(j=0; j<totalNumOfCurves; j++)
				{
					if(check[j]==0)
					{
						check[j]=1;
						if(distType==1)
							dist[j].distance=discreteFrechetDistance(centroidCurves->curveArray[i],crvs->curveArray[j]);
						else
							dist[j].distance=dynamicTimeWarping(centroidCurves->curveArray[i],crvs->curveArray[j]);
						dist[j].centroid=i;
						added=1;
					}
					else if(check[j]==1)
					{
						if(distType==1)
							distance=discreteFrechetDistance(centroidCurves->curveArray[i],crvs->curveArray[j]);
						else
							distance=dynamicTimeWarping(centroidCurves->curveArray[i],crvs->curveArray[j]);
						if(distance<dist[j].distance)
						{
							dist[j].distance=distance;
							dist[j].centroid=i;
							added=1;
						}
					}
				}
				free(rangeArray);
				centroids_added=centroids_added+added;
			}
			for(j=0; j<totalNumOfCurves; j++)
			{
				if(check[j]==1)
					check[j]=2;
			}
			r2=r;
			r=2*r;
		}
		for(i=0; i<totalNumOfCurves; i++)
		{
			if(check[i]==0)
			{
				if(distType==1)
					distance=discreteFrechetDistance(centroidCurves->curveArray[0],crvs->curveArray[j]);
				else
					distance=dynamicTimeWarping(centroidCurves->curveArray[0],crvs->curveArray[j]);
				jmin=0;
				minDistance=distance;
				for(j=1; j<numOfCentroids; j++)
				{
					if(distType==1)
						distance=discreteFrechetDistance(centroidCurves->curveArray[i],crvs->curveArray[j]);
					else
						distance=dynamicTimeWarping(centroidCurves->curveArray[i],crvs->curveArray[j]);
					if(minDistance>distance)
					{
						jmin=j;
						minDistance=distance;
					}
				}
			}
		}
		free(check);
		return dist;
	}
}

cluster* createCluster(curves* crvs,int *centroids, minDist* centroidsArray, int totalNumOfCurves, int numOfCentroids)
{
	int i, j;
	int* numArray;
	cluster* cls;
	if(centroids!=NULL)
	{
		numArray=malloc(numOfCentroids*sizeof(int));
		memset(numArray,0,numOfCentroids*sizeof(int));
		for( i=0; i<totalNumOfCurves; i++)
		{
			for(j=0; j<numOfCentroids; j++)
			{
				if(centroids[j]==centroidsArray[i].centroid)
				{
					numArray[j]++;
					break;
				}
			}
		}
		cls=(cluster*)malloc(numOfCentroids*sizeof(cluster));
		for(i=0; i<numOfCentroids; i++)
		{
			cls[i].id=i;
			cls[i].centroidId=centroids[i];
			cls[i].numOfCurves=numArray[i];
			cls[i].crvs=malloc(numArray[i]*sizeof(int));
		}
		memset(numArray,0,numOfCentroids*sizeof(int));
		for(i=0; i<totalNumOfCurves; i++)
		{
			for(j=0; j<numOfCentroids; j++)
			{
				if(centroids[j]==centroidsArray[i].centroid)
				{
					cls[j].crvs[numArray[j]]=i;
					numArray[j]++;
				}
			}
		}
		free(numArray);
		return cls;
	}
	else
	{
		numArray=malloc(numOfCentroids*sizeof(int));
		memset(numArray,0,numOfCentroids*sizeof(int));
		for( i=0; i<totalNumOfCurves; i++)
		{
			for(j=0; j<numOfCentroids; j++)
			{
				if(j==centroidsArray[i].centroid)
				{
					numArray[j]++;
					break;
				}
			}
		}
		cls=(cluster*)malloc(numOfCentroids*sizeof(cluster));
		for(i=0; i<numOfCentroids; i++)
		{
			cls[i].id=i;
			cls[i].centroidId=i;
			cls[i].numOfCurves=numArray[i];
			cls[i].crvs=malloc(numArray[i]*sizeof(int));
		}
		memset(numArray,0,numOfCentroids*sizeof(int));
		for(i=0; i<totalNumOfCurves; i++)
		{
			for(j=0; j<numOfCentroids; j++)
			{
				if(j==centroidsArray[i].centroid)
				{
					cls[j].crvs[numArray[j]]=i;
					numArray[j]++;
				}
			}
		}
		free(numArray);
		return cls;
	}
}

void deleteCluster(cluster* cls, int numOfCentroids)
{
	int i;
	for(i=0; i<numOfCentroids; i++)
	{
		free(cls[i].crvs);
	}
	free(cls);
}
