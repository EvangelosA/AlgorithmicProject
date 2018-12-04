#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "initial_curves.h"
#include "grid.h"
#include "hash.h"
#include "range.h"
#include "assignment.h"
#include "traversal_computation.h"
#include "mean_frechet_tree.h"

int pickCurveRandomly(int *checkArray,int numOfCurves)
{
	int num,start;
	num=rand();
	num= num % numOfCurves;
	start=num;
		do{
		if(checkArray[num]==0)
		{
			checkArray[num]=1;
			return num;
		}
		num++;
		num= num % numOfCurves;
	}
	while(num!=start);
	return -1;
}

int *initializeStartLevel(cluster cls)
{
	int i;
	int *startArray;
	int *checkArray;
	checkArray=malloc(sizeof(int)*cls.numOfCurves);
	for(i=0;i<cls.numOfCurves;i++)
	{
		checkArray[i]=0;
	}
	startArray=malloc(sizeof(int)*cls.numOfCurves);
	for(i=0 ; i<cls.numOfCurves ; i++)
	{
		startArray[i]=cls.crvs[pickCurveRandomly(checkArray,cls.numOfCurves)];
	}
	free(checkArray);
	return startArray;
}

curves *initializeFirstLevel(curves *crvs,int * startLevel,int numOfStartCurves)
{
	curves *firstLevel;
	int numFirstCurves, i, j;
	curveInfo *temp;
	if(numOfStartCurves%2==0)
		numFirstCurves=numOfStartCurves/2;
	else
		numFirstCurves=(numOfStartCurves+1)/2;
	firstLevel=malloc(sizeof(curves));
	firstLevel->dimension=crvs->dimension;

	firstLevel->totalNumOfCurves=numFirstCurves;
	firstLevel->curveArray=malloc(numFirstCurves*sizeof(curveInfo));
	if(numOfStartCurves%2==0)
	{
		for(i=0; i<numFirstCurves; i++)
		{
			temp=meanCurve(&crvs->curveArray[startLevel[i*2]],&crvs->curveArray[startLevel[i*2+1]]);
			firstLevel->curveArray[i].points=temp->points;
			firstLevel->curveArray[i].dimension=temp->dimension;
			firstLevel->curveArray[i].numOfPoints=temp->numOfPoints;
			firstLevel->curveArray[i].curveID=i;
			free(temp);
		}
		return firstLevel;
	}
	else
	{
		for(i=0; i<numFirstCurves-1; i++)
		{
			temp=meanCurve(&crvs->curveArray[startLevel[i*2]],&crvs->curveArray[startLevel[i*2+1]]);
			firstLevel->curveArray[i].points=temp->points;
			firstLevel->curveArray[i].dimension=temp->dimension;
			firstLevel->curveArray[i].numOfPoints=temp->numOfPoints;
			firstLevel->curveArray[i].curveID=i;
			free(temp);
		}
		firstLevel->curveArray[numFirstCurves-1].points=malloc(sizeof(double)*crvs->curveArray[startLevel[numOfStartCurves-1]].numOfPoints*crvs->dimension);
		for(j=0; j<crvs->curveArray[startLevel[numOfStartCurves-1]].numOfPoints*crvs->dimension; j++)
			firstLevel->curveArray[numFirstCurves-1].points[j]=crvs->curveArray[startLevel[numOfStartCurves-1]].points[j];
		firstLevel->curveArray[numFirstCurves-1].dimension=crvs->dimension;
		firstLevel->curveArray[numFirstCurves-1].numOfPoints=crvs->curveArray[startLevel[numOfStartCurves-1]].numOfPoints;
		return firstLevel;
	}
}

curves *produceLevel(curves *prevLevel)
{
	curves *firstLevel;
	int numFirstCurves, i, j;
	curveInfo *temp;
	if(prevLevel->totalNumOfCurves%2==0)
		numFirstCurves=prevLevel->totalNumOfCurves/2;
	else
		numFirstCurves=(prevLevel->totalNumOfCurves+1)/2;
	firstLevel=malloc(sizeof(curves));
	firstLevel->dimension=prevLevel->dimension;

	firstLevel->totalNumOfCurves=numFirstCurves;
	firstLevel->curveArray=malloc(numFirstCurves*sizeof(curveInfo));
	if(prevLevel->totalNumOfCurves%2==0)
	{
		for(i=0; i<numFirstCurves; i++)
		{
			temp=meanCurve(&prevLevel->curveArray[i*2],&prevLevel->curveArray[i*2+1]);
			firstLevel->curveArray[i].points=temp->points;
			firstLevel->curveArray[i].dimension=temp->dimension;
			firstLevel->curveArray[i].numOfPoints=temp->numOfPoints;
			firstLevel->curveArray[i].curveID=i;
			free(temp);
		}
		return firstLevel;
	}
	else
	{
		for(i=0; i<numFirstCurves-1; i++)
		{
			temp=meanCurve(&prevLevel->curveArray[i*2],&prevLevel->curveArray[i*2+1]);
			firstLevel->curveArray[i].points=temp->points;
			firstLevel->curveArray[i].dimension=temp->dimension;
			firstLevel->curveArray[i].numOfPoints=temp->numOfPoints;
			firstLevel->curveArray[i].curveID=i;
			free(temp);
		}
		firstLevel->curveArray[numFirstCurves-1].points=malloc(sizeof(double)*prevLevel->curveArray[prevLevel->totalNumOfCurves-1].numOfPoints*prevLevel->dimension);
		for(j=0; j<prevLevel->curveArray[prevLevel->totalNumOfCurves-1].numOfPoints*prevLevel->dimension; j++)
			firstLevel->curveArray[numFirstCurves-1].points[j]=prevLevel->curveArray[prevLevel->totalNumOfCurves-1].points[j];
		firstLevel->curveArray[numFirstCurves-1].dimension=prevLevel->dimension;
		firstLevel->curveArray[numFirstCurves-1].numOfPoints=prevLevel->curveArray[prevLevel->totalNumOfCurves-1].numOfPoints;
		return firstLevel;
	}
}

curveInfo *chooseCentroid(curves *crvs,cluster cls)
{
	curveInfo *centroid;
	int *startLevel;
	int i;
	curves *level1;
	curves *level2;
	if (cls.numOfCurves==0)
	{
		centroid=malloc(sizeof(curveInfo));
		centroid->dimension=crvs->dimension;
		centroid->numOfPoints=0;
		centroid->points=NULL;
		return centroid;
	}
	else if(cls.numOfCurves==1)
	{
		centroid=malloc(sizeof(curveInfo));
		centroid->dimension=crvs->dimension;
		centroid->numOfPoints=crvs->curveArray[cls.crvs[0]].numOfPoints;
		centroid->points=malloc(sizeof(double)*centroid->dimension*centroid->numOfPoints);
		for(i=0; i<centroid->numOfPoints*centroid->dimension; i++)
			centroid->points[i]=crvs->curveArray[cls.crvs[0]].points[i];
		return centroid;
	}
	else
	{
		startLevel=initializeStartLevel(cls);
		level1=initializeFirstLevel(crvs,startLevel,cls.numOfCurves);
		while(level1->totalNumOfCurves>1)
		{
			level2=produceLevel(level1);
			destroyInitialStructures(level1);
			level1=level2;
		}
		centroid=level1->curveArray;
		free(level1);
		free(startLevel);
		return centroid;
	}
}

curves *meanFrechetCentroids(curves *crvs,cluster *cls,int numOfClusters)
{
	int i, maxNumOfPoints;
	curveInfo *temp;
	curves *centroids;
	maxNumOfPoints=0;
	centroids=malloc(sizeof(curves));
	centroids->dimension=crvs->dimension;
	centroids->totalNumOfCurves=numOfClusters;
	centroids->curveArray_size=numOfClusters;
	centroids->curveArray=malloc(numOfClusters*sizeof(curveInfo));
	for(i=0; i<numOfClusters; i++)
	{
		temp=chooseCentroid(crvs,cls[i]);
		centroids->curveArray[i].curveID=i;
		centroids->curveArray[i].dimension=crvs->dimension;
		centroids->curveArray[i].numOfPoints=temp->numOfPoints;
		centroids->curveArray[i].points=temp->points;
		free(temp);
		if(maxNumOfPoints<centroids->curveArray[i].numOfPoints)
		{
			maxNumOfPoints=centroids->curveArray[i].numOfPoints;
		}
	}
	centroids->maxNumOfPoints=maxNumOfPoints;
	return centroids;
}

