#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "initial_curves.h"
#include "traversal_computation.h"


traversalArrayInfo* DFD_traversal( curveInfo queryCurve, curveInfo curve)
{
	int numOfPoints_curve, numOfPoints_query, i, j, k, m, dimension, traversalArraySize, nextAvailablePos, minIdx, nopQ;
	double sum, euclideanDistance, max, min, dist_DFD, minDist;
	double *p, *q;
	double **c;
	indx minIndex;
	indx *traversalArray;
	traversalArrayInfo *info;

	numOfPoints_curve = curve.numOfPoints;
	p = curve.points;

	numOfPoints_query = queryCurve.numOfPoints;
	nopQ = numOfPoints_query;
	q = queryCurve.points;

	dimension = curve.dimension;

	/*** KATASKEUI 2D ARRAY c ***/
	c = malloc( (numOfPoints_query+1) * sizeof(double *));	

	for( i=0; i<numOfPoints_query+1; i++ ){
		c[i] = malloc( (numOfPoints_curve+1) * sizeof(double));
	}

	/*** EFARMOGI ALGORITHMOU DFD ***/
	for( i=1; i<numOfPoints_query+1; i++ )
	{
		for( j=1; j<numOfPoints_curve+1; j++ )
		{
			/*** IPOLOGISMOS EUCLIDEAN DISTANCE ***/
			sum = 0;
			for( k=0; k<dimension; k++ ){
				sum = sum + pow( q[ ((i-1)*dimension) + k] - p[ ((j-1)*dimension) + k], 2 );	// pow(x,y) apo math.h
			}
			euclideanDistance = sqrt(sum);

			/*** ANATHESI TIMIS SE THESI TOU PINAKA c ***/
			if( i == 1 && j == 1 )
			{
				c[i][j] = euclideanDistance;
			}
			else if( i == 1 )
			{

				if( euclideanDistance > c[1][j-1] )
					c[i][j] = euclideanDistance;
				else
					c[i][j] = c[1][j-1];
			}
			else if( j == 1 )
			{
				if( euclideanDistance > c[i-1][1] )
					c[i][j] = euclideanDistance;
				else
					c[i][j] = c[i-1][1];
			}
			else
			{
				max = euclideanDistance;
				min = INFINITY;

				if( c[i-1][j] < min )
					min = c[i-1][j];

				if( c[i][j-1] < min )
					min = c[i][j-1];

				if( c[i-1][j-1] < min )
					min = c[i-1][j-1];

				if( euclideanDistance > min )
					max = euclideanDistance;
				else
					max = min;

				c[i][j] = max;
			}
		}
	}



	/******************************************/
	/*** EDO KSEKINAEI TO BACK-PROPAGATION ***/
	/****************************************/


	/*** DIMIOURGIA STRUCT TYPOU traversalArrayInfo POY PERIEXEI TON TRAVERSAL ARRAY + PLIROFORIES GI' AUTON ***/
	info = malloc( sizeof(traversalArrayInfo) );

	/*** DIMIOURGIA TOY ARRAY POY THA PERIEXEI TO TRAVERSAL ***/
	traversalArraySize = numOfPoints_query + numOfPoints_curve - 1;
	traversalArray = malloc( traversalArraySize * sizeof(indx) ); 
	nextAvailablePos = 0;

	/*** STIN PROTI THESI VAZOUME APEFTHIAS TO TELEUTAIO STOIXEIO TOU PINAKA c ***/
	traversalArray[nextAvailablePos].i = numOfPoints_query;
	traversalArray[nextAvailablePos].j = numOfPoints_curve;
	nextAvailablePos++;

	/*** EURESI TRAVERSAL ***/
	while(numOfPoints_query != 1 && numOfPoints_curve != 1)
	{
		minDist = INFINITY;
		if( c[numOfPoints_query-1][numOfPoints_curve] < minDist )
		{
			minDist = c[numOfPoints_query-1][numOfPoints_curve];
			minIdx = 0;
		}

		if( c[numOfPoints_query][numOfPoints_curve-1] < minDist )
		{
			minDist = c[numOfPoints_query][numOfPoints_curve-1];
			minIdx = 1;
		}

		if( c[numOfPoints_query-1][numOfPoints_curve-1] < minDist )
		{
			minDist = c[numOfPoints_query-1][numOfPoints_curve-1];
			minIdx = 2;
		}

		if( minIdx == 0 )
		{
			numOfPoints_query--;
			traversalArray[nextAvailablePos].i = numOfPoints_query;
			traversalArray[nextAvailablePos].j = numOfPoints_curve;
		}
		else if( minIdx == 1 )
		{
			numOfPoints_curve--;
			traversalArray[nextAvailablePos].i = numOfPoints_query;
			traversalArray[nextAvailablePos].j = numOfPoints_curve;
		}
		else if( minIdx == 2 )
		{
			numOfPoints_query--;
			numOfPoints_curve--;
			traversalArray[nextAvailablePos].i = numOfPoints_query;
			traversalArray[nextAvailablePos].j = numOfPoints_curve;
		}
		else
		{
			printf("Error: minIdx value not 0,1,2\n");
			exit(EXIT_FAILURE);
		}

		nextAvailablePos++;
	}

	
	/*** CORNER CASES ***/
	if( numOfPoints_query != 1 )
	{
		while( numOfPoints_query != 1 )
		{
			numOfPoints_query--;
			traversalArray[nextAvailablePos].i = numOfPoints_query;
			traversalArray[nextAvailablePos].j = numOfPoints_curve;
			nextAvailablePos++;
		}
	}
	else if( numOfPoints_curve != 1 )
	{
		while( numOfPoints_curve != 1 )
		{
			numOfPoints_curve--;
			traversalArray[nextAvailablePos].i = numOfPoints_query;
			traversalArray[nextAvailablePos].j = numOfPoints_curve;
			nextAvailablePos++;
		}
	}


	/*** ANATHESI TIMON STA MELI TOU STRUCT info ***/
	info -> traversalArray = traversalArray; 
	info -> traversalArraySize = traversalArraySize;
	info -> nextAvailablePos = nextAvailablePos;


	/*** KATASTROFI AXREIASTON DOMON (PINAKAS c) ***/
	for( m=0; m<nopQ+1; m++ ){
		free(c[m]);
	}
	free(c);

	return info;
}

curveInfo *meanCurve(curveInfo *c1, curveInfo *c2)
{
	traversalArrayInfo *info;
	curveInfo *mean;
	int i,j,pos1,pos2,same,delay;
	double p1,p2;
	double *temp_points;
	delay=0;
	info=DFD_traversal(*c1,*c2);
	mean=malloc(sizeof(curveInfo));
	mean->numOfPoints=info->nextAvailablePos;
	mean->points=malloc(c1->dimension*mean->numOfPoints*sizeof(double));
	mean->dimension=c1->dimension;
	for(i=0; i<mean->numOfPoints; i++)
	{
		for(j=0; j<mean->dimension; j++)
		{
			pos1=info->traversalArray[info->nextAvailablePos-i-1].i-1;
			pos2=info->traversalArray[info->nextAvailablePos-i-1].j-1;
			p1=c1->points[(pos1)*c1->dimension+j];
			p2=c2->points[(pos2)*c2->dimension+j];
			mean->points[(i-delay)*c1->dimension+j]=(p1+p2)/2;
		}
		if(i>0)
		{
			same=1;
			for(j=0; j<mean->dimension; j++)
			{
				if(mean->points[(i-delay)*c1->dimension+j]!=mean->points[(i-delay-1)*c1->dimension+j])
				{
					same=0;
					break;
				}
			}
			if(same==1)
			{
				delay++;
			}
			if(i%3==0 && i-delay>20)
			{
				delay++;																																																																																																																																																																																																																																																																																																																																																																																																						
			}
		}
	}
	mean->numOfPoints=mean->numOfPoints-delay;
	temp_points=realloc(mean->points,mean->dimension*mean->numOfPoints*sizeof(double));
	if(temp_points!=NULL)
	{
		mean->points=temp_points;
	}
	free(info->traversalArray);
	free(info);
	return mean;
}