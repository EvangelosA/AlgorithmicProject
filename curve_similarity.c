#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "initial_curves.h"
#include "curve_similarity.h"

//opos stis diafaneies: p kai q pinakes me ta dimeia ton 2 kampilon

double discreteFrechetDistance( curveInfo queryCurve, curveInfo curve)
{
	int numOfPoints_curve, numOfPoints_query, i, j, k, m, dimension;
	double sum, euclideanDistance, max, min, dist_DFD;
	double *p, *q;
	double **c;

	numOfPoints_curve = curve.numOfPoints;
	p = curve.points;

	numOfPoints_query = queryCurve.numOfPoints;
	q = queryCurve.points;

	dimension = curve.dimension;

	/* KATASKEUI 2D ARRAY */
	c = malloc( (numOfPoints_query+1) * sizeof(double *));	

	for( i=0; i<numOfPoints_query+1; i++ ){
		c[i] = malloc( (numOfPoints_curve+1) * sizeof(double));
	}
	/*___________________*/


	/* efarmogi algorithmou DFD */
	for( i=1; i<numOfPoints_query+1; i++ )
	{
		for( j=1; j<numOfPoints_curve+1; j++ )
		{
			/* ipologismos euclidean distance */
			sum = 0;
			for( k=0; k<dimension; k++ ){
				sum = sum + pow( q[ ((i-1)*dimension) + k] - p[ ((j-1)*dimension) + k], 2 );	// pow(x,y) apo math.h
			}
			euclideanDistance = sqrt(sum);


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

	dist_DFD = c[numOfPoints_query][numOfPoints_curve];


	/* katastrefoume o,ti domes eftiakse i sinartisi DTW */
	for( m=0; m<numOfPoints_query+1; m++ ){
		free(c[m]);
	}
	free(c);

	return dist_DFD;


}

double dynamicTimeWarping( curveInfo queryCurve, curveInfo curve )
{
	int numOfPoints_curve, numOfPoints_query, i, j, k, m, dimension;
	double sum, euclideanDistance, min, dist_DTW;
	double *p, *q;
	double **c;

	numOfPoints_curve = curve.numOfPoints;
	p = curve.points;

	numOfPoints_query = queryCurve.numOfPoints;
	q = queryCurve.points;

	dimension = curve.dimension;

	/* KATASKEUI 2D ARRAY */
	c = malloc( (numOfPoints_query+1) * sizeof(double *));	

	for( i=0; i<numOfPoints_query+1; i++ ){
		c[i] = malloc( (numOfPoints_curve+1) * sizeof(double));
	}
	/*___________________*/


	/* efarmogi algorithmou DTW */
	for( i=1; i<numOfPoints_query+1; i++){
		c[i][0] = INFINITY;
	}

	for( j=1; j<numOfPoints_curve+1; j++){
		c[0][j] = INFINITY;
	}

	c[0][0] = 0;
	
	for( i=1; i<numOfPoints_query+1; i++ )	// i ~> q (query)
	{
		for( j=1; j<numOfPoints_curve+1; j++ )	// j ~> p (curve)
		{
			/* ipologismos euclidean distance */
			sum = 0;
			for( k=0; k<dimension; k++ ){
				//sum = sum + pow( q[ ((i-1)*dimension) + k] - p[ ((j-1)*dimension) + k] ,2 );	// pow(x,y) apo math.h
				sum = sum + (( q[ ((i-1)*dimension) + k] - p[ ((j-1)*dimension) + k] )*( q[ ((i-1)*dimension) + k] - p[ ((j-1)*dimension) + k] ) );	// pow(x,y) apo math.h
			}
			euclideanDistance = sqrt(sum);
			//euclideanDistance = sum/2;

			/* ipologismos min{ c[i-1][j], c[i][j-1], c[i-1][j-1] } */	
			min = c[i-1][j];

			if( c[i][j-1] < min ){
				min = c[i][j-1];
			}

			if( c[i-1][j-1] < min ){
				min = c[i-1][j-1];
			}

			c[i][j] = euclideanDistance + min;
		}
	}

	dist_DTW = c[numOfPoints_query][numOfPoints_curve];


	/* katastrefoume o,ti domes eftiakse i sinartisi DTW */
	for( m=0; m<numOfPoints_query+1; m++ ){
		free(c[m]);
	}
	free(c);

	return dist_DTW;
}