#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "initial_curves.h"
#include "grid.h"

double change_coordinate_to_grid(double delta, double t, double point_coord) //Here i calculate for a specific coordinate, delta and t, what is the nearest new coordinate in the grid
{
	double low_fragment, high_fragment;
	low_fragment = ( (int) ((point_coord-t)/delta) ) * delta + t;//the type that calculates the low fragment
	high_fragment = low_fragment + delta;
	if( fabs(low_fragment-point_coord) < fabs(high_fragment-point_coord) )
	{
		return low_fragment;
	}
	else
	{
		return high_fragment;
	}
}

int array_equal(double* p1, double* p2, int length)
{
	int i;
	for( i=0; i<=length-1; i++ )
	{
		if( p1[i]!=p2[i] )
		{
			return 0;
		}
	}
	return 1;
}

double* create_grid_instance(double delta, double** t, double* curve, int num_of_points, int dimensions, int max_points)//we create the grid instance from a single curve
{
	double* new_instance;
	int i, j, k, w;
	double coord;
	new_instance = malloc(max_points * KAPPA * dimensions * sizeof(double));
	for( i=0; i<=KAPPA-1; i++ )
	{
		k=0;
		if( num_of_points == 1 )
		{
			for( j=0; j<=dimensions-1; j++ ){
				new_instance[i * max_points * dimensions + j] = change_coordinate_to_grid(delta, t[i][j], curve[j]);
			}
		}
		else if( num_of_points > 1 )
		{
			for( j=0; j<=dimensions*2-1; j++)
			{
				w = j % dimensions;
				coord = curve[j];
				new_instance[i*max_points*dimensions+j] = change_coordinate_to_grid(delta, t[i][w], coord);
			}
			for( j=dimensions*2; j<=num_of_points*dimensions-1; j++ )
			{
				if( j % dimensions == 0 )
				{
					double* p1;
					double* p2;
					p1 = &(new_instance[i*max_points*dimensions + j - dimensions*2 - k]);
					p2 = &(new_instance[i*max_points*dimensions + j - dimensions - k]);
					if( array_equal(p2, p1, dimensions) ){
						k = k + dimensions;
					}
				}
				new_instance[i*max_points*dimensions+j-k] = change_coordinate_to_grid(delta, t[i][j%dimensions], curve[j]);
			}
			if( num_of_points>1 )
			{
				double* p1;
				double* p2;
				p1 = &(new_instance[i*max_points*dimensions+j-dimensions*2-k]);
				p2 = &(new_instance[i*max_points*dimensions+j-dimensions-k]);
				if( array_equal(p1, p2, dimensions) ){
					k = k + dimensions;
				}
			}
		}
		for( j=num_of_points*dimensions; j<=max_points*dimensions+k-1; j++ ){
			new_instance[i*max_points*dimensions + j - k] = 0.0;
		}
	}
	return new_instance;
}

void delete_grid_instance(double* instance) //Delete the grid instance
{
	free(instance);
}

double*** create_random_ts(int l, int dimensions) //this functions creates a 3D array for t's
{
	double*** new_array;
	int i, j, k;
	new_array = (double***)malloc(l * sizeof(double**));
	for(i=0;i<=l-1;i++)
	{
		new_array[i] = (double**)malloc(KAPPA * sizeof(double*));
		for(j=0; j<=KAPPA-1; j++)
		{
			new_array[i][j] = (double*)malloc(dimensions * sizeof(double));
			for(k=0; k<=dimensions-1; k++)
			{
				new_array[i][j][k] = (((double)(rand()%RAND_CONST))/RAND_CONST) * ((double)dimensions);
			}
		}
	}
	return new_array;
}

void delete_ts(double*** t, int l) //Deleting t's
{
	int i, j, k;
	for(i=0; i<=l-1; i++)
	{
		for(j=0; j<=KAPPA-1; j++)
		{
			free(t[i][j]);
		}
		free(t[i]);
	}
	free(t);
}

gridInfo* create_gridInfo(double delta, double*** t, curves* crvs, int dimensions, int l) //Creates the gridInfo for every curve
{
	gridInfo* new_grid;
	int i, j;
	new_grid = (gridInfo*)malloc(crvs->totalNumOfCurves * sizeof(gridInfo));
	for(i=0; i<=(crvs->totalNumOfCurves)-1; i++)
	{
		(new_grid[i]).grid_instances = (double**)malloc(l * sizeof(double*));
		for(j=0; j<=l-1; j++)
		{
			(new_grid[i]).grid_instances[j] = create_grid_instance(delta, t[j], (crvs->curveArray[i]).points, (crvs->curveArray[i]).numOfPoints, dimensions, crvs->maxNumOfPoints_forGrids);
		}
		(new_grid[i]).numOfPoints = (crvs->curveArray[i]).numOfPoints;
	}
	return new_grid;
}

void delete_gridInfo(gridInfo* grd, int totalNumOfCurves, int l) //Deletes all the gridInfo array
{
	int i, j;
	for( i=0; i<=totalNumOfCurves-1; i++ )
	{
		for( j=0; j<=l-1; j++ )
		{
			delete_grid_instance((grd[i]).grid_instances[j]);
		}
		free((grd[i]).grid_instances);
	}
	free(grd);
}

grids* create_grids(double delta, double*** t, curves* crvs, int dimensions, int l) //this creates the grids struct (use only this in main)
{
	grids* new_grid;
	new_grid = (grids*)malloc(sizeof(grids));
	new_grid -> numOfGrids = crvs -> totalNumOfCurves;
	new_grid -> gridArray = create_gridInfo(delta, t, crvs, dimensions, l);
	return new_grid;
}

void delete_grids(grids* grd, int l) //Deleting grids struct
{
	delete_gridInfo(grd->gridArray, grd->numOfGrids, l);
	free(grd);
}
