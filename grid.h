#ifndef GRID__H__
#define GRID__H__

#include "initial_curves.h"
#define RAND_CONST 1000

int KAPPA;

typedef struct gridInfo
{
	int numOfPoints;
	double** grid_instances;
}gridInfo;

typedef struct grids
{
        int numOfGrids; //Number of grids is the number of curves
        gridInfo* gridArray;
}grids;


double change_coordinate_to_grid(double , double , double );

int array_equal(double* , double* , int );
double* create_grid_instance(double , double** , double* , int , int , int );
void delete_grid_instance(double* );

double*** create_random_ts(int , int );
void delete_ts(double*** , int );

gridInfo* create_gridInfo(double , double*** , curves* , int , int );
void delete_gridInfo(gridInfo* , int , int );

grids* create_grids(double , double*** , curves* , int , int );
void delete_grids(grids* , int );

#endif
