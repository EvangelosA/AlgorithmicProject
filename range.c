#include <stdio.h>
#include <stdlib.h>
#include "grid.h"
#include "initial_curves.h"
#include "curve_similarity.h"
#include "hash.h"
#include "range.h"
#define MIN_VALUE 10000000

int same_instances(gridInfo p, gridInfo q, int l, int max_points, int dim) // Check if two instances are the same
{
	if(p.numOfPoints != q.numOfPoints)
	{
		return 0;
	}
	else
	{
		int i;
		for(i=0; i<=max_points*KAPPA*dim-1; i++)
		{
			if(p.grid_instances[l][i] != q.grid_instances[l][i])
			{
				return 0;
			}
		}
		return 1;
	}
}

int* bucket_range(bucket b, curves crvs, grids g_c, curveInfo q, gridInfo q_i, int l, double r, double r2, int dist_type, int* range_table, int dim) //Finding the R-nearest in a bucket
{
	int exist = 0;
	int i;
	double dist;

	for(i=0; i<=b.num_of_curves-1; i++)
	{
		if(same_instances(g_c.gridArray[b.curves[i]], q_i, l, crvs.maxNumOfPoints_forGrids, dim))
		{
			exist = 1;
			if(dist_type == 1)
			{
				dist=discreteFrechetDistance(q, crvs.curveArray[b.curves[i]]);
				if(dist <= r && dist > r2)
				{
					range_table[b.curves[i]] = 1;
				}
			}
			else if(dist_type == 2)
			{
				dist=dynamicTimeWarping(q, crvs.curveArray[b.curves[i]]);
				if(dist <= r && dist > r2)
				{
					range_table[b.curves[i]] = 1;
				}
			}
		}
	}
	if(exist == 0)
	{
		for(i=0; i<=b.num_of_curves-1; i++)
		{
                        if(dist_type == 1)
                        {
                        	dist=discreteFrechetDistance(q, crvs.curveArray[b.curves[i]]);
                                if(dist <= r && dist > r2)
                                {
                                        range_table[b.curves[i]] = 1;
                                }
                        }
                        else if(dist_type == 2)
                        {
                        	dist=dynamicTimeWarping(q, crvs.curveArray[b.curves[i]]);
                                if(dist <= r && dist >r2)
                                {
                                        range_table[b.curves[i]] = 1;
                                }
                        }
		}
	}
	return range_table;
}

int* hash_structure_range(curves *crvs, grids *crvs_grds, curveInfo q, gridInfo q_g, ls* structures, int structures_num, double r, double r2, int distType,int k, long long int m, int table_size)//Finding all the nearest of a Q
{
	int *range_array;
	int i, j, value;
	range_array = (int*)malloc(crvs->totalNumOfCurves * sizeof(int));
	for(i=0; i<=crvs->totalNumOfCurves-1; i++)
	{
		range_array[i] = 0;
	}
	for(j=0; j<=structures_num-1; j++)
	{
		value = general_hash_function(q_g.grid_instances[j], structures[j], k, crvs->maxNumOfPoints_forGrids, crvs->dimension, m, table_size);
		range_array=bucket_range(structures[j].table[value], *crvs, *crvs_grds, q, q_g, j, r, r2, distType, range_array, crvs->dimension);
	}
	return range_array;
}

int** total_range(curves* crvs, grids* grds, curves* qrs, grids* query_grds, ls* structures, int l, double r, int distType, int k, long long int m, int table_size)//Creating all the results of all queries
{
	int** total;
	int *check;
	total = (int**)malloc(qrs->totalNumOfCurves * sizeof(int*));
	int i;
	for(i=0; i<=qrs->totalNumOfCurves-1; i++)
	{
		check = hash_structure_range(crvs, grds, qrs->curveArray[i], query_grds->gridArray[i], structures, l, r, 0.0, distType, k, m, table_size);
		total[i]=flatten_order(check, crvs);
		free(check);
	}
	return total;
}

void delete_range_results(int** array, int length)//Deleting results
{
	int i;
	for(i=0; i<=length-1; i++)
	{
		delete_check_array(array[i]);
	}
	free(array);
}

int* flatten_order(int* array, curves* crvs)//Flat and order function for check array
{
	int* ordered_array;
	int i, count;
	count = 0;
	for(i=0; i<=crvs->totalNumOfCurves-1; i++)
	{
		if(array[i] == 1)
		{
			count++;
		}
	}
	int j, min, jmin;
	ordered_array = (int*)malloc((count+1) * sizeof(int));
	for(i=0; i<=count-1; i++)
	{
		min = MIN_VALUE;
		jmin = -1;
		for(j=0; j<=crvs->totalNumOfCurves-1; j++)
		{
			if(crvs->curveArray[j].curveID < min && array[j] == 1)
			{
				min = crvs->curveArray[j].curveID;
				jmin = j;
			}
		}
		array[jmin] = 0;
		ordered_array[i] = min;
	}
	ordered_array[i] = -1;
	return ordered_array;
}

int* range_all_function(curves crvs, curveInfo q, double r, int dist_type)//All range function - Not in use
{
	int i;
	int* check_curves;
	check_curves = (int*)malloc(crvs.totalNumOfCurves * sizeof(int));
	if(dist_type == 1)
	{
		for(i=0; i<=crvs.totalNumOfCurves-1; i++)
		{
			if(discreteFrechetDistance(q, crvs.curveArray[i]) <= r)
			{
				check_curves[i] = 1;
			}
		}
	}
	else if(dist_type == 2)
	{
		for(i=0; i<=crvs.totalNumOfCurves-1; i++)
		{
			if(dynamicTimeWarping(q, crvs.curveArray[i]) <= r)
			{
				check_curves[i] = 1;
			}
		}
	}
	return check_curves;
}

void print_range(curves crvs, int* check)//Printing the R-nearest of a centroid
{
	int i;
	i=0;
	while(check[i] != -1)
	{
		fprintf(stdout, "%d:%d\n", i, check[i]);
		i++;
	}
}

void delete_check_array(int* array)//Deleting the check array of a centroid
{
	free(array);
}
