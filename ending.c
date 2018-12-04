#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "input.h"
#include "initial_curves.h"
#include "grid.h"
#include "hash.h"
#include "range.h"
#include "traversal_computation.h"
#include "initialization.h"
#include "assignment.h"
#include "mean_frechet_tree.h"
#include "pam.h"
#include "ending.h"

int meanFrechetCondition(cluster *currentCls,cluster *prevCls,int numOfCulsters)
{
	int i, j;
	for(i=0; i<numOfCulsters; i++)
	{
		if(currentCls[i].numOfCurves!=prevCls[i].numOfCurves)
			return 1;
	}
	for(i=0; i<numOfCulsters; i++)
	{
		for(j=0; j<currentCls[i].numOfCurves; j++)
		{
			if(currentCls[i].crvs[j]!=prevCls[i].crvs[j])
				return 1;
		}
	}
	return 0;
}

int PAM_Condition(cluster *currentCls,cluster *prevCls,int numOfCulsters)
{
	int i;
	for(i=0; i<numOfCulsters; i++)
	{
		if(currentCls[i].centroidId!=prevCls[i].centroidId)
			return 1;
	}
	return 0;
}