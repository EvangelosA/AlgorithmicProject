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
//#include "curve_similarity.h"
#include "traversal_computation.h"
#include "initialization.h"
#include "assignment.h"
#include "mean_frechet_tree.h"
#include "pam.h"
#include "ending.h"
#include "silhouette.h"
#include "output_file.h"

#define KVEC 3
#define M 4294967291

extern int KAPPA;

int main( int argc, char *argv[] )
{
	int i, j, length, hashType, table_size, update_alg, check_end, iter, numOfPoints_mean;
	int *centroids=NULL;;
	double *silhouetteResults;
	double **distanceStorage=NULL, **distanceStorage_a=NULL, **distanceStorage_b=NULL;
	double ***t;
	FILE *outputFile;
	inputStruct *input;
	inputStruct_user *input_user; 
	curves *crvs;
	curves *centroidCurves=NULL;
	config *conf;
	minDist *centroidsArray;
	cluster *cls;
	cluster *prevCls;
	grids  *grds;
	grids *centroidGrids;
	ls *hashStructures;
	curveInfo *curveArray_mean;
	clock_t start, end;

	cls=NULL;
	prevCls=NULL;
	centroidCurves=NULL;
	centroidGrids=NULL;
	centroidsArray=NULL;
	check_end=1;
	srand( time(NULL) );

	/*** PARSE TO INPUT APO TO COMMAND LINE KAI APOTHIKEUSI ***/
	input = commandLineParser(argc, argv);

	/*** PARSE TO INPUT TOU XRISTI KAI APOTHIKEUSI ***/
	input_user = userParser();

	/*** CREATE TOU ARXEIOU EKSODOU ***/
	outputFile = fopen(input -> outputFilePath, "a+");

	/*** PARSE TO ARXEIO EISODOU KAI APOTHIKEUSI KAMPILON ***/
	crvs = parseAndSave_input( input -> inputFilePath );

	/*** PARSE TO CONFIG FILE KAI APOTHIKEUSI ***/
	conf = parseAndSave_config( input -> configurationFilePath );



	/*** CREATING HASH TABLES ***/
	printf("CREATING HASH TABLES\n");

	KAPPA = conf -> num_of_grid_curves;

	crvs -> maxNumOfPoints_forGrids = crvs -> maxNumOfPoints;

	t = create_random_ts(conf->num_of_hash_tables, crvs->dimension);

	grds = create_grids(0.005, t, crvs, crvs->dimension, conf->num_of_hash_tables);

	if(strcmp(input -> metric, "Frechet") == 0)
	{
		hashType=2;
		update_alg=1;
	}
	else
	{
		hashType=2;
		update_alg=2;
	}
	
	table_size = crvs -> totalNumOfCurves / 8;

	hashStructures = create_ls_array(hashType, time(NULL), KVEC, crvs->dimension, conf->num_of_hash_tables, crvs->maxNumOfPoints_forGrids, table_size);
	hashStructures = create_all_hash_tables(*grds, hashStructures, conf->num_of_hash_tables, KVEC, crvs->maxNumOfPoints_forGrids, crvs->dimension, M, table_size);


	/*** INITIALIZATION ***/
	printf("INITIALIZATION\n");
	if( input_user -> initializationMethod == 1 )
		centroids = randomSelection(crvs, crvs->totalNumOfCurves, conf->num_of_clusters);	// pinakas k (osa kai ta clusters) akeraion, kathe akeraios = thesi curve sto curveArray
	else
		centroids = kmeansplusplus(crvs, crvs->totalNumOfCurves, conf->num_of_clusters, input->metric);

/*
	printf("PRINTING CENTROIDS\n");
	for(i=0; i<conf->num_of_clusters; i++)
		printf("%d\n", centroids[i]);
*/
	iter=0;
	start = clock();
	do
	{
		if(prevCls)
		{
			deleteCluster(prevCls, conf -> num_of_clusters);
		}
		prevCls=cls;

		if(centroidsArray)
		{
			free(centroidsArray);
			centroidsArray=NULL;
		}

		/*** ASSIGNMENT ***/
		printf("ASSIGNMENT\n");
		if( input_user -> assignmentMethod == 1 )
			centroidsArray=loydsAssignment(crvs, crvs->totalNumOfCurves, centroids, conf->num_of_clusters, input->metric, centroidCurves );
		else
			centroidsArray = rangeAssignment(crvs, crvs->totalNumOfCurves, grds, centroids, conf->num_of_clusters, input->metric, hashStructures, conf->num_of_hash_tables, 0.05, KVEC, M, table_size, centroidCurves, centroidGrids);

		cls = createCluster(crvs, centroids, centroidsArray, crvs->totalNumOfCurves, conf->num_of_clusters);
		
		if(iter>0)
		{
			// CONDITION FOR ENDING
			if(update_alg==1) // Frechet mean case
				check_end = meanFrechetCondition( cls, prevCls, conf->num_of_clusters);
			else	// PAM case
				check_end = PAM_Condition( cls, prevCls, conf->num_of_clusters);

			if( check_end == 0 )
				break;
		}

		if(centroids)
		{
			free(centroids);
			centroids=NULL;
		}

		/*** UPDATE ***/
		if(update_alg==1)
		{
			if(centroidCurves)
			{
				destroyInitialStructures(centroidCurves);
				centroidCurves=NULL;
			}
			if(centroidGrids)
			{
				delete_grids(centroidGrids,conf->num_of_hash_tables);
				centroidGrids=NULL;
			}

			printf("MEAN FRECHET UPDATE\n");
			centroidCurves=meanFrechetCentroids(crvs,cls,conf->num_of_clusters);
			//printf("CREATING GRIDS FOR CENTROIDS\n");
			centroidCurves->maxNumOfPoints_forGrids=crvs->maxNumOfPoints;
			crvs->maxNumOfPoints_forGrids=crvs->maxNumOfPoints;
			centroidGrids = create_grids(0.005, t, centroidCurves, centroidCurves->dimension, conf->num_of_hash_tables);
		}
		else
		{
			printf("PARTITION AROUND MEDOIDS UPDATE\n");
			if( distanceStorage )
			{
				delete2D(distanceStorage, crvs->totalNumOfCurves);
				distanceStorage = NULL;
			}
			distanceStorage = createAndInitialize2D(crvs->totalNumOfCurves);
			centroids = PartitionAroundMedoids(cls, conf->num_of_clusters, centroidsArray, crvs->totalNumOfCurves, crvs->curveArray, distanceStorage);
/*
			int cc;
			for(cc=0; cc<conf->num_of_clusters; cc++)
			{
				printf("%d\n", centroids[cc]);
			}
*/
		}

		iter++;
	}while( check_end );
	end = clock();

	double ellapsedTime = ((double)(end - start)) / CLOCKS_PER_SEC;

	// HERE CODE FOR SILHOUET
	printf("SILHOUETTE\n");
	distanceStorage_a = createAndInitialize2D(crvs->totalNumOfCurves);	//gia ta kalitera clusters
	distanceStorage_b = createAndInitialize2D(crvs->totalNumOfCurves);	//gia ta deutera kalitera clusters
	silhouetteResults = silhouetteCalculation(cls, conf->num_of_clusters, centroidsArray, crvs->curveArray, centroidCurves, crvs->totalNumOfCurves, update_alg, distanceStorage_a, distanceStorage_b);

/*
	for(i=0; i<conf->num_of_clusters; i++)
		printf("Cluster%d: %f\n", i, silhouetteResults[i]);

	printf("overall: %f\n", silhouetteResults[i]);
*/

	/*** OUTPUT ***/
	writeToFile(outputFile, input->metric, input_user->initializationMethod, input_user->assignmentMethod, conf->num_of_clusters, update_alg, cls, silhouetteResults, centroidCurves, ellapsedTime);

	/*** FREE ALLOCATED MEMMORY ***/
	if( update_alg == 2)
		delete2D(distanceStorage, crvs->totalNumOfCurves);
	
	free(silhouetteResults);
	delete2D(distanceStorage_a, crvs->totalNumOfCurves);
	delete2D(distanceStorage_b, crvs->totalNumOfCurves);

	delete_ls_array(hashStructures, conf->num_of_hash_tables, KVEC, table_size);
	delete_grids(grds, conf->num_of_hash_tables);
	if(update_alg==1)
	{
		delete_grids(centroidGrids, conf->num_of_hash_tables);
	}
	delete_ts(t, conf->num_of_hash_tables);
	if(update_alg==1)
	{
		destroyInitialStructures(centroidCurves);
	}
	deleteCluster(cls,conf->num_of_clusters);
	deleteCluster(prevCls,conf->num_of_clusters);
	free(centroidsArray);

	free(input -> inputFilePath);
	free(input -> outputFilePath);
	free(input -> configurationFilePath);
	free(input -> metric);
	free(input);
	free(input_user);
	free(conf);
	if(centroids)
	{
		free(centroids);
		centroids=NULL;
	}

	destroyInitialStructures(crvs);
	fclose(outputFile);

	return 0;	
}