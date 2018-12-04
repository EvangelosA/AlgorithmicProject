#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "initial_curves.h"

curves* parseAndSave_input(char *path)
{
	int curveID, i, size, totalCurves, nextAvailablePos, nextAvailablePos_points, maxNumOfPoints, numOfPoints, count, dimension, dimFlag;
	double *points;
	char *split, datasetLine[lineSize_dataset], *d="@dimension";
	FILE *dataset, *dataset_firstRun;
	curves *crvs;
	curveInfo *curveArray;

	totalCurves = 0;
	nextAvailablePos = 0;
	size = 1;

	crvs = malloc( sizeof(curves) );							//dimiourgia struct tipou curves pou krataei plirofories gia ton pinaka curveArray pou periexei oles tis kampiles
	crvs -> curveArray_size = size;
	crvs -> nextAvailablePos = 0;
	crvs -> totalNumOfCurves = 0;
	crvs -> maxNumOfPoints_forGrids = 0;
	crvs -> curveArray = malloc( size * sizeof(curveInfo) );	//dimiourgia pinaka apo curves pou krataei plirofories gia kathe curve tou dataset
	curveArray = crvs -> curveArray;

	dataset = fopen(path, "r");
	dataset_firstRun = dataset;	//antigrafo

	maxNumOfPoints = 0;
	count = 0;
	dimFlag = 0;
	while( fgets(datasetLine, lineSize_dataset, dataset_firstRun) )		//kanoume ena run sthn arxh gia na vroume ton max arithmo simeion pou mporei na exei mia curve kaina paroume to dimension
	{
		if( count == 0 )
		{
			split = strtok(datasetLine, " \t\r\n");
			if( strcmp(d, split) == 0 )
			{
				split = strtok(NULL, " \t\r\n");
				dimension = atoi(split);
				dimFlag = 1;					// simainei oti exoume 1h grammi pou deixnei to dimension
			}
			else
			{
				dimension = 2;
				split = strtok(NULL, " \t\r\n");
				numOfPoints = atoi(split);
				if( numOfPoints > maxNumOfPoints )
					maxNumOfPoints = numOfPoints;
			}
			crvs -> dimension = dimension;
		}
		else
		{
			split = strtok(datasetLine, " \t\r\n");
			split = strtok(NULL, " \t\r\n");
			//printf("%s\n", split);
			numOfPoints = atoi(split);
			if( numOfPoints > maxNumOfPoints )
				maxNumOfPoints = numOfPoints;
		}
		count++;
	}

	//printf("max number of points: %d\n", maxNumOfPoints);

	fclose(dataset);
	dataset = fopen(path, "r");	//reopen

	memset(datasetLine, 0, lineSize_dataset);

	count = 0;
	while( fgets(datasetLine, lineSize_dataset, dataset) )	//deutero perasma gia na kataskeuasoume oses extra domes xreiazetai kai na tis gemisoume
	{
		if( count == 0 && dimFlag == 1)
		{
			memset(datasetLine, 0, lineSize_dataset);	// discard the first line
		}
		else
		{
			split = strtok(datasetLine, " ,()\t\r\n");
			if(split)
			{
				//printf("curve: %d\n", totalCurves);

				if( nextAvailablePos == size )	//tote realloc giati gemise o pinakas me tis curves
				{
					crvs -> curveArray = realloc(crvs -> curveArray, (2 * size) * sizeof(curveInfo) );
					curveArray = crvs -> curveArray;
					size = 2 * size;
				}

			 
				curveArray[nextAvailablePos].curveID = atoi(split);
				curveArray[nextAvailablePos].dimension = dimension;

				split = strtok(NULL, " ,()\t\r\n");
				curveArray[nextAvailablePos].numOfPoints = atoi(split);
				numOfPoints = curveArray[nextAvailablePos].numOfPoints;	//antigrafo
				//printf("numOfPoints: %d\n", numOfPoints);

				curveArray[nextAvailablePos].points = malloc( (curveArray[nextAvailablePos].numOfPoints * dimension) * sizeof(double) );
				points = curveArray[nextAvailablePos].points;	//antigrafo
				
				nextAvailablePos_points = 0;

				for( i=0; i<numOfPoints; i++ )	//gemisma pinaka points me simeia
				{
					split = strtok(NULL, " ,()\t\r\n");
					points[nextAvailablePos_points] = atof(split);

					nextAvailablePos_points++;

					split = strtok(NULL, " ,()\t\r\n");
					points[nextAvailablePos_points] = atof(split);

					nextAvailablePos_points++;
					//printf("nextAvailablePos_points: %d\n", nextAvailablePos_points);
				}

				totalCurves++;
				nextAvailablePos++;
			}
			else
				printf("split is NULL\n");
		}
		count++;
	}
	
	crvs -> totalNumOfCurves = totalCurves;
	crvs -> nextAvailablePos = nextAvailablePos;
	crvs -> curveArray_size = size;
	crvs -> maxNumOfPoints = maxNumOfPoints;

	fclose(dataset);
	
	return crvs;
}

config* parseAndSave_config( char *path )
{
	int num_of_clusters, num_of_grid_curves, num_of_hash_tables;
	char buffer[lineSize_configFile];
	char *split, *clusters="number_of_clusters", *grid_curves="number_of_grid_curves", *hash_tables="number_of_hash_tables";
	FILE *configFile;
	config *configValues;

	configValues = malloc( sizeof(config) );
	configValues -> num_of_clusters = -1;
	configValues -> num_of_grid_curves = -1;
	configValues -> num_of_hash_tables = -1;

	configFile = fopen(path, "r");

	while( fgets(buffer, lineSize_configFile, configFile) )
	{
		split = strtok( buffer, ":\r\n" );
		if( strcmp(buffer, clusters) == 0 )
		{
			split = strtok( NULL, ":\r\n" );
			configValues -> num_of_clusters = atoi(split);
		}
		else if( strcmp(buffer, grid_curves) == 0 )	
		{
			split = strtok( NULL, ":\r\n" );
			configValues -> num_of_grid_curves = atoi(split);
		}
		else if( strcmp(buffer, hash_tables) == 0 )	
		{
			split = strtok( NULL, ":\r\n" );
			configValues -> num_of_hash_tables = atoi(split);
		}
		else
		{
			printf("Error: config file has different format than expected\n");
			exit(EXIT_FAILURE);
		}

		memset(buffer, 0, lineSize_configFile);
	}

	//an kapoio den exei parei timi apo to arxeio dose default times
	if( configValues -> num_of_grid_curves == -1 )
		configValues -> num_of_grid_curves = 2;

	if( configValues -> num_of_hash_tables == -1 )
		configValues -> num_of_hash_tables = 3;

	fclose(configFile);
	
	return configValues;

}

void destroyInitialStructures(curves *crvs)
{
	int i;
	
	for( i=0; i<crvs->totalNumOfCurves; i++ )
		free(crvs -> curveArray[i].points);
	
	free(crvs -> curveArray);
	free(crvs);
}