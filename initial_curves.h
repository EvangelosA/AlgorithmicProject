#ifndef INITIAL_CURVES__H__
#define INITIAL_CURVES__H__

#define lineSize_dataset 30000
#define lineSize_configFile 50

typedef struct curveInfo
{
	int curveID;
	int dimension;
	int numOfPoints;
	int nextAvailablePos;
	double *points;

}curveInfo;

typedef struct curves
{
	int dimension;
	int curveArray_size;
	int nextAvailablePos;
	int totalNumOfCurves;
	int maxNumOfPoints;
	int maxNumOfPoints_forGrids;
	curveInfo *curveArray;
}curves;

typedef struct config
{
	int num_of_clusters;
	int num_of_grid_curves;
	int num_of_hash_tables;
}config;

curves* parseAndSave_input(char *);
config* parseAndSave_config(char *);
void destroyInitialStructures(curves *);

#endif