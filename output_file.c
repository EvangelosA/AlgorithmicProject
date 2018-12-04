#include <stdio.h>
#include "grid.h"
#include "hash.h"
#include "assignment.h"
#include "output_file.h"

void writeToFile(FILE* outputFile, char *metric, int initializationMethod, int assignmentMethod, int num_of_clusters, int update_alg, cluster* cls, double* silhouetteResults, curves* centroidCurves, double ellapsedTime)
{
	int i, j, numOfPoints_mean;
	curveInfo *curveArray_mean;

	fprintf(outputFile, "Algorithm: %dx%dx%dx\n", initializationMethod, assignmentMethod, update_alg);
	fprintf(outputFile, "Metric: %s\n", metric);
	if( update_alg == 1 )
	{
		for( i=0; i<num_of_clusters; i++ )
		{
			fprintf(outputFile, "CLUSTER-%d (size: %d, centroid: [", i, centroidCurves->totalNumOfCurves);
			curveArray_mean = centroidCurves -> curveArray;
			numOfPoints_mean = curveArray_mean[i].numOfPoints;
			for( j=0; j<numOfPoints_mean; j++ )
			{
				if( j != numOfPoints_mean-1 )
					fprintf(outputFile, "%f, ", curveArray_mean[i].points[j]);
				else
					fprintf(outputFile, "%f]\n", curveArray_mean[i].points[j]);
			}
		}
	}
	else
	{
		for( i=0; i<num_of_clusters; i++ )
			fprintf(outputFile, "CLUSTER-%d (size: %d, centroid: %d)\n", i, cls[i].numOfCurves, cls[i].centroidId);
	}

	fprintf(outputFile, "clustering_time: %f seconds\n", ellapsedTime);

	fprintf(outputFile, "Silhouette: [");
	for( i=0; i<num_of_clusters; i++  )
		fprintf(outputFile, "%f, ", silhouetteResults[i]);

	fprintf(outputFile, "%f]\n", silhouetteResults[num_of_clusters]);
	fprintf(outputFile, "\n" );

	return;
}
