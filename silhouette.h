#ifndef SILHOUETTE__H__
#define SILHOUETTE__H__

void secondBestClusters(cluster* , int , minDist* , curveInfo* , curves* , int , int );
double meanDistanceOfObject(curveInfo* , int* , int , int , double** , int);
double* silhouetteCalculation(cluster* , int , minDist* , curveInfo*, curves *, int , int , double**, double** );

#endif