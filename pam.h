#ifndef PAM__H__
#define PAM__H__

double** createAndInitialize2D(int );
void delete2D(double** , int );
double calculateObjectiveFunction(curveInfo* , int* , int , int , double** );
int* PartitionAroundMedoids(cluster*, int, minDist*, int, curveInfo*, double** );

#endif
