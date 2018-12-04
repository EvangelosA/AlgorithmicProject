#ifndef INITIALIZATION__H__
#define INITIALIZATION__H__

#define ERR_NO_NUM -1
#define ERR_NO_MEM -2

typedef struct plusplusInfo
{
	double minDistanceToCentroid;		//D(i)
	double minDistanceToCentroidSqrt;	//D(i)^2
	int closestCentroid;
	int isCentroid; 					//0:NO, 1:YES
	double partialSum;					//P(r)
}plusplusInfo;

typedef struct pSum
{
	int originalPos;
	double partialSum;
}pSum;

int myRandom(int, int);
int* randomSelection(curves*, int, int );
int binarySearch(pSum* , double , int , int);
int* kmeansplusplus(curves*, int, int , char* );

#endif
