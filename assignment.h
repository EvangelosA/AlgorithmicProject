typedef struct minDist{
	double distance;
	int centroid;
	double second_distance;
	int second_centroid;
	int second_centroid_Pos;
} minDist;

typedef struct cluster{
	int id;
	int centroidId;
	int numOfCurves;
	int *crvs;
} cluster;

minDist* loydsAssignment(curves* , int, int*, int, char*, curves*);

minDist* rangeAssignment(curves*, int, grids*, int*, int, char* , ls* , int , double , int , long long int , int , curves*, grids*);


cluster* createCluster(curves*, int*, minDist*, int, int);
void deleteCluster(cluster* , int);
