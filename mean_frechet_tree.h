int pickCurveRandomly(int *, int );
int *initializeStartLevel(cluster );
curves *initializeFirstLevel(curves *,int *,int );
curves *produceLevel(curves *);
curveInfo *chooseCentroid(curves *crvs,cluster cls);
curves *meanFrechetCentroids(curves *,cluster *,int);
