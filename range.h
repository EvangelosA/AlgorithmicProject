#ifndef RANGE__H__
#define RANGE__H__

extern int KAPPA;

int same_instances(gridInfo, gridInfo, int, int, int);

int* bucket_range(bucket, curves, grids, curveInfo, gridInfo, int, double, double, int, int*, int);
int* hash_structure_range(curves*, grids*, curveInfo, gridInfo, ls*, int, double, double, int, int, long long int, int);
int** total_range(curves*, grids*, curves*, grids*, ls*, int, double, int, int, long long int, int);
void delete_range_results(int**, int);
int* flatten_order(int*, curves*);

int* range_all_function(curves, curveInfo, double, int);
void print_range(curves, int*);
void delete_check_array(int*);

#endif
