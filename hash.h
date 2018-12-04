#ifndef HASH__H__
#define HASH__H__

#define HASH_EXTEND 5
#define CONST_RAND 1000
#define LSH_W 4

extern int KAPPA;

typedef struct bucket
{
	int num_of_curves;
	int* curves;
}bucket;

typedef bucket* hash_table;

typedef struct ls
{
	int* rs;
	double* t;
	double** v;
	int w;
	int type;
	hash_table table;
}ls;

hash_table create_hash_table(int);
bucket add_to_bucket(bucket, int);
void print_hash_table(hash_table, int);
void delete_hash_table(hash_table, int);

int *create_rs(int, int);
void delete_rs(int*);

int classic_hash_function(double*, int *, int, long long int m, int);

double* create_lsh_ts(int, int, int);
void delete_lsh_ts(double*);

double* create_normal_vector(int, int);
void delete_normal_vector(double*);

double** create_k_vectors(int, int, int);
void delete_k_vectors(double**, int);

ls* create_ls_array(int, int, int, int, int, int, int);
void delete_ls_array(ls*,int,int,int);

int lsh_h_function(double*, double*, int, double, int);
int lsh_hash_function(double*, ls, int, int, int, long long int, int);

int general_hash_function(double*, ls, int, int, int, long long int, int);

hash_table create_hash_table_with_instances(grids, ls, int, int, int, int, long long int, int );
ls *create_all_hash_tables(grids, ls* , int, int, int, int, long long int, int );

#endif
