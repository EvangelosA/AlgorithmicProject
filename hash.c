#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "grid.h"
#include "hash.h"

hash_table create_hash_table(int length) //Creating the hash table
{
	hash_table new_h;
	new_h = (hash_table)malloc(length * sizeof(bucket));
	int i;
	for( i=0 ; i<=length-1 ; i++)
	{
		new_h[i].curves = NULL;
		new_h[i].num_of_curves = 0;
	}
	return new_h;
}

bucket add_to_bucket(bucket b, int new_num) //Adding to hash table a new curve
{
	if(b.num_of_curves == 0)
	{
		b.curves = (int *)malloc(HASH_EXTEND * sizeof(int));
	}
	else if(b.num_of_curves%HASH_EXTEND == 0)
	{
		b.curves = (int *)realloc(b.curves,(HASH_EXTEND + b.num_of_curves) * sizeof(int));
		if(b.curves == NULL)
		{
			fprintf(stderr,"Realloc for hash bucket extension failed\n");
			b.num_of_curves = 0;
			b.curves = NULL;
			return b;
		}
	}
	b.num_of_curves++;
	b.curves[b.num_of_curves-1] = new_num;
	return b;
}

void print_hash_table(hash_table table, int len) //Print hash table - Test function
{
	int i,j;
	printf("Printing Hash table:\n");
	for(i=0; i<=len-1; i++)
	{
		printf("\tPos %d:\n\t\t",i);
		for(j=0; j<=table[i].num_of_curves-1; j++)
		{
			printf("%d,",table[i].curves[j]);
		}
	}
}

void delete_hash_table(hash_table t, int length) //Deleting the table
{
	int i;
	for( i=0; i<=length-1; i++)
	{
		if(t[i].curves != NULL)
		{
			free(t[i].curves);
		}
	}
	free(t);
}

int* create_rs(int seed, int size) //Creating an array of random rs. For both LSH and classic
{
	int* new_rs;
	int i;
	new_rs = (int*)malloc(size*sizeof(int));
	for(i=0; i<=size-1; i++)
	{
		new_rs[i] = rand()%CONST_RAND;
	}
	return new_rs;
}

void delete_rs(int *rs)
{
	free(rs);
}

int classic_hash_function(double* instance, int* rs, int length, long long int m, int table_size) //This is the classic hash function
{
	double sum;
	int i;
	sum = 0.0;
	for(i=0; i<=length-1; i++)
	{
		if(instance[i]>0)
		{
			sum = sum + instance[i]*rs[i];
		}
		else
		{
			sum = sum - instance[i]*rs[i];
		}
	}
	return ((int)sum%m)%table_size;
}

double* create_lsh_ts(int seed, int k, int w) //Creating random t's for lsh functions
{
	double* new_ts;
	new_ts = (double*)malloc(k * sizeof(double));
	int i;
	for(i=0; i<=k-1; i++)
	{
		new_ts[i] = (((double)(rand()%CONST_RAND))/CONST_RAND)*w;
	}
	return new_ts;
}

void delete_lsh_ts(double* ts)
{
	if(ts!=NULL)
	{
		free(ts);
	}
}

double* create_normal_vector(int seed, int dim) //Creating a normal vector for h of LSH
{
	double u,v,s;
	int i, check;
	double* vector;
	vector = (double*)malloc(dim * sizeof(double));
	for(i=0; i<=dim-1; i=i+2)
	{
		check=1;
		do
		{
			do
			{
				u = ((double)(rand()%CONST_RAND)/CONST_RAND)*2-1;
				v = ((double)(rand()%CONST_RAND)/CONST_RAND)*2-1;
				s = u*u+v*v;
			}while(s>=1 || s<=0 || u<=-1 || u>=1 || v<=-1 || v>=1);

			vector[i] = (u*sqrt(-2*log(s)/s)+3)/6;
			if(vector[i]>=1 || vector[i]<=0)
				check=1;
			else
				check=0;
			if(dim%2 == 0 || i != dim-1)
			{
				vector[i+1] = (v*sqrt(-2*log(s)/s)+3)/6;
				if(vector[i+1]>=1 || vector[i+1]<=0)
					check=1;
			}
		}while(check);
	}
	return vector;
}

void delete_normal_vector(double* v)
{
	free(v);
}

double** create_k_vectors(int seed, int k, int dim) //Create all the k vectors for LSH
{
	double** vectors;
	vectors = (double**)malloc(k * sizeof(double*));
	int i;
	for(i=0; i<=k-1; i++)
	{
		vectors[i] = create_normal_vector(seed, dim);
	}
	return vectors;
}

void delete_k_vectors(double** vectors, int k)
{
	int i;
	if(vectors != NULL)
	{
		for(i=0; i<=k-1; i++)
		{
			delete_normal_vector(vectors[i]);
		}
		free(vectors);
	}
}

ls* create_ls_array(int type, int seed, int k, int dim, int l, int max_points, int table_size) //Creating the structure for the ls struct
{
	ls* new_array;
	new_array = (ls*)malloc(l * sizeof(ls));
	int i;
	for(i=0; i<=l-1; i++)
	{
		if(type == 1)
		{
			new_array[i].type = 1;
			new_array[i].t = NULL;
			new_array[i].v = NULL;
			new_array[i].w = 0;
			new_array[i].rs = create_rs(seed, dim*max_points*KAPPA);
			new_array[i].table = create_hash_table(table_size);
		}
		else
		{
			new_array[i].type = 2;
			new_array[i].t = create_lsh_ts(seed, k, LSH_W);
			new_array[i].v = create_k_vectors(seed, k, max_points*dim*KAPPA);
			new_array[i].w = LSH_W;
			new_array[i].rs = create_rs(seed, k);
			new_array[i].table = create_hash_table(table_size);
		}
	}
	return new_array;
}

void delete_ls_array(ls* array, int l, int k,int table_size)
{
	int i;
	for(i=0;i<=l-1;i++)
	{
		delete_lsh_ts(array[i].t);
		delete_k_vectors(array[i].v, k);
		delete_rs(array[i].rs);
		delete_hash_table(array[i].table, table_size);
	}
	free(array);
}

int lsh_h_function(double* instance, double* v, int length, double t, int w) //The h function for LSH
{
	int i;
	double sum;
	sum=0;
	for(i=0; i<=length-1; i++)
	{
		sum = sum + instance[i]*v[i];
	}
	sum = sum + t;
	return (int)sum/w;
}

int lsh_hash_function(double* instance, ls l, int k, int max_points, int dim, long long int m, int table_size) //The total LSH hash function
{
	int i, sum;
	sum = 0;
	for(i=0; i<=k-1; i++)
	{
		sum = sum + lsh_h_function(instance, l.v[i], max_points*KAPPA*dim, l.t[i], l.w);
	}
	sum = (sum%m)%table_size;
	while(sum<0)
	{
		sum=sum+table_size;
	}
	return sum;
}

int general_hash_function(double* instance, ls l, int k, int maxpoints, int dim, long long int m, int table_size) //General function for both LSH and Classic
{
	if(l.type == 1)
	{
		return classic_hash_function(instance, l.rs, maxpoints*dim*KAPPA, m, table_size);
	}
	else if(l.type == 2)
	{
		return lsh_hash_function(instance, l, k, maxpoints, dim, m, table_size);
	}
	else
	{
		fprintf(stderr, "Hash function must be 1 or 2\n");
		return -1;
	}
}

hash_table create_hash_table_with_instances(grids grd, ls l, int k, int max_points, int dim, int specific_l, long long int m, int table_size) //Creating one hash table
{
	int i, value, asd;
	hash_table n_t;
	n_t = l.table;
	for(i=0; i<=grd.numOfGrids-1; i++)
	{
		value = general_hash_function(grd.gridArray[i].grid_instances[specific_l], l, k, max_points, dim, m, table_size);

		if(value >= 1870){
			asd = 0;
		}

		n_t[value] = add_to_bucket(n_t[value], i);
	}
	return n_t;
}

ls* create_all_hash_tables(grids grd, ls* l_array, int l, int k, int max_points, int dim, long long int m, int table_size) //Creating all hash tables
{
	int i;
	for(i=0; i<=l-1; i++)
	{
		l_array[i].table = create_hash_table_with_instances(grd, l_array[i], k,max_points, dim, i, m, table_size);
	}
	return l_array;
}
