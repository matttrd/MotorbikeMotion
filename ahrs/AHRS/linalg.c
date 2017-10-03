#include "linalg.h"
#include "stdlib.h"

// Functions body

float **malloc_matrix(int rows, int cols){
	
	float **mat = (float **) malloc(sizeof(float *)*rows);
	if (NULL == mat){
		printf("Could not allocate matrix");
		abort();
	};
	
	int i=0;
	for(i=0; i<rows; i++){
	/* Allocate array, store pointer  */
		mat[i] = (float *) malloc(sizeof(float)*cols);
		if (NULL == mat[i]){
			printf("Could not allocate matrix");
			abort();
		};
	};
	return mat;
}

void free_matrix(int rows, float **mat){
	int i=0;
	for(i=0;i<rows;i++)
		free(mat[i]);
	free(mat);
}



void sum_vectors(float *u, float *v, int l, float *w){
	int i;
	for (i = 0; i < l; i++){
		w[i] = u[i] + v[i];
	}
}

void multiply_matrix_vector(float **matrix, float *v, int ncols, int nrows, int l_v, float *w){
	
	if (ncols != l_v){
		printf("The dimensions of the vector and the number of columns must agree");
		abort();
	}

	int i = 0, j = 0;
	for (i = 0; i < nrows; i++){
		w[i] = 0;
		for (j = 0; j < ncols; j++){
			w[i] += matrix[i][j] * v[j];
		}
	} 
}

float **transpose(float **matrix, int nrows, int ncols){
	float **transpose = malloc_matrix(ncols, nrows);
	
	int i = 0; int j = 0;
	for (i = 0; i < ncols; i++){
		for (j  = 0; j < nrows; j++){
			transpose[i][j] = matrix[j][i];
		}
	}
	return transpose;
}

void multiply_matrices(float **A, int nrows_A, float **B, int ncols_B , int ncols_A_nrows_B, float **output){
	int i = 0; int j = 0; int k = 0;
	for (i = 0; i < nrows_A; i++){
		output[i][j] = 0;
		for (j = 0; j < ncols_B; j++){
			for (k = 0; k < ncols_A_nrows_B; k++){
				output[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	
}