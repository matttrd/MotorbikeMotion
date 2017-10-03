#include <stdlib.h>
#include <stdio.h>

// Functions declarations
float **malloc_matrix(int rows, int cols);

void free_matrix(int rows, float **mat);

void sum_vectors(float *u, float *v, int l, float *w);

void multiply_matrix_vector(float **matrix, float *v, int nrows, int ncols, int l_v, float *w);

void multiply_matrices(float **A, int nrows_A, float **B, int ncols_B , int ncols_A_nrows_B, float **output);