//
//  EKF.c
//  EKF

#include <stdlib.h>
//#include "cblas.h"
#include <Accelerate/Accelerate.h>

#include "EKF.h"

//----------------------------------------- LINEAR ALGEBRA UTILITIES --------------------------------------

static inline void transpose(precision *Ht, int m, int n, precision *H){
	int j, i;
	for (i = 0; i < m; i++){
		for (j = 0; j < n; j++){
			Ht[m * j + i] = H[n * i + j];
		}
	}
}

static inline void symmetrize(int n, precision *M){
	// symmetrizes an n x n square Matrix
	int nn = n*n;
	// allocate space for M^t
	precision *Mtemp = malloc(sizeof(precision) * nn);
	
#if PRECISION == DOUBLE
	cblas_dcopy(nn, M, 1, Mtemp, 1);
#else
	cblas_scopy(nn, M, 1, Mtemp, 1);
#endif
	// compute 0.5*(M + M^t)
	int i,j;
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			M[i*n + j] = 0.5 * (M[i*n + j] + Mtemp[j*n + i]);
		}
	}
	
	free(Mtemp);
	
}
//--------------------------------------------------------------------------------------------------------



EKF *init_ekf(int n, int m, int p, precision *init_state, precision *init_Q, precision *init_R, precision *init_P){
	EKF *ekf = malloc(sizeof(EKF));
	
	ekf->x = calloc(n,sizeof(precision));
	ekf->z = calloc(m, sizeof(precision));
	ekf->u = calloc(p, sizeof(precision));
	//ekf->u_old = calloc(p, sizeof(precision));
	
	ekf->Q = malloc(sizeof(precision) * n * n);
	ekf->R = malloc(sizeof(precision) * m * m);
	ekf->P = malloc(sizeof(precision) * n * n);
	
	ekf->f = malloc(sizeof(precision) * (n + n * n));
	ekf->F = ekf->f + n;
	ekf->h = malloc(sizeof(precision) * (m + m * n));
	ekf->H = ekf->h + m;
	
	ekf->n = n;
	ekf->m = m;
	ekf->p = p;
	
	// initial values
	// state initialization
#if PRECISION == DOUBLE
	cblas_dcopy(n, init_state, 1, ekf -> x, 1);
#else
	cblas_scopy(n, init_state, 1, ekf -> x, 1);
#endif
	
	// initialize Q
#if PRECISION == DOUBLE
	cblas_dcopy(n*n, init_Q, 1, ekf -> Q, 1);
#else
	cblas_scopy(n*n, init_Q, 1, ekf -> Q, 1);
#endif
	
	// initialize R
#if PRECISION == DOUBLE
	cblas_dcopy(m*m, init_R, 1, ekf -> R, 1);
#else
	cblas_scopy(m*m, init_R, 1, ekf -> R, 1);
#endif
	
	// initialize P
#if PRECISION == DOUBLE
	cblas_dcopy(n*n, init_P, 1, ekf -> P, 1);
#else
	cblas_scopy(n*n, init_P, 1, ekf -> P, 1);
#endif
	return ekf;
}

void free_ekf(EKF *ekf){
	
	// free matrices
	free(ekf -> f);
	free(ekf -> h);
	free(ekf -> Q);
	free(ekf -> R);
	free(ekf -> P);
	
	// free states
	free(ekf -> x);
	free(ekf ->z);
	free(ekf->u);
	//free(ekf->u_old);
	
	// free struct
	free(ekf);
}

void ekf_prediction(EKF *ekf, precision dt){
	// local variables to optimize access to the struct
	const int n = ekf -> n;
	const int nn = n * n;
	
	// compute f and F
	ekf->compute_f_F(ekf, dt);

//--------------------------------------- CONDITIONAL UPDATE OF Q ------------------------------------------
#ifndef FIX_Q
	ekf->update_Q(ekf->Q, dt);
#endif
//------------------------------------------------------------------------------------------------------------
	
	// state prediction
#if PRECISION == DOUBLE
	cblas_daxpy(ekf -> n, 1.0, ekf -> f, 1, ekf -> x, 1);
#else
	cblas_saxpy(ekf -> n, 1.0, ekf -> f, 1, ekf -> x, 1);
#endif

	
	
	
	// Pt|t-1 = F Pt-1|t-1 F^t + Q
	// compute FP
	precision *FP = malloc(sizeof(precision) * nn);
#if PRECISION == DOUBLE
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, ekf -> F, n, ekf -> P, n, 0.0, FP, n);
	
#else
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, ekf -> F, n, ekf -> P, n, 0.0, FP, n);
#endif


	// compute FPF^t
#if PRECISION == DOUBLE
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0, FP, n, ekf->F, n, 0.0, ekf->P, n);
#else
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0, FP, n, ekf->F, n, 0.0, ekf->P, n);
#endif
	
	free(FP);

	// add Q
#if PRECISION == DOUBLE
	cblas_daxpy(nn, 1.0, ekf->Q, 1, ekf->P, 1);
#else
	cblas_saxpy(nn, 1.0, ekf->Q, 1, ekf->P, 1);
#endif
	
	// make sure P is still symmetric
	symmetrize(n, ekf->P);
}

void ekf_update(EKF *ekf){
	// local variables to minimize access to the struct
	int m = ekf -> m;
	int n = ekf -> n;
	int mn = m * n;
	int mm = m * m;


	// compute h and H
	ekf->compute_h_H(ekf);
	
//--------------------------------------- CONDITIONAL UPDATE OF R -------------------------------------------
#ifndef FIX_R
	ekf->update_R(ekf);
#endif
//-----------------------------------------------------------------------------------------------------------
	precision *y_tilde = malloc(sizeof(precision) * m);
	precision *z_nav = ekf->z;

	
#if PRECISION == DOUBLE
	cblas_dcopy(m, z_nav, 1, y_tilde, 1);
	cblas_daxpy(m, -1.0, ekf->h, 1, y_tilde, 1);
#else
	cblas_scopy(m, ekf->num_parameters, 1, y_tilde, 1);
	cblas_saxpy(m, -1.0, ekf->h, 1, z_nav, 1);
#endif
	
	// Kalman gain computation
	
	// compute S = HPH^t + R
	// allocate space for HP
	precision *HP = malloc(sizeof(precision) * mn);
	
	// compute HP
#if PRECISION == DOUBLE
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, n, 1.0, ekf -> H, n, ekf -> P, n, 0.0, HP, n);
#else
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, n, 1.0, ekf -> H, n, ekf -> P, n, 0.0, HP, n);
#endif
	
	
	precision *S = malloc(sizeof(precision) * mm);
	
#if PRECISION == DOUBLE
	// S = HPH^t
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, m, n, 1.0, HP, n, ekf -> H, n, 0.0, S, m);
	// S = S + R
	cblas_daxpy(mm, 1.0, ekf->R, 1, S, 1);
#else
	// S = HPH^t
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, m, n, 1.0, HP, n, ekf -> H, n, 0.0, S, m);
	// S = S + R
	cblas_saxpy(mm, 1.0, ekf->R, 1, S, 1);
#endif
	
	free(HP);
	
	// invert S
	int i;
	int *ipivot = malloc(sizeof(int) * (m+1));
	precision *work = malloc(sizeof(precision) * mm);
#if PRECISION == DOUBLE
	// LU decomoposition of a general matrix
	dgetrf_(&m, &m, S, &m, ipivot, &i);
	// generate inverse of a matrix given its LU decomposition
	dgetri_(&m, S, &m, ipivot, work, &mm, &i);

#else
	// LU decomoposition of a general matrix
	sgetrf_(&m, &m, S, &m, ipivot, &i);
	// generate inverse of a matrix given its LU decomposition
	sgetri_(&m, S, &m, ipivot, work, &mm, &i);
#endif
	free(ipivot);
	free(work);

	// K = PH^t * S^-1 WARNING ekf -> S is already S^-1
	precision *PHt = malloc(sizeof(precision) * mn);
	memset(PHt, 0, mn);
	precision *K = malloc(sizeof(precision) * mn);
	
	// transpose H
	precision *Ht = malloc(sizeof(precision) * mn);
	transpose(Ht, m, n, ekf->H);
	
#if PRECISION == DOUBLE
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, m, n, 1.0, ekf->P, n, Ht, m, 0.0, PHt, m);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, m, m, 1.0, PHt, m, S, m, 0.0, K, m);

#else
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, m, n, 1.0, ekf->P, n, Ht, m, 0.0, PHt, m);
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, m, m, 1.0, PHt, m, S, m, 0.0, K, m);
#endif
	free(Ht);
	free(PHt);
	free(S);
	
	// state update xt|t = xt|t - 1 + K y_tilde
#if PRECISION == DOUBLE
	cblas_dgemv(CblasRowMajor, CblasNoTrans, n, m, 1.0, K, m, y_tilde, 1, 1.0, ekf->x, 1);
#else
	cblas_sgemv(CblasRowMajor, CblasNoTrans, n, m, 1.0, K, m, y_tilde, 1, 1.0, ekf->x, 1);
#endif


	free(y_tilde);
	
	// update P
	//allocate space for -KH
	int nn = n*n;
	precision *KH = malloc(sizeof(precision)*nn);
#if PRECISION == DOUBLE
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, m, -1.0, K, m, ekf->H, n, 0.0, KH, n);
	free(K);
	for (i = 0; i < n; i++){
		KH[(n + 1)* i] = KH[(n + 1)* i] + 1.0;
	}
#else
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, m, -1.0, K, m, ekf->H, n, 0.0, KH, n);
	free(K);
	for (i = 0; i < n; i++){
		KH[(n + 1)* i] = KH[(n + 1)* i] + 1.0;
	}
#endif
	// allocate space for (I - KH)P 
	precision * tmp = malloc(sizeof(precision) * nn);
#if PRECISION == DOUBLE
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, KH, n, ekf->P, n, 0.0, tmp, n);
	cblas_dcopy(nn, tmp, 1, ekf->P, 1);
#else
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, KH, n, ekf->P, n, 0.0, tmp, n);
	cblas_scopy(nn, tmp, 1, ekf->P, 1);
#endif
	free(KH);
	free(tmp);
	symmetrize(n, ekf->P);
}
