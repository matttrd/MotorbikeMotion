//
//  INS_matrix_fun.h

#ifndef INS_matrix_fun_h
#define INS_matrix_fun_h


#include <stdio.h>
#include <math.h>
//#include <stdbool.h>
#include <stdlib.h>
#include "quaternion.h"
#include "EKFstruct.h"



//---------------------------------------- KALMAN FILTER PARAMETERS -----------------------------------------

#define EARTH_RADIUS 6378137
#define FLATTENING 1/298.257223563 //not necessary
#define ECC (FLATTENING*(2-FLATTENING))
#define OMEGA_IE 7.292115e-5
#define ge 9.81
// magnetic flux reference
#define	MAG_REF_N	0.9983
#define	MAG_REF_E	-0.0584

// Kalman filter dimensions
#define N_STATES 13
#define N_INPUTS 15
#define N_OBS	12

#define N_STATES_2 169
#define N_OBS_2	144

// Matrix Q diagonal entries
#define Q_LAT	0.0001
#define Q_LON	0.0001
#define Q_H		0.1

#define Q_V		1.0

#define Q_QUAT	1e-6

#define Q_BIAS	1e-9


// Matrix R diagonal entries


#define R_LAT	1e-10
#define R_LON	1e-10
#define R_H		1e-6

#define R_GPS_BAD	1e+30

#define	R_V		1.0
#define R_V_BAD 1e+30

#define R_ACC_X 0.1/*0.0857*/
#define R_ACC_Y 0.05/*0.0266*/
#define R_ACC_Z 0.1610
#define R_ACC_BAD 10e+30

#define R_MAG_X 0.0052
#define R_MAG_Y 2.1423
#define R_MAG_Z 0.8757

// Bad acc check
#define BETA_MU	0.5 //acc
#define BETA_GYRO 5.0 //gyro

precision invSqrt(precision x);

void compute_f_and_F(EKF *ekf, precision dt);

void compute_INS_h_and_H(EKF *ekf);

void normalize_quaternion(EKF *ekf);

void update_R(EKF *ekf);

#endif /* INS_matrix_fun_h */
