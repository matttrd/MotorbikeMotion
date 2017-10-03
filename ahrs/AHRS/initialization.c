#include <stdlib.h>
#include <Accelerate/Accelerate.h>

#include “../../EKF/EKFstruct.h"
#include "driver_functions.h"
#include "conversions.h"
//#include "linalg.h"
#include <math.h>


precision bar_g_b[3] = { 0.0 };

precision bar_m_b[3] = { 0.0 };

int gps_counter = 1;

void initialization_stage_estimates(EKF *ekf, precision* measurements, int init_counter){
	
	
	precision gx = measurements[6];
	precision gy = measurements[7];
	precision gz = measurements[8];
	precision ax = measurements[9];
	precision ay = measurements[10];
	precision az = measurements[11];
	precision mx = measurements[12];
	precision my = measurements[13];
	precision mz = measurements[14];
	
	const int nsamples = init_counter + 1;
	const precision recip_nsamples = 1.0 / nsamples;
	
	// Update average
	precision *hat_x_g = ekf->x + 10;
	hat_x_g[0] = ((hat_x_g[0] * (nsamples - 1)) + gx) * recip_nsamples;
	hat_x_g[1] = ((hat_x_g[1] * (nsamples - 1)) + gy) * recip_nsamples;
	hat_x_g[2] = ((hat_x_g[2] * (nsamples - 1)) + gz) * recip_nsamples;
	
	bar_g_b[0] = ((bar_g_b[0] * (nsamples - 1)) + ax) * recip_nsamples;
	bar_g_b[1] = ((bar_g_b[1] * (nsamples - 1)) + ay) * recip_nsamples;
	bar_g_b[2] = ((bar_g_b[2] * (nsamples - 1)) + az) * recip_nsamples;
	
	bar_m_b[0] = ((bar_m_b[0] * (nsamples - 1)) + mx) * recip_nsamples;
	bar_m_b[1] = ((bar_m_b[1] * (nsamples - 1)) + my) * recip_nsamples;
	bar_m_b[2] = ((bar_m_b[2] * (nsamples - 1)) + mz) * recip_nsamples;
	
	precision *euler = malloc(3 * sizeof(precision));
	euler_angles_from_acc_and_mag(bar_g_b[0], bar_g_b[1], bar_g_b[2], bar_m_b[0], bar_m_b[1], bar_m_b[2], euler);
	
	precision *rotation = malloc(sizeof(precision) * 9);
	euler_angles_to_rotation(euler[0], euler[1], euler[2], rotation);
	free(euler);
	
	precision *quat = malloc(sizeof(precision) * 4);
	rotation_to_quaternion(rotation, quat);
	free(rotation);
	
	ekf->x[6] = quat[0];
	ekf->x[7] = quat[1];
	ekf->x[8] = quat[2];
	ekf->x[9] = quat[3];
	
	normalize_quaternion(ekf);
	free(quat);
	
	if (measurements[0] < 1000){
		precision recip_gps_counter = 1.0/gps_counter;
#if PRECISION == DOUBLE
		cblas_dscal(6, gps_counter - 1, ekf->x, 1);
		cblas_daxpy(6, 1.0, measurements, 1, ekf->x, 1);
		cblas_dscal(6, recip_gps_counter, ekf->x, 1);
#else
		cblas_dscal(6, gps_counter - 1, ekf->x, 1);
		cblas_daxpy(6, 1.0, measurements, 1, ekf->x, 1);
		cblas_dscal(6, recip_gps_counter, ekf->x, 1);
#endif
		gps_counter++;
	}
};

void GNSS_INS_initialization(EKF *ekf,GNSS_INS_t *ins_struc){

	precision *Q = ekf->Q;
	precision *P = ekf->P;
	
	Q[0] = Q_LAT;
	Q[13 + 1] = Q_LON;
	Q[26 + 2] = Q_H;
	Q[39 + 3] = Q_V;
	Q[52 + 4] = Q_V;
	Q[65 + 5] = Q_V;
	Q[78 + 6] = Q_QUAT;
	Q[91 + 7] = Q_QUAT;
	Q[104 + 8] = Q_QUAT;
	Q[117+ 9] = Q_QUAT;
	Q[130 + 10] = Q_BIAS;
	Q[143 + 11] = Q_BIAS;
	Q[156 + 12] = Q_BIAS;
	
	
	precision *R = ekf->R;
	
	R[0] = R_LAT;
	R[12 + 1] = R_LON;
	R[24 + 2] = R_H;
	R[36 + 3] = R_V;
	R[48 + 4] = R_V;
	R[60 + 5] = R_V;
	
	
	R[72 + 6] = R_ACC_X;
	R[84 + 7] = R_ACC_Y;
	R[96 + 8] = R_ACC_Z;
	R[108 + 9] = R_MAG_X;
	R[120 + 10] = R_MAG_Y;
	R[132 + 11] = R_MAG_Z;
	
	//=================================== FORCE INITIALIZATION========================================
	//precision x0[] = {45.4170, 11.8669, 172.0583, 0, 0, 0, 1, 0, 0, 0, 0.0038, 0.0061, -0.0219};
	//cblas_dcopy(13, &x0[0], 1, ekf->x, 1);
	//================================================================================================

	
	int i;
	for (i = 0; i < N_STATES - 3; i++){
		P[i * (N_STATES + 1)] += 1;
	}
	P[10] = 0.00001;
	P[11] = 0.00001;
	P[12] = 0.00001;
    
    
    
    ekf->z[9] = cos(ins_struc->mag_variation*M_PI/180);
    ekf->z[10] = sin(ins_struc->mag_variation*M_PI/180);
}