#include "conversions.h"
#include <stdlib.h>
#include <math.h>
#include <Accelerate/Accelerate.h>

void euler_angles_from_acc_and_mag(precision ax, precision ay, precision az, precision mx, precision my, precision mz, precision *ret){
	// this function assumes that the navigation frame is NED
	// the input is IMU data and the output is specified by the pointer ret
#if PRECISION == DOUBLE
	const precision roll = atan(ay/az);
	const precision pitch = atan2(ax, sqrtf(powf(ay, 2) + powf(az, 2)));
	
	// store sine and cosine of roll and pitch to optimize computations
	const precision c_roll = cos(roll);
	const precision s_roll = sin(roll);
	
	const precision c_pitch = cos(pitch);
	const precision s_pitch = sin(pitch);
#else
	const precision roll = atanf(ay/az);
	const precision pitch = atan2f(ax, sqrtf(powf(ay, 2) + powf(az, 2)));
	
	// store sine and cosine of roll and pitch to optimize computations
	const precision c_roll = cosf(roll);
	const precision s_roll = sinf(roll);
	
	const precision c_pitch = cosf(pitch);
	const precision s_pitch = sinf(pitch);
#endif
	
	precision *mat = malloc(sizeof(precision)*9);
	
	mat[0]= c_pitch;	 mat[1] = s_pitch * s_roll;   mat[2] = s_pitch * c_roll;
	mat[3] = 0.0;		 mat[4] = c_roll;             mat[5] = - s_roll;
	mat[6] = - s_pitch;  mat[7] = c_pitch * s_roll;   mat[8] = c_pitch * c_roll;

	// allocate pointer for mult matrix
	precision *mag = malloc(sizeof(precision) * 3);
	mag[0] = mx;
	mag[1] = my;
	mag[2] = mz;
	
	precision *bar_m_w = malloc(sizeof(precision) * 3);
	
#if PRECISION == DOUBLE
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, mat, 3, mag, 1, 0.0, bar_m_w, 1);
	free(mag);
	free(mat);
	
	// compute yaw
	ret[2] = atan2(- bar_m_w[0], bar_m_w[1]);
	free(bar_m_w);
#else
	cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, mat, 3, mag, 1, 0.0, bar_m_w, 1);
	free(mag);
	free(mat);
	
	// compute yaw
	ret[2] = atan2f(- bar_m_w[0], bar_m_w[1]);
	free(bar_m_w);
#endif
	
	// prepare roll and pitch for return
	ret[0] = roll;
	ret[1] = pitch;
}

void euler_angles_to_rotation(precision phi, precision theta, precision psi, precision *mat){
	
	const precision c_phi = cosf(phi);
	const precision s_phi = sinf(phi);
	
	const precision c_theta = cosf(theta);
	const precision s_theta = sinf(theta);
	
	const precision c_psi = cosf(psi);
	const precision s_psi = sinf(psi);
	
	mat[0] = c_psi * c_theta;						mat[1] = s_psi*c_theta;						  mat[2] = -s_theta;
	mat[3] = - s_psi * c_phi + c_psi*s_theta*s_phi;  mat[4] = c_psi*c_phi + s_psi*s_theta*s_phi;    mat[5] =  c_theta*s_phi;
	mat[6] = s_psi * s_phi + c_psi*s_theta*c_phi;	mat[7] = -c_psi*s_phi + s_psi*s_theta*c_phi;   mat[8] = c_theta*c_phi;
	
}

void rotation_to_quaternion(precision *R, precision *ret){
	precision trace = R[0] + R[4] + R[8];

	if (trace < 0){
		ret[1] = 0.5 * sqrtf(1 + R[0] - R[4] - R[8]);
		ret[0] = (R[2*3 +1] - R[1*3 + 2]) / (4 * ret[1]);
		ret[2] = (R[0 * 3 + 1] + R[1 * 3 + 0]) / (4 * ret[1]);
		ret[3] = (R[2 * 3 + 0] + R[0 * 3 + 2]) / (4 * ret[1]);
	} else{
		ret[0] = 0.5 * sqrtf(1 + trace);
		ret[1] = (R[2*3 + 1] - R[1*3 + 2]) / (4 * ret[0]);
		ret[2] = (R[0 * 3 + 2] - R[2*3 + 0]) / (4 * ret[0]);
		ret[3] = (R[1 * 3 + 0] - R[0 * 3 + 1]) / (4 * ret[0]);
	}
}
