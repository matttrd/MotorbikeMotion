#include "AHRSstruct.h"

void update_kf_fdk(AHRSstatus *AHRS, float gx, float gy, float gz, float t){
	// remove bias
	gx = AHRS -> hat_x_g[0] - gx;
	gy = AHRS -> hat_x_g[1] - gy;
	gz = AHRS -> hat_x_g[2] - gz;
	
	// Quaternion real part
	AHRS -> KF_Fdk[0][0] = 1.0;
	
	// Quaternion part
	AHRS -> KF_Fdk[1][1] = 1.0;          AHRS -> KF_Fdk[1][2] = - gz * t;		  AHRS -> KF_Fdk[1][3] =  gy * t;
	AHRS -> KF_Fdk[2][1] =  gz * t;      AHRS -> KF_Fdk[2][2] =  1.0;			  AHRS -> KF_Fdk[2][3] =  - gx * t;
	AHRS -> KF_Fdk[3][1] = - gy * t;     AHRS -> KF_Fdk[3][2] = gx * t;           AHRS -> KF_Fdk[3][3] =  1.0;
	
	// Quaternion gyro bias part
	AHRS -> KF_Fdk[1][4] = 1.0 + 0.5 * t;
	AHRS -> KF_Fdk[2][5] = 1.0 + 0.5 * t;
	AHRS -> KF_Fdk[3][6] = 1.0 + 0.5 * t;
	
	// Gyro bias part
	AHRS -> KF_Fdk[4][4] = 1.0 + lambda_g_x * t;
	AHRS -> KF_Fdk[5][5] = 1.0 + lambda_g_y * t;
	AHRS -> KF_Fdk[6][6] = 1.0 + lambda_g_z * t;
}


void update_kf_qdk(AHRSstatus *AHRS, float t){
	
	int i = 0;
	for (i = 0; i < n_states; i++){
		AHRS -> KF_Qdk[i][i] = (AHRS -> KF_Gdk[i]) * (AHRS ->KF_Q[i]) * (AHRS -> KF_Gdk[i]) * t;
	}
}