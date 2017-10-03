#include <math.h>
#include "AHRSstruct.h"
#include "AHRS.h"

void mechanization_step(AHRSstatus *AHRS, float gx, float gy, float gz, float t){
	gx = 0.5 * t *(AHRS -> hat_x_g[0] - gx);
	gy = 0.5 * t *(AHRS -> hat_x_g[1] - gy);
	gz = 0.5 * t *(AHRS -> hat_x_g[2] - gz);
	
	// propagate guaternion
	AHRS -> hat_b[0] = (AHRS -> prev_hat_b[0]) - (AHRS -> prev_hat_b[1]) * gx - (AHRS -> prev_hat_b[2]) * gy - (AHRS -> prev_hat_b[3]) * gz;
	AHRS -> hat_b[1] = (AHRS -> prev_hat_b[0]) * gx + (AHRS -> prev_hat_b[1]) - (AHRS -> prev_hat_b[2]) * gz + (AHRS -> prev_hat_b[3]) * gy;
	AHRS -> hat_b[2] = (AHRS -> prev_hat_b[0]) * gy + (AHRS -> prev_hat_b[1]) * gz + (AHRS -> prev_hat_b[2]) - (AHRS -> prev_hat_b[3]) * gx;
	AHRS -> hat_b[3] = (AHRS -> prev_hat_b[0]) * gz - (AHRS -> prev_hat_b[1]) * gy + (AHRS -> prev_hat_b[2]) * gx + (AHRS -> prev_hat_b[3]);
	
	// normalize quaternion
	const float recip_norm = 1.0 / sqrtf(powf(AHRS -> hat_b[0], 2.0) + powf(AHRS -> hat_b[1], 2.0) + powf(AHRS -> hat_b[2], 2.0) + powf(AHRS -> hat_b[3], 2.0));
	
	AHRS -> hat_b[0] *= recip_norm;
	AHRS -> hat_b[1] *= recip_norm;
	AHRS -> hat_b[2] *= recip_norm;
	AHRS -> hat_b[3] *= recip_norm;
	
	// propagate gyro bias
	
	AHRS -> hat_x_g[0] = (1 + t * lambda_g_x) * (AHRS -> hat_x_g[0]) + c_vector_x;
	AHRS -> hat_x_g[1] = (1 + t * lambda_g_y) * (AHRS -> hat_x_g[1]) + c_vector_y;
	AHRS -> hat_x_g[2] = (1 + t * lambda_g_z) * (AHRS -> hat_x_g[2]) + c_vector_z;
}