
#include "AHRS.h"

//----------------------------------------------------------------------------------------------------------------------------------
// EKF struct
struct ekf_ahrs_status {
	// Initialization
	float bar_g_b[3];
	float bar_m_b[3];
	float bar_x_g[3];			// gyro bias estimate
	float nu_g[3];					// gyro measurements variance of initial samples
	
	float gyro_init_list[init_stop_counter][3];
	float acc_init_list[init_stop_counter][3];
	float mag_init_list[init_stop_counter][3];
	
	// Mechanization matrices
	float KF_Fdk[n_states][n_states];
	
	// AR1 models
	//		for gyro
	const float sigma_omega_g[3];
	const float F_g[3][3];
	const float c_vector[3];

	
	// Kalman filter matrices
	float KF_P_pred[n_states][n_states];
	float KF_P_update[n_states][n_states];
	float KF_Q[n_states];				// SIMP we save only the diagonal of Q
	float KF_Gdk[n_states];			// SIMP we save only the diagonal of Gdk
	float KF_Qdk[n_states][n_states];
	float KF_R[6][6];
	
	// Kalman filter state (indirect)
	float delta_x_g[3];
	float delta_x[4];
	
	// Estimated quaternions and gyro bias
	float hat_b[4];				// Current est. quaternion
	float hat_x_g[3];			// Current est. gyroscope bias
	float prev_hat_b[4];		// Previous est. quaternion
	float prev_hat_x_g[3];		// Previous est. gyroscope bias
	};
typedef struct ekf_ahrs_status AHRSstatus;