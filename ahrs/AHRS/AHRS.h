#ifndef AHRS_h
#define AHRS_h
// Definitions
# define dt							0.05
# define init_samples				100
# define INIT_STOP_COUNTER			init_samples - 1
#define recipNsamples				1.0f / init_samples;
#define recip_INIT_STOP_COUNTER		1.0f / INIT_STOP_COUNTER;


// Magnetic flux components
# define m_n_x				22.2797
# define m_n_y				1.0626
# define m_n_z				42.0277

// AR1 for gyro

# define sigma_omega_g_x	1.917E-07
# define sigma_omega_g_y	2.1793E-07
# define sigma_omega_g_z	2.8224E-07

# define lambda_g_x			-49.6392
# define lambda_g_y			-48.4244
# define lambda_g_z			-49.5950

# define c_vector_x			3.3763E-04
# define c_vector_y			2.1855E-04
# define c_vector_z			-6.6856E-05

// Validity thesholds
# define beta_mu			0.5
# define beta_gyro			5.0
# define INIT_MU			1.0

// Kalman Filter design
# define n_states			7

// Function declarations


void EKF_AHRS_update(float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz);
#endif