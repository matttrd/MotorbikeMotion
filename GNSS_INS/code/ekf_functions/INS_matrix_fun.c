//
//  INS_matrix_fun.c

#include "INS_matrix_fun.h"



precision invSqrt(precision x) {
	precision y;
#if PRECISION == DOUBLE
	y = 1.0/ sqrt(x);
#else
	float halfx = 0.5f * x;
	float y = x;
	long i = *(long*)&y;
	i = 0x5f3759df - (i>>1);
	y = *(float*)&i;
	y = y * (1.5f - (halfx * y * y));
#endif
    return y;
}

// here remember the sign of ge!!!!!!!!!!!!!!!!!!! It must be ge = 9.81.

void compute_f_and_F(EKF *ekf, precision dt)
{
	
	
	precision *fF = ekf->f;
	precision *x = ekf->x;
	precision *z = ekf->u + 6;
   // precision *fF = (precision*) malloc(sizeof(precision)*182);
    
    // Define Matrix as an array
    
    // Here z represents u(t), that is, [w_x,w_y,w_z],[a_x,a_y,a_z],[lambda_gx,lambda_gy,lambda_gz] in this order
    precision *prev_measurement_gyro = &z[0];
    precision *prev_measurement_acc  = &z[3];
	precision *lambda_g              = calloc(3, sizeof(precision));
		
    precision omega_ie = OMEGA_IE;
    precision phi = x[0];
    precision lam = x[1];
    precision h   = x[2];
    precision v_n = x[3];
    precision v_e = x[4];
    precision v_d = x[5];
    
    precision q0 = x[6];
    precision q1 = x[7];
    precision q2 = x[8];
    precision q3 = x[9];
    precision x_g1 = x[10];
    precision x_g2 = x[11];
    precision x_g3 = x[12];
    
    precision omega_ib_b1 = prev_measurement_gyro[0];
    precision omega_ib_b2 = prev_measurement_gyro[1];
    precision omega_ib_b3 = prev_measurement_gyro[2];
    
    
    precision lambda_xg_1 = lambda_g[0];
    precision lambda_xg_2 = lambda_g[1];
    precision lambda_xg_3 = lambda_g[2];
    
    quaternion_t *q = quaternion_new_set(q0,q1,q2,q3);
    quaternion_t *q_conj = quaternion_new();
    quaternConj(q_conj,q);
    quaternion_t *f_nav = quaternion_new();
    precision accx = prev_measurement_acc[0];
    precision accy = prev_measurement_acc[1];
    precision accz = prev_measurement_acc[2];
    quaternion_t *q_acc = quaternion_new_set(0,
                                             accx,
                                             accy,
                                             accz);
    quaternProd(f_nav,q_conj,q_acc);
    quaternion_t *f_nav_1 = quaternion_new();
    quaternProd(f_nav_1,f_nav,q);
    precision f_n = f_nav_1->q[1];
    precision f_e = f_nav_1->q[2];
    precision f_d = f_nav_1->q[3];
	

    
    // Compute R_M,R_N
    precision sin_lam = sin(M_PI/180.0 * lam);
    precision sin_phi = sin(M_PI/180.0 * phi);
    precision cos_phi = cos(M_PI/180.0 * phi);
    precision cos_phi_2 = cos_phi*cos_phi;
    precision sin_phi_2 = sin_phi*sin_phi;
    //precision cos_lam = cos(lam);
    
    precision ECC_2 = ECC*ECC;
    precision R_M = EARTH_RADIUS*(1-ECC_2)/pow((1-ECC_2*sin_lam*sin_lam),1.5);
    precision R_N = EARTH_RADIUS*invSqrt(1-ECC_2*sin_lam*sin_lam);
    
    precision R_N_plus_h = (R_N + h);
    precision R_M_plus_h = (R_M + h);
    precision R_N_plus_h_2 = R_N_plus_h*R_N_plus_h;
    precision R_M_plus_h_2 = R_M_plus_h*R_M_plus_h;
    precision v_e_2 = v_e*v_e;
    precision R_N_plus_h_dot_cos_phi = R_N_plus_h*cos_phi;
    precision R_N_plus_h_dot_cos_phi_2 = (R_N_plus_h*cos_phi_2);
    precision R_N_plus_h_2_dot_cos_phi = (R_N_plus_h_2*cos_phi);
    precision omega_ie_dot_2 = 2*omega_ie;
    precision vn_frac_R_M_plus_h = v_n/R_M_plus_h;
    precision v_e_frac_R_N_plus_h_dot_cos_phi = v_e/R_N_plus_h_dot_cos_phi;
    precision diff1 = -omega_ie_dot_2 - v_e_frac_R_N_plus_h_dot_cos_phi;
    precision sum1 = omega_ie_dot_2 + v_e_frac_R_N_plus_h_dot_cos_phi;
    precision q0_dot_sum_1_dot_cos_phi = q0*sum1*cos_phi;
    precision diff1_dot_sin_phi = (diff1)*sin_phi;
    precision q2_dot_diff_1_dot_sin_phi = q2*diff1_dot_sin_phi;
    precision sum2 = q0_dot_sum_1_dot_cos_phi + q2_dot_diff_1_dot_sin_phi + q3*vn_frac_R_M_plus_h;
    precision sum6 =(-q1*sum1*cos_phi + q2*vn_frac_R_M_plus_h - q3*diff1_dot_sin_phi);
    precision sum7 = (q0*diff1_dot_sin_phi - q1*vn_frac_R_M_plus_h - q2*sum1*cos_phi);
    precision sum3 = sum6 + q2*sum7;
    precision sum4 = (-q0*vn_frac_R_M_plus_h - q1*diff1_dot_sin_phi + q3*sum1*cos_phi);
    precision sum5 = -omega_ib_b1 + q0*sum2 - q1*sum3 - q3*sum4 + x_g1;
    //precision q0_2 = q0*q0; // qui secondo me siamo al limite, meglio fare conti che usare memoria per poche operazione
	precision inv_PI_180 = 180.0/M_PI;
	
    // fF[0]-[12] = f;
    //fF[13]-[end] = F;
    
    fF[0] = inv_PI_180 * vn_frac_R_M_plus_h;
    fF[1] = inv_PI_180 * v_e_frac_R_N_plus_h_dot_cos_phi;
    fF[2] = -v_d;
    fF[3] = f_n + v_d*vn_frac_R_M_plus_h + v_e*diff1_dot_sin_phi;
    fF[4] = f_e + v_d*sum1*cos_phi - v_n*diff1_dot_sin_phi;
    fF[5] = f_d + ge - v_e*sum1*cos_phi - v_n*vn_frac_R_M_plus_h;
    fF[6] = -0.5*q1*(sum5) - 0.5*q2*(-omega_ib_b2 + q0*sum4 - q1*sum7 - q2*sum6 + q3*sum2 + x_g2) - 0.5*q3*(-omega_ib_b3 + q0*sum7 + q1*sum4 - q2*sum2 - q3*sum6 + x_g3);
    
    fF[7] =  0.5*q0*(sum5) - 0.5*q2*(-omega_ib_b3 + q0*sum7 + q1*sum4 - q2*sum2 - q3*sum6 + x_g3) + 0.5*q3*(-omega_ib_b2 + q0*sum4 - q1*sum7 - q2*sum6 + q3*sum2 + x_g2);
    
    fF[8] = 0.5*q0*(-omega_ib_b2 + q0*sum4 - q1*sum7 - q2*sum6 + q3*sum2 + x_g2) + 0.5*q1*(-omega_ib_b3 + q0*sum7 + q1*sum4 - q2*sum2 - q3*sum6 + x_g3) - 0.5*q3*(sum5);
    
    fF[9] = 0.5*q0*(-omega_ib_b3 + q0*sum7 + q1*sum4 - q2*sum2 - q3*sum6 + x_g3) - 0.5*q1*(-omega_ib_b2 + q0*sum4 - q1*sum7 - q2*sum6 + q3*sum2 + x_g2) + 0.5*q2*(sum5);
    fF[10] = lambda_xg_1*x_g1;
    fF[11] = lambda_xg_2*x_g2;
    fF[12] = lambda_xg_3*x_g3;
    
    
    precision sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 = sin_phi_2/R_N_plus_h_dot_cos_phi_2;
    precision diff1_dot_cos_phi = (diff1)*cos_phi;
    precision sin_phi_frac_R_N_plus_h_dot_cos_phi = sin_phi/R_N_plus_h_dot_cos_phi;
    precision sin_phi_frac_R_N_plus_h_2_dot_cos_phi = sin_phi/R_N_plus_h_2_dot_cos_phi;
    precision v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi = v_e*sin_phi_frac_R_N_plus_h_dot_cos_phi;
    precision sum1_dot_sin_phi = sum1*sin_phi;
    precision vn_frac_R_M_plus_h_2 = vn_frac_R_M_plus_h/R_M_plus_h;
    // row 1
    fF[13] = 0;
    fF[14] = 0;
    fF[15] = inv_PI_180 *(-v_n/R_M_plus_h_2);
    fF[16] = inv_PI_180 /R_M_plus_h;
    fF[17] = 0;
    fF[18] = 0;
    fF[19] = 0;
    fF[20] = 0;
    fF[21] = 0;
    fF[22] = 0;
    fF[23] = 0;
    fF[24] = 0;
    fF[25] = 0;
    // row 2
    fF[26] = inv_PI_180 *v_e*sin_phi/R_N_plus_h_dot_cos_phi_2;
    fF[27] = 0;
    fF[28] = -inv_PI_180 *v_e/R_N_plus_h_2_dot_cos_phi;
    fF[29] = 0;
    fF[30] = inv_PI_180 /R_N_plus_h_dot_cos_phi;
    fF[31] = 0;
    fF[32] = 0;
    fF[33] = 0;
    fF[34] = 0;
    fF[35]= 0;
    fF[36]= 0;
    fF[37]= 0;
    fF[38]= 0;
    // row 3
    fF[39] = 0;
    fF[40] = 0;
    fF[41]= 0;
    fF[42]= 0;
    fF[43]= 0;
    fF[44] = -1;
    fF[45]= 0;
    fF[46]= 0;
    fF[47]= 0;
    fF[48]= 0;
    fF[49]= 0;
    fF[50]= 0;
    fF[51]= 0;
    // row 4
    fF[52] = -v_e_2*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + v_e*diff1_dot_cos_phi;
    fF[53] = 0;
    fF[54] = -v_d*vn_frac_R_M_plus_h_2 + v_e_2*sin_phi_frac_R_N_plus_h_2_dot_cos_phi;
    fF[55] = v_d/R_M_plus_h;
    fF[56] = -v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + diff1_dot_sin_phi;
    fF[57] = vn_frac_R_M_plus_h;
    fF[58] = 0;
    fF[59] =0;
    fF[60] =0;
    fF[61] =0;
    fF[62] =0;
    fF[63] =0;
    fF[64] =0;
    // row 5
    fF[65] =v_d*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - v_d*sum1_dot_sin_phi + v_e*v_n*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - v_n*diff1_dot_cos_phi;
    fF[66] = 0;
    fF[67] = -v_d*v_e/R_N_plus_h_2 - v_e*v_n*sin_phi_frac_R_N_plus_h_2_dot_cos_phi;
    fF[68] = -diff1_dot_sin_phi;
    fF[69] = v_d/R_N_plus_h + v_n*sin_phi_frac_R_N_plus_h_dot_cos_phi;
    fF[70] = sum1*cos_phi;
    fF[71] =0;
    fF[72] =0;
    fF[73] =0;
    fF[74] =0;
    fF[75] =0;
    fF[76] =0;
    fF[77] =0;
    // row 6
    fF[78] =-v_e_2*sin_phi_frac_R_N_plus_h_dot_cos_phi + v_e*sum1_dot_sin_phi;
    fF[79] = 0;
    fF[80] = v_e_2/R_N_plus_h_2 + v_n*vn_frac_R_M_plus_h_2;
    fF[81] = -2*vn_frac_R_M_plus_h;
    fF[82] =-v_e/R_N_plus_h - sum1*cos_phi;
    fF[83] = 0;
    fF[84] =0;
    fF[85] =0;
    fF[86] =0;
    fF[87] =0;
    fF[88] =0;
    fF[89] =0;
    fF[90] =0;
    // row 7
    fF[91] = -0.5*q1*(q0*(q0*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q0*sum1_dot_sin_phi - q2*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q2*diff1_dot_cos_phi) - q1*(-q1*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q1*sum1_dot_sin_phi + q3*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q3*diff1_dot_cos_phi) + q2*(-q0*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q0*diff1_dot_cos_phi - q2*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q2*sum1_dot_sin_phi) - q3*(q1*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q1*diff1_dot_cos_phi + q3*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q3*sum1_dot_sin_phi)) - 0.5*q2*(q0*(q1*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q1*diff1_dot_cos_phi + q3*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q3*sum1_dot_sin_phi) - q1*(-q0*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q0*diff1_dot_cos_phi - q2*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q2*sum1_dot_sin_phi) - q2*(-q1*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q1*sum1_dot_sin_phi + q3*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q3*diff1_dot_cos_phi) + q3*(q0*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q0*sum1_dot_sin_phi - q2*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q2*diff1_dot_cos_phi)) - 0.5*q3*(q0*(-q0*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q0*diff1_dot_cos_phi - q2*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q2*sum1_dot_sin_phi) + q1*(q1*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q1*diff1_dot_cos_phi + q3*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q3*sum1_dot_sin_phi) - q2*(q0*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q0*sum1_dot_sin_phi - q2*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q2*diff1_dot_cos_phi) - q3*(-q1*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q1*sum1_dot_sin_phi + q3*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q3*diff1_dot_cos_phi));
    fF[92] = 0;
    fF[93] = -0.5*q1*(q0*(-q0*v_e/R_N_plus_h_2 + q2*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*vn_frac_R_M_plus_h_2) - q1*(q1*v_e/R_N_plus_h_2 - q2*vn_frac_R_M_plus_h_2 - q3*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi) + q2*(q0*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi + q1*vn_frac_R_M_plus_h_2 + q2*v_e/R_N_plus_h_2) - q3*(q0*vn_frac_R_M_plus_h_2 - q1*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*v_e/R_N_plus_h_2)) - 0.5*q2*(q0*(q0*vn_frac_R_M_plus_h_2 - q1*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*v_e/R_N_plus_h_2) - q1*(q0*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi + q1*vn_frac_R_M_plus_h_2 + q2*v_e/R_N_plus_h_2) - q2*(q1*v_e/R_N_plus_h_2 - q2*vn_frac_R_M_plus_h_2 - q3*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi) + q3*(-q0*v_e/R_N_plus_h_2 + q2*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*vn_frac_R_M_plus_h_2)) - 0.5*q3*(q0*(q0*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi + q1*vn_frac_R_M_plus_h_2 + q2*v_e/R_N_plus_h_2) + q1*(q0*vn_frac_R_M_plus_h_2 - q1*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*v_e/R_N_plus_h_2) - q2*(-q0*v_e/R_N_plus_h_2 + q2*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*vn_frac_R_M_plus_h_2) - q3*(q1*v_e/R_N_plus_h_2 - q2*vn_frac_R_M_plus_h_2 - q3*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi));
    fF[94] = -0.5*q1*(2*q0*q3/R_M_plus_h - 2*q1*q2/R_M_plus_h) - 0.5*q2*(-q0*q0/R_M_plus_h + q1*q1/R_M_plus_h - q2*q2/R_M_plus_h + q3*q3/R_M_plus_h) - 0.5*q3*(-2*q0*q1/R_M_plus_h - 2*q2*q3/R_M_plus_h);
    fF[95] = -0.5*q1*(q0*(q0/R_N_plus_h - q2*sin_phi_frac_R_N_plus_h_dot_cos_phi) - q1*(-q1/R_N_plus_h + q3*sin_phi_frac_R_N_plus_h_dot_cos_phi) + q2*(-q0*sin_phi_frac_R_N_plus_h_dot_cos_phi - q2/R_N_plus_h) - q3*(q1*sin_phi_frac_R_N_plus_h_dot_cos_phi + q3/R_N_plus_h)) - 0.5*q2*(q0*(q1*sin_phi_frac_R_N_plus_h_dot_cos_phi + q3/R_N_plus_h) - q1*(-q0*sin_phi_frac_R_N_plus_h_dot_cos_phi - q2/R_N_plus_h) - q2*(-q1/R_N_plus_h + q3*sin_phi_frac_R_N_plus_h_dot_cos_phi) + q3*(q0/R_N_plus_h - q2*sin_phi_frac_R_N_plus_h_dot_cos_phi)) - 0.5*q3*(q0*(-q0*sin_phi_frac_R_N_plus_h_dot_cos_phi - q2/R_N_plus_h) + q1*(q1*sin_phi_frac_R_N_plus_h_dot_cos_phi + q3/R_N_plus_h) - q2*(q0/R_N_plus_h - q2*sin_phi_frac_R_N_plus_h_dot_cos_phi) - q3*(-q1/R_N_plus_h + q3*sin_phi_frac_R_N_plus_h_dot_cos_phi));
    fF[96] = 0;
    fF[97] = -0.5*q1*(2*q0_dot_sum_1_dot_cos_phi + 2*q2_dot_diff_1_dot_sin_phi + 2*q3*vn_frac_R_M_plus_h) - 0.5*q2*(-2*q0*vn_frac_R_M_plus_h - 2*q1*diff1_dot_sin_phi + 2*q3*sum1*cos_phi) - 0.5*q3*(2*q0*diff1_dot_sin_phi - 2*q1*vn_frac_R_M_plus_h - 2*q2*sum1*cos_phi);
    fF[98] = 0.5*omega_ib_b1 - 0.5*q0*sum2 + 0.5*q1*sum6 - 0.5*q1*(2*q1*sum1*cos_phi - 2*q2*vn_frac_R_M_plus_h + 2*q3*diff1_dot_sin_phi) - 0.5*q2*(-2*q0*diff1_dot_sin_phi + 2*q1*vn_frac_R_M_plus_h + 2*q2*sum1*cos_phi) - 0.5*q2*sum7 - 0.5*q3*(-2*q0*vn_frac_R_M_plus_h - 2*q1*diff1_dot_sin_phi + 2*q3*sum1*cos_phi) + 0.5*q3*sum4 - 0.5*x_g1;
    fF[99] = 0.5*omega_ib_b2 - 0.5*q0*sum4 + 0.5*q1*sum7 - 0.5*q1*(2*q0*diff1_dot_sin_phi - 2*q1*vn_frac_R_M_plus_h - 2*q2*sum1*cos_phi) + 0.5*q2*sum6 - 0.5*q2*(2*q1*sum1*cos_phi - 2*q2*vn_frac_R_M_plus_h + 2*q3*diff1_dot_sin_phi) - 0.5*q3*(-2*q0_dot_sum_1_dot_cos_phi - 2*q2_dot_diff_1_dot_sin_phi - 2*q3*vn_frac_R_M_plus_h) - 0.5*q3*sum2 - 0.5*x_g2;
    
    fF[100] = 0.5*omega_ib_b3 - 0.5*q0*sum7 - 0.5*q1*sum4 - 0.5*q1*(2*q0*vn_frac_R_M_plus_h + 2*q1*diff1_dot_sin_phi - 2*q3*sum1*cos_phi) + 0.5*q2*sum2 - 0.5*q2*(2*q0_dot_sum_1_dot_cos_phi + 2*q2_dot_diff_1_dot_sin_phi + 2*q3*vn_frac_R_M_plus_h) + 0.5*q3*sum6 - 0.5*q3*(2*q1*sum1*cos_phi - 2*q2*vn_frac_R_M_plus_h + 2*q3*diff1_dot_sin_phi) - 0.5*x_g3;
    fF[101] =  -0.5*q1;
    fF[102] = -0.5*q2;
    fF[103] = -0.5*q3;
    // row 7
    fF[104] = 0.5*q0*(q0*(q0*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q0*sum1_dot_sin_phi - q2*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q2*diff1_dot_cos_phi) - q1*(-q1*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q1*sum1_dot_sin_phi + q3*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q3*diff1_dot_cos_phi) + q2*(-q0*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q0*diff1_dot_cos_phi - q2*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q2*sum1_dot_sin_phi) - q3*(q1*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q1*diff1_dot_cos_phi + q3*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q3*sum1_dot_sin_phi)) - 0.5*q2*(q0*(-q0*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q0*diff1_dot_cos_phi - q2*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q2*sum1_dot_sin_phi) + q1*(q1*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q1*diff1_dot_cos_phi + q3*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q3*sum1_dot_sin_phi) - q2*(q0*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q0*sum1_dot_sin_phi - q2*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q2*diff1_dot_cos_phi) - q3*(-q1*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q1*sum1_dot_sin_phi + q3*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q3*diff1_dot_cos_phi)) + 0.5*q3*(q0*(q1*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q1*diff1_dot_cos_phi + q3*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q3*sum1_dot_sin_phi) - q1*(-q0*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q0*diff1_dot_cos_phi - q2*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q2*sum1_dot_sin_phi) - q2*(-q1*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q1*sum1_dot_sin_phi + q3*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q3*diff1_dot_cos_phi) + q3*(q0*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q0*sum1_dot_sin_phi - q2*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q2*diff1_dot_cos_phi));
    fF[105] = 0;
    fF[106] = 0.5*q0*(q0*(-q0*v_e/R_N_plus_h_2 + q2*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*vn_frac_R_M_plus_h_2) - q1*(q1*v_e/R_N_plus_h_2 - q2*vn_frac_R_M_plus_h_2 - q3*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi) + q2*(q0*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi + q1*vn_frac_R_M_plus_h_2 + q2*v_e/R_N_plus_h_2) - q3*(q0*vn_frac_R_M_plus_h_2 - q1*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*v_e/R_N_plus_h_2)) - 0.5*q2*(q0*(q0*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi + q1*vn_frac_R_M_plus_h_2 + q2*v_e/R_N_plus_h_2) + q1*(q0*vn_frac_R_M_plus_h_2 - q1*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*v_e/R_N_plus_h_2) - q2*(-q0*v_e/R_N_plus_h_2 + q2*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*vn_frac_R_M_plus_h_2) - q3*(q1*v_e/R_N_plus_h_2 - q2*vn_frac_R_M_plus_h_2 - q3*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi)) + 0.5*q3*(q0*(q0*vn_frac_R_M_plus_h_2 - q1*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*v_e/R_N_plus_h_2) - q1*(q0*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi + q1*vn_frac_R_M_plus_h_2 + q2*v_e/R_N_plus_h_2) - q2*(q1*v_e/R_N_plus_h_2 - q2*vn_frac_R_M_plus_h_2 - q3*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi) + q3*(-q0*v_e/R_N_plus_h_2 + q2*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*vn_frac_R_M_plus_h_2));
    fF[107] = 0.5*q0*(2*q0*q3/R_M_plus_h - 2*q1*q2/R_M_plus_h) - 0.5*q2*(-2*q0*q1/R_M_plus_h - 2*q2*q3/R_M_plus_h) + 0.5*q3*(-q0*q0/R_M_plus_h + q1*q1/R_M_plus_h - q2*q2/R_M_plus_h + q3*q3/R_M_plus_h);
    fF[108] = 0.5*q0*(q0*(q0/R_N_plus_h - q2*sin_phi_frac_R_N_plus_h_dot_cos_phi) - q1*(-q1/R_N_plus_h + q3*sin_phi_frac_R_N_plus_h_dot_cos_phi) + q2*(-q0*sin_phi_frac_R_N_plus_h_dot_cos_phi - q2/R_N_plus_h) - q3*(q1*sin_phi_frac_R_N_plus_h_dot_cos_phi + q3/R_N_plus_h)) - 0.5*q2*(q0*(-q0*sin_phi_frac_R_N_plus_h_dot_cos_phi - q2/R_N_plus_h) + q1*(q1*sin_phi_frac_R_N_plus_h_dot_cos_phi + q3/R_N_plus_h) - q2*(q0/R_N_plus_h - q2*sin_phi_frac_R_N_plus_h_dot_cos_phi) - q3*(-q1/R_N_plus_h + q3*sin_phi_frac_R_N_plus_h_dot_cos_phi)) + 0.5*q3*(q0*(q1*sin_phi_frac_R_N_plus_h_dot_cos_phi + q3/R_N_plus_h) - q1*(-q0*sin_phi_frac_R_N_plus_h_dot_cos_phi - q2/R_N_plus_h) - q2*(-q1/R_N_plus_h + q3*sin_phi_frac_R_N_plus_h_dot_cos_phi) + q3*(q0/R_N_plus_h - q2*sin_phi_frac_R_N_plus_h_dot_cos_phi));
    fF[109] = 0;
    fF[110] = -0.5*omega_ib_b1 + 0.5*q0*sum2 + 0.5*q0*(2*q0_dot_sum_1_dot_cos_phi + 2*q2_dot_diff_1_dot_sin_phi + 2*q3*vn_frac_R_M_plus_h) - 0.5*q1*sum6 + 0.5*q2*sum7 - 0.5*q2*(2*q0*diff1_dot_sin_phi - 2*q1*vn_frac_R_M_plus_h - 2*q2*sum1*cos_phi) + 0.5*q3*(-2*q0*vn_frac_R_M_plus_h - 2*q1*diff1_dot_sin_phi + 2*q3*sum1*cos_phi) - 0.5*q3*sum4 + 0.5*x_g1;
    fF[111] = 0.5*q0*(2*q1*sum1*cos_phi - 2*q2*vn_frac_R_M_plus_h + 2*q3*diff1_dot_sin_phi) - 0.5*q2*(-2*q0*vn_frac_R_M_plus_h - 2*q1*diff1_dot_sin_phi + 2*q3*sum1*cos_phi) + 0.5*q3*(-2*q0*diff1_dot_sin_phi + 2*q1*vn_frac_R_M_plus_h + 2*q2*sum1*cos_phi);
    fF[112] = 0.5*omega_ib_b3 - 0.5*q0*sum7 + 0.5*q0*(2*q0*diff1_dot_sin_phi - 2*q1*vn_frac_R_M_plus_h - 2*q2*sum1*cos_phi) - 0.5*q1*sum4 - 0.5*q2*(-2*q0_dot_sum_1_dot_cos_phi - 2*q2_dot_diff_1_dot_sin_phi - 2*q3*vn_frac_R_M_plus_h) + 0.5*q2*sum2 + 0.5*q3*sum6 + 0.5*q3*(2*q1*sum1*cos_phi - 2*q2*vn_frac_R_M_plus_h + 2*q3*diff1_dot_sin_phi) - 0.5*x_g3;
    fF[113] = -0.5*omega_ib_b2 + 0.5*q0*sum4 + 0.5*q0*(2*q0*vn_frac_R_M_plus_h + 2*q1*diff1_dot_sin_phi - 2*q3*sum1*cos_phi) - 0.5*q1*sum7 - 0.5*q2*sum6 - 0.5*q2*(2*q1*sum1*cos_phi - 2*q2*vn_frac_R_M_plus_h + 2*q3*diff1_dot_sin_phi) + 0.5*q3*sum2 + 0.5*q3*(2*q0_dot_sum_1_dot_cos_phi + 2*q2_dot_diff_1_dot_sin_phi + 2*q3*vn_frac_R_M_plus_h) + 0.5*x_g2;
    fF[114] = 0.5*q0;
    fF[115] = 0.5*q3;
    fF[116] =  -0.5*q2;
    // row 8
    fF[117] = 0.5*q0*(q0*(q1*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q1*diff1_dot_cos_phi + q3*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q3*sum1_dot_sin_phi) - q1*(-q0*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q0*diff1_dot_cos_phi - q2*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q2*sum1_dot_sin_phi) - q2*(-q1*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q1*sum1_dot_sin_phi + q3*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q3*diff1_dot_cos_phi) + q3*(q0*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q0*sum1_dot_sin_phi - q2*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q2*diff1_dot_cos_phi)) + 0.5*q1*(q0*(-q0*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q0*diff1_dot_cos_phi - q2*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q2*sum1_dot_sin_phi) + q1*(q1*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q1*diff1_dot_cos_phi + q3*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q3*sum1_dot_sin_phi) - q2*(q0*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q0*sum1_dot_sin_phi - q2*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q2*diff1_dot_cos_phi) - q3*(-q1*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q1*sum1_dot_sin_phi + q3*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q3*diff1_dot_cos_phi)) - 0.5*q3*(q0*(q0*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q0*sum1_dot_sin_phi - q2*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q2*diff1_dot_cos_phi) - q1*(-q1*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q1*sum1_dot_sin_phi + q3*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q3*diff1_dot_cos_phi) + q2*(-q0*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q0*diff1_dot_cos_phi - q2*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q2*sum1_dot_sin_phi) - q3*(q1*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q1*diff1_dot_cos_phi + q3*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q3*sum1_dot_sin_phi));
    fF[118] = 0;
    fF[119] = 0.5*q0*(q0*(q0*vn_frac_R_M_plus_h_2 - q1*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*v_e/R_N_plus_h_2) - q1*(q0*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi + q1*vn_frac_R_M_plus_h_2 + q2*v_e/R_N_plus_h_2) - q2*(q1*v_e/R_N_plus_h_2 - q2*vn_frac_R_M_plus_h_2 - q3*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi) + q3*(-q0*v_e/R_N_plus_h_2 + q2*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*vn_frac_R_M_plus_h_2)) + 0.5*q1*(q0*(q0*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi + q1*vn_frac_R_M_plus_h_2 + q2*v_e/R_N_plus_h_2) + q1*(q0*vn_frac_R_M_plus_h_2 - q1*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*v_e/R_N_plus_h_2) - q2*(-q0*v_e/R_N_plus_h_2 + q2*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*vn_frac_R_M_plus_h_2) - q3*(q1*v_e/R_N_plus_h_2 - q2*vn_frac_R_M_plus_h_2 - q3*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi)) - 0.5*q3*(q0*(-q0*v_e/R_N_plus_h_2 + q2*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*vn_frac_R_M_plus_h_2) - q1*(q1*v_e/R_N_plus_h_2 - q2*vn_frac_R_M_plus_h_2 - q3*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi) + q2*(q0*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi + q1*vn_frac_R_M_plus_h_2 + q2*v_e/R_N_plus_h_2) - q3*(q0*vn_frac_R_M_plus_h_2 - q1*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*v_e/R_N_plus_h_2));
    fF[120] = 0.5*q0*(-q0*q0/R_M_plus_h + q1*q1/R_M_plus_h - q2*q2/R_M_plus_h + q3*q3/R_M_plus_h) + 0.5*q1*(-2*q0*q1/R_M_plus_h - 2*q2*q3/R_M_plus_h) - 0.5*q3*(2*q0*q3/R_M_plus_h - 2*q1*q2/R_M_plus_h);
    fF[121] = 0.5*q0*(q0*(q1*sin_phi_frac_R_N_plus_h_dot_cos_phi + q3/R_N_plus_h) - q1*(-q0*sin_phi_frac_R_N_plus_h_dot_cos_phi - q2/R_N_plus_h) - q2*(-q1/R_N_plus_h + q3*sin_phi_frac_R_N_plus_h_dot_cos_phi) + q3*(q0/R_N_plus_h - q2*sin_phi_frac_R_N_plus_h_dot_cos_phi)) + 0.5*q1*(q0*(-q0*sin_phi_frac_R_N_plus_h_dot_cos_phi - q2/R_N_plus_h) + q1*(q1*sin_phi_frac_R_N_plus_h_dot_cos_phi + q3/R_N_plus_h) - q2*(q0/R_N_plus_h - q2*sin_phi_frac_R_N_plus_h_dot_cos_phi) - q3*(-q1/R_N_plus_h + q3*sin_phi_frac_R_N_plus_h_dot_cos_phi)) - 0.5*q3*(q0*(q0/R_N_plus_h - q2*sin_phi_frac_R_N_plus_h_dot_cos_phi) - q1*(-q1/R_N_plus_h + q3*sin_phi_frac_R_N_plus_h_dot_cos_phi) + q2*(-q0*sin_phi_frac_R_N_plus_h_dot_cos_phi - q2/R_N_plus_h) - q3*(q1*sin_phi_frac_R_N_plus_h_dot_cos_phi + q3/R_N_plus_h));
    fF[122] = 0;
    fF[123] = -0.5*omega_ib_b2 + 0.5*q0*(-2*q0*vn_frac_R_M_plus_h - 2*q1*diff1_dot_sin_phi + 2*q3*sum1*cos_phi) + 0.5*q0*sum4 - 0.5*q1*sum7 + 0.5*q1*(2*q0*diff1_dot_sin_phi - 2*q1*vn_frac_R_M_plus_h - 2*q2*sum1*cos_phi) - 0.5*q2*sum6 + 0.5*q3*sum2 - 0.5*q3*(2*q0_dot_sum_1_dot_cos_phi + 2*q2_dot_diff_1_dot_sin_phi + 2*q3*vn_frac_R_M_plus_h) + 0.5*x_g2;
    fF[124] = -0.5*omega_ib_b3 + 0.5*q0*(-2*q0*diff1_dot_sin_phi + 2*q1*vn_frac_R_M_plus_h + 2*q2*sum1*cos_phi) + 0.5*q0*sum7 + 0.5*q1*(-2*q0*vn_frac_R_M_plus_h - 2*q1*diff1_dot_sin_phi + 2*q3*sum1*cos_phi) + 0.5*q1*sum4 - 0.5*q2*sum2 - 0.5*q3*sum6 - 0.5*q3*(2*q1*sum1*cos_phi - 2*q2*vn_frac_R_M_plus_h + 2*q3*diff1_dot_sin_phi) + 0.5*x_g3;
    fF[125] = 0.5*q0*(2*q1*sum1*cos_phi - 2*q2*vn_frac_R_M_plus_h + 2*q3*diff1_dot_sin_phi) + 0.5*q1*(-2*q0_dot_sum_1_dot_cos_phi - 2*q2_dot_diff_1_dot_sin_phi - 2*q3*vn_frac_R_M_plus_h) - 0.5*q3*(2*q0*diff1_dot_sin_phi - 2*q1*vn_frac_R_M_plus_h - 2*q2*sum1*cos_phi);
    fF[126] = 0.5*omega_ib_b1 - 0.5*q0*sum2 + 0.5*q0*(2*q0_dot_sum_1_dot_cos_phi + 2*q2_dot_diff_1_dot_sin_phi + 2*q3*vn_frac_R_M_plus_h) + 0.5*q1*sum6 + 0.5*q1*(2*q1*sum1*cos_phi - 2*q2*vn_frac_R_M_plus_h + 2*q3*diff1_dot_sin_phi) - 0.5*q2*sum7 + 0.5*q3*sum4 - 0.5*q3*(2*q0*vn_frac_R_M_plus_h + 2*q1*diff1_dot_sin_phi - 2*q3*sum1*cos_phi) - 0.5*x_g1;
    fF[127] =  -0.5*q3;
    fF[128] = 0.5*q0;
    fF[129] = 0.5*q1;
    // row 9
    fF[130] = 0.5*q0*(q0*(-q0*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q0*diff1_dot_cos_phi - q2*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q2*sum1_dot_sin_phi) + q1*(q1*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q1*diff1_dot_cos_phi + q3*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q3*sum1_dot_sin_phi) - q2*(q0*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q0*sum1_dot_sin_phi - q2*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q2*diff1_dot_cos_phi) - q3*(-q1*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q1*sum1_dot_sin_phi + q3*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q3*diff1_dot_cos_phi)) - 0.5*q1*(q0*(q1*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q1*diff1_dot_cos_phi + q3*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q3*sum1_dot_sin_phi) - q1*(-q0*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q0*diff1_dot_cos_phi - q2*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q2*sum1_dot_sin_phi) - q2*(-q1*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q1*sum1_dot_sin_phi + q3*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q3*diff1_dot_cos_phi) + q3*(q0*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q0*sum1_dot_sin_phi - q2*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q2*diff1_dot_cos_phi)) + 0.5*q2*(q0*(q0*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q0*sum1_dot_sin_phi - q2*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q2*diff1_dot_cos_phi) - q1*(-q1*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q1*sum1_dot_sin_phi + q3*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q3*diff1_dot_cos_phi) + q2*(-q0*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 + q0*diff1_dot_cos_phi - q2*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi + q2*sum1_dot_sin_phi) - q3*(q1*v_e*sin_phi_2_frac_R_N_plus_h_dot_cos_phi_2 - q1*diff1_dot_cos_phi + q3*v_e_prod_sin_phi_frac_R_N_plus_h_dot_cos_phi - q3*sum1_dot_sin_phi));
    fF[131] = 0;
    fF[132] = 0.5*q0*(q0*(q0*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi + q1*vn_frac_R_M_plus_h_2 + q2*v_e/R_N_plus_h_2) + q1*(q0*vn_frac_R_M_plus_h_2 - q1*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*v_e/R_N_plus_h_2) - q2*(-q0*v_e/R_N_plus_h_2 + q2*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*vn_frac_R_M_plus_h_2) - q3*(q1*v_e/R_N_plus_h_2 - q2*vn_frac_R_M_plus_h_2 - q3*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi)) - 0.5*q1*(q0*(q0*vn_frac_R_M_plus_h_2 - q1*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*v_e/R_N_plus_h_2) - q1*(q0*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi + q1*vn_frac_R_M_plus_h_2 + q2*v_e/R_N_plus_h_2) - q2*(q1*v_e/R_N_plus_h_2 - q2*vn_frac_R_M_plus_h_2 - q3*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi) + q3*(-q0*v_e/R_N_plus_h_2 + q2*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*vn_frac_R_M_plus_h_2)) + 0.5*q2*(q0*(-q0*v_e/R_N_plus_h_2 + q2*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*vn_frac_R_M_plus_h_2) - q1*(q1*v_e/R_N_plus_h_2 - q2*vn_frac_R_M_plus_h_2 - q3*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi) + q2*(q0*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi + q1*vn_frac_R_M_plus_h_2 + q2*v_e/R_N_plus_h_2) - q3*(q0*vn_frac_R_M_plus_h_2 - q1*v_e*sin_phi_frac_R_N_plus_h_2_dot_cos_phi - q3*v_e/R_N_plus_h_2));
    fF[133] = 0.5*q0*(-2*q0*q1/R_M_plus_h - 2*q2*q3/R_M_plus_h) - 0.5*q1*(-q0*q0/R_M_plus_h + q1*q1/R_M_plus_h - q2*q2/R_M_plus_h + q3*q3/R_M_plus_h) + 0.5*q2*(2*q0*q3/R_M_plus_h - 2*q1*q2/R_M_plus_h);
    fF[134] = 0.5*q0*(q0*(-q0*sin_phi_frac_R_N_plus_h_dot_cos_phi - q2/R_N_plus_h) + q1*(q1*sin_phi_frac_R_N_plus_h_dot_cos_phi + q3/R_N_plus_h) - q2*(q0/R_N_plus_h - q2*sin_phi_frac_R_N_plus_h_dot_cos_phi) - q3*(-q1/R_N_plus_h + q3*sin_phi_frac_R_N_plus_h_dot_cos_phi)) - 0.5*q1*(q0*(q1*sin_phi_frac_R_N_plus_h_dot_cos_phi + q3/R_N_plus_h) - q1*(-q0*sin_phi_frac_R_N_plus_h_dot_cos_phi - q2/R_N_plus_h) - q2*(-q1/R_N_plus_h + q3*sin_phi_frac_R_N_plus_h_dot_cos_phi) + q3*(q0/R_N_plus_h - q2*sin_phi_frac_R_N_plus_h_dot_cos_phi)) + 0.5*q2*(q0*(q0/R_N_plus_h - q2*sin_phi_frac_R_N_plus_h_dot_cos_phi) - q1*(-q1/R_N_plus_h + q3*sin_phi_frac_R_N_plus_h_dot_cos_phi) + q2*(-q0*sin_phi_frac_R_N_plus_h_dot_cos_phi - q2/R_N_plus_h) - q3*(q1*sin_phi_frac_R_N_plus_h_dot_cos_phi + q3/R_N_plus_h));
    fF[135] = 0;
    fF[136] = -0.5*omega_ib_b3 + 0.5*q0*sum7 + 0.5*q0*(2*q0*diff1_dot_sin_phi - 2*q1*vn_frac_R_M_plus_h - 2*q2*sum1*cos_phi) - 0.5*q1*(-2*q0*vn_frac_R_M_plus_h - 2*q1*diff1_dot_sin_phi + 2*q3*sum1*cos_phi) + 0.5*q1*sum4 - 0.5*q2*sum2 + 0.5*q2*(2*q0_dot_sum_1_dot_cos_phi + 2*q2_dot_diff_1_dot_sin_phi + 2*q3*vn_frac_R_M_plus_h) - 0.5*q3*sum6 + 0.5*x_g3;
    fF[137] = 0.5*omega_ib_b2 + 0.5*q0*(-2*q0*vn_frac_R_M_plus_h - 2*q1*diff1_dot_sin_phi + 2*q3*sum1*cos_phi) - 0.5*q0*sum4 - 0.5*q1*(-2*q0*diff1_dot_sin_phi + 2*q1*vn_frac_R_M_plus_h + 2*q2*sum1*cos_phi) + 0.5*q1*sum7 + 0.5*q2*sum6 + 0.5*q2*(2*q1*sum1*cos_phi - 2*q2*vn_frac_R_M_plus_h + 2*q3*diff1_dot_sin_phi) - 0.5*q3*sum2 - 0.5*x_g2;
    fF[138] = -0.5*omega_ib_b1 + 0.5*q0*(-2*q0_dot_sum_1_dot_cos_phi - 2*q2_dot_diff_1_dot_sin_phi - 2*q3*vn_frac_R_M_plus_h) + 0.5*q0*sum2 - 0.5*q1*sum6 - 0.5*q1*(2*q1*sum1*cos_phi - 2*q2*vn_frac_R_M_plus_h + 2*q3*diff1_dot_sin_phi) + 0.5*q2*sum7 + 0.5*q2*(2*q0*diff1_dot_sin_phi - 2*q1*vn_frac_R_M_plus_h - 2*q2*sum1*cos_phi) - 0.5*q3*sum4 + 0.5*x_g1;
    fF[139] = 0.5*q0*(2*q1*sum1*cos_phi - 2*q2*vn_frac_R_M_plus_h + 2*q3*diff1_dot_sin_phi) - 0.5*q1*(2*q0_dot_sum_1_dot_cos_phi + 2*q2_dot_diff_1_dot_sin_phi + 2*q3*vn_frac_R_M_plus_h) + 0.5*q2*(2*q0*vn_frac_R_M_plus_h + 2*q1*diff1_dot_sin_phi - 2*q3*sum1*cos_phi);
    fF[140] = 0.5*q2;
    fF[141] =  -0.5*q1;
    fF[142] = 0.5*q0;
    // row 10
    fF[143] = 0;
    fF[144] = 0;
    fF[145] = 0;
    fF[146] = 0;
    fF[147] = 0;
    fF[148] = 0;
    fF[149] = 0;
    fF[150] = 0;
    fF[151] = 0;
    fF[152] = 0;
    fF[153] = lambda_xg_1;
    fF[154] = 0;
    fF[155] = 0;
    // row 11
    fF[156] = 0;
    fF[157] = 0;
    fF[158] = 0;
    fF[159] = 0;
    fF[160] = 0;
    fF[161] = 0;
    fF[162] = 0;
    fF[163] = 0;
    fF[164] = 0;
    fF[165] = 0;
    fF[166] = 0;
    fF[167] = lambda_xg_2;
    fF[168] = 0;
    // row 12
    fF[169] = 0;
    fF[170] = 0;
    fF[171] = 0;
    fF[172] = 0;
    fF[173] = 0;
    fF[174] = 0;
    fF[175] = 0;
    fF[176] = 0;
    fF[177] = 0;
    fF[178] = 0;
    fF[179] = 0;
    fF[180] = 0;
    fF[181] = lambda_xg_3;
	
	for (int i = 0; i < 182; i++)
	{precision tmp = fF[i];
		if (tmp!=0)
			fF[i] = tmp*dt;
	}
	
    for (int i = 0; i < 13; i ++)
    {
        size_t idx = i*14+13;
        fF[idx] = fF[idx] + 1;
    }
	
	free(q);
	free(q_conj);
	free(f_nav);
	free(q_acc);
	free(f_nav_1);
	free(lambda_g);
   // return fF;
}

//precision* compute_INS_h_and_H(precision *x,precision *measurement_gyro,precision *measurement_acc,precision *measurement_mag, precision *lambda_g)
void compute_INS_h_and_H(EKF *ekf)
{
	precision *hH_INS = ekf->h;
	precision *x = ekf->x;
	precision *z = ekf->u + 6;
    //precision* hH_INS = (precision*) malloc((sizeof(precision)*168));
    
    // Here z represents u(t), that is, [w_x,w_y,w_z],[a_x,a_y,a_z],[magx,magy,magz] in this order
    //precision *measurement_gyro = &z[0];
    precision *measurement_acc  = &z[3];
    precision *measurement_mag  = &z[6];
    
    precision phi = x[0];
    precision lam = x[1];
    precision h   = x[2];
    precision v_n = x[3];
    precision v_e = x[4];
    precision v_d = x[5];
    
    precision q0 = x[6];
    precision q1 = x[7];
    precision q2 = x[8];
    precision q3 = x[9];
    
    quaternion_t *q = quaternion_new_set(q0,q1,q2,q3);
    quaternion_t *q_conj = quaternion_new();
    quaternConj(q_conj,q);
    quaternion_t *f_nav = quaternion_new();
	
    
    precision normMag = measurement_mag[0]*measurement_mag[0] + measurement_mag[1]*measurement_mag[1] + measurement_mag[2]*measurement_mag[2];
    precision m_1 = measurement_mag[0]*invSqrt(normMag);
    precision m_2 = measurement_mag[1]*invSqrt(normMag);
    precision m_3 = measurement_mag[2]*invSqrt(normMag);
    
    precision a_1 = measurement_acc[0];
    precision a_2 = measurement_acc[1];
    precision a_3 = measurement_acc[2];
    
    quaternion_t *q_acc = quaternion_new_set(0,a_1,a_2,a_3);
    
    quaternProd(f_nav,q_conj,q_acc);
    quaternion_t *f_nav_1 = quaternion_new();
    quaternProd(f_nav_1,f_nav,q);
	

    
    
    hH_INS[0]  = phi;
    hH_INS[1]  = lam;
    hH_INS[2]  = h;
    hH_INS[3]  = v_n;
    hH_INS[4]  = v_e;
    hH_INS[5]  = v_d;
    hH_INS[6]  = q0*(a_1*q0 + a_2*q3 - a_3*q2) + q1*(a_1*q1 + a_2*q2 + a_3*q3) - q2*(a_1*q2 - a_2*q1 + a_3*q0) + q3*(-a_1*q3 + a_2*q0 + a_3*q1);
    hH_INS[7]  = q0*(-a_1*q3 + a_2*q0 + a_3*q1) + q1*(a_1*q2 - a_2*q1 + a_3*q0) + q2*(a_1*q1 + a_2*q2 + a_3*q3) - q3*(a_1*q0 + a_2*q3 - a_3*q2);
    hH_INS[8]  = q0*(a_1*q2 - a_2*q1 + a_3*q0) - q1*(-a_1*q3 + a_2*q0 + a_3*q1) + q2*(a_1*q0 + a_2*q3 - a_3*q2) + q3*(a_1*q1 + a_2*q2 + a_3*q3);
    hH_INS[9]  = q0*(m_1*q0 + m_2*q3 - m_3*q2) + q1*(m_1*q1 + m_2*q2 + m_3*q3) - q2*(m_1*q2 - m_2*q1 + m_3*q0) + q3*(-m_1*q3 + m_2*q0 + m_3*q1);
    hH_INS[10] = q0*(-m_1*q3 + m_2*q0 + m_3*q1) + q1*(m_1*q2 - m_2*q1 + m_3*q0) + q2*(m_1*q1 + m_2*q2 + m_3*q3) - q3*(m_1*q0 + m_2*q3 - m_3*q2);
    hH_INS[11] = q0*(m_1*q2 - m_2*q1 + m_3*q0) - q1*(-m_1*q3 + m_2*q0 + m_3*q1) + q2*(m_1*q0 + m_2*q3 - m_3*q2) + q3*(m_1*q1 + m_2*q2 + m_3*q3);
    
    // H is 13x12 = 156 elements
    
    precision a_1_tim_2 = 2*a_1;
    precision a_2_tim_2 = 2*a_2;
    precision a_3_tim_2 = 2*a_3;
    
    precision m_1_tim_2 = 2*m_1;
    precision m_2_tim_2 = 2*m_2;
    precision m_3_tim_2 = 2*m_3;
	
	if (ekf->H[0] != 1.0){
		hH_INS[12] = 1;
		hH_INS[13] = 0;
		hH_INS[14] = 0;
		hH_INS[15] = 0;
		hH_INS[16] = 0;
		hH_INS[17] = 0;
		hH_INS[18] = 0;
		hH_INS[19] = 0;
		hH_INS[20] = 0;
		hH_INS[21] = 0;
		hH_INS[22] = 0;
		hH_INS[23] = 0;
		hH_INS[24] = 0;
		
		
		// row 2
		hH_INS[25] = 0;
		hH_INS[26] = 1;
		hH_INS[27] = 0;
		hH_INS[28] = 0;
		hH_INS[29] = 0;
		hH_INS[30] = 0;
		hH_INS[31] = 0;
		hH_INS[32] = 0;
		hH_INS[33] = 0;
		hH_INS[34] = 0;
		hH_INS[35] = 0;
		hH_INS[36] = 0;
		hH_INS[37] = 0;
		// row 3
		
		hH_INS[38] = 0;
		hH_INS[39] = 0;
		hH_INS[40] = 1;
		hH_INS[41] = 0;
		hH_INS[42] = 0;
		hH_INS[43] = 0;
		hH_INS[44] = 0;
		hH_INS[45] = 0;
		hH_INS[46] = 0;
		hH_INS[47] = 0;
		hH_INS[48] = 0;
		hH_INS[49] = 0;
		hH_INS[50] = 0;
		
		// row 4
		hH_INS[51] = 0;
		hH_INS[52] = 0;
		hH_INS[53] = 0;
		hH_INS[54] = 1;
		hH_INS[55] = 0;
		hH_INS[56] = 0;
		hH_INS[57] = 0;
		hH_INS[58] = 0;
		hH_INS[59] = 0;
		hH_INS[60] = 0;
		hH_INS[61] = 0;
		hH_INS[62] = 0;
		hH_INS[63] = 0;
		
		// row 5
		hH_INS[64] = 0;
		hH_INS[65] = 0;
		hH_INS[66] = 0;
		hH_INS[67] = 0;
		hH_INS[68] = 1;
		hH_INS[69] = 0;
		hH_INS[70] = 0;
		hH_INS[71] = 0;
		hH_INS[72] = 0;
		hH_INS[73] = 0;
		hH_INS[74] = 0;
		hH_INS[75] = 0;
		hH_INS[76] = 0;
		
		//row 6
		hH_INS[77] = 0;
		hH_INS[78] = 0;
		hH_INS[79] = 0;
		hH_INS[80] = 0;
		hH_INS[81] = 0;
		hH_INS[82] = 1;
		hH_INS[83] = 0;
		hH_INS[84] = 0;
		hH_INS[85] = 0;
		hH_INS[86] = 0;
		hH_INS[87] = 0;
		hH_INS[88] = 0;
		hH_INS[89] = 0;
	}
		
    // row 7
    hH_INS[90]  = 0;
    hH_INS[91]  = 0;
    hH_INS[92]  = 0;
    hH_INS[93]  = 0;
    hH_INS[94]  = 0;
    hH_INS[95]  = 0;
    hH_INS[96]  = a_1_tim_2*q0 + a_2_tim_2*q3 - a_3_tim_2*q2;
    hH_INS[97]  = a_1_tim_2*q1 + a_2_tim_2*q2 + a_3_tim_2*q3;
    hH_INS[98]  = -a_1_tim_2*q2 + a_2_tim_2*q1 - a_3_tim_2*q0;
    hH_INS[99]  = -a_1_tim_2*q3 + a_2_tim_2*q0 + a_3_tim_2*q1;
    hH_INS[100] = 0;
    hH_INS[101] = 0;
    hH_INS[102] = 0;
    
    // row 8
    hH_INS[103] = 0;
    hH_INS[104] = 0;
    hH_INS[105] = 0;
    hH_INS[106] = 0;
    hH_INS[107] = 0;
    hH_INS[108] = 0;
    hH_INS[109] = -a_1_tim_2*q3 + a_2_tim_2*q0 + a_3_tim_2*q1;
    hH_INS[110] = a_1_tim_2*q2 - a_2_tim_2*q1 + a_3_tim_2*q0;
    hH_INS[111] = a_1_tim_2*q1 + a_2_tim_2*q2 + a_3_tim_2*q3;
    hH_INS[112] = -a_1_tim_2*q0 - a_2_tim_2*q3 + a_3_tim_2*q2;
    hH_INS[113] = 0;
    hH_INS[114] = 0;
    hH_INS[115] = 0;
    
    // row 9
    hH_INS[116] = 0;
    hH_INS[117] = 0;
    hH_INS[118] = 0;
    hH_INS[119] = 0;
    hH_INS[120] = 0;
    hH_INS[121] = 0;
    hH_INS[122] = a_1_tim_2*q2 - a_2_tim_2*q1 + a_3_tim_2*q0;
    hH_INS[123] = a_1_tim_2*q3 - a_2_tim_2*q0 - a_3_tim_2*q1;
    hH_INS[124] = a_1_tim_2*q0 + a_2_tim_2*q3 - a_3_tim_2*q2;
    hH_INS[125] = a_1_tim_2*q1 + a_2_tim_2*q2 + a_3_tim_2*q3;
    hH_INS[126] = 0;
    hH_INS[127] = 0;
    hH_INS[128] = 0;
    
    // row 10
    hH_INS[129] = 0;
    hH_INS[130] = 0;
    hH_INS[131] = 0;
    hH_INS[132] = 0;
    hH_INS[133] = 0;
    hH_INS[134] = 0;
    hH_INS[135] = m_1_tim_2*q0 + m_2_tim_2*q3 - m_3_tim_2*q2;
    hH_INS[136] = m_1_tim_2*q1 + m_2_tim_2*q2 + m_3_tim_2*q3;
    hH_INS[137] = -m_1_tim_2*q2 + m_2_tim_2*q1 - m_3_tim_2*q0;
    hH_INS[138] = -m_1_tim_2*q3 + m_2_tim_2*q0 + m_3_tim_2*q1;
    hH_INS[139] = 0;
    hH_INS[140] = 0;
    hH_INS[141] = 0;
    
    // row 11
    hH_INS[142] = 0;
    hH_INS[143] = 0;
    hH_INS[144] = 0;
    hH_INS[145] = 0;
    hH_INS[146] = 0;
    hH_INS[147] = 0;
    hH_INS[148] = -m_1_tim_2*q3 + m_2_tim_2*q0 + m_3_tim_2*q1;
    hH_INS[149] = m_1_tim_2*q2 - m_2_tim_2*q1 + m_3_tim_2*q0;
    hH_INS[150] = m_1_tim_2*q1 + m_2_tim_2*q2 + m_3_tim_2*q3;
    hH_INS[151] = -m_1_tim_2*q0 - m_2_tim_2*q3 + m_3_tim_2*q2;
    hH_INS[152] = 0;
    hH_INS[153] = 0;
    hH_INS[154] = 0;
    
    // row 12
    hH_INS[155] = 0;
    hH_INS[156] = 0;
    hH_INS[157] = 0;
    hH_INS[158] = 0;
    hH_INS[159] = 0;
    hH_INS[160] = 0;
    hH_INS[161] = m_1_tim_2*q2 - m_2_tim_2*q1 + m_3_tim_2*q0;
    hH_INS[162] = m_1_tim_2*q3 - m_2_tim_2*q0 - m_3_tim_2*q1;
    hH_INS[163] = m_1_tim_2*q0 + m_2_tim_2*q3 - m_3_tim_2*q2;
    hH_INS[164] = m_1_tim_2*q1 + m_2_tim_2*q2 + m_3_tim_2*q3;
    hH_INS[165] = 0;
    hH_INS[166] = 0;
    hH_INS[167] = 0;
    
	free(q);
	free(q_conj);
	free(f_nav);
	free(q_acc);
	free(f_nav_1);
	
 //   return hH_INS;
}

void normalize_quaternion(EKF *ekf){
	precision recipNorm = invSqrt(ekf->x[6] * ekf->x[6] + ekf->x[7] * ekf->x[7] + ekf->x[8] * ekf->x[8] + ekf->x[9] * ekf->x[9]);
	ekf->x[6] *= recipNorm;
	ekf->x[7] *= recipNorm;
	ekf->x[8] *= recipNorm;
	ekf->x[9] *= recipNorm;
}

void update_R(EKF *ekf){
	precision *R = ekf->R;
	int i;
	if (ekf->u[0] > 1000){
		for (i = 0; i < 6; i++){
			R[12*i + i] = R_GPS_BAD;
		}
	} else {
		R[0] = R_LAT;
		R[12 + 1] = R_LON;
		R[24 + 2] = R_H;
		R[36 + 3] = R_V;
		R[48 + 4] = R_V;
		R[60 + 5] = R_V;
	}
	
	precision *gyro = ekf->u + 6;
	precision *acc = ekf->u + 9;
	precision normAcc_2 = acc[0]*acc[0] + acc[1]*acc[1] + acc[2]*acc[2];
	precision normGyro_2 = gyro[0]*gyro[0] + gyro[1]*gyro[1] + gyro[2]*gyro[2];
	precision normAcc, normGyro, mu;
#if PRECISION == DOUBLE
	normAcc = sqrt(normAcc_2);
	normGyro = sqrt(normGyro_2);
	mu = fabs(normAcc - ge);
#else
	normAcc = sqrtf(normAcc_2);
	normGyro = sqrtf(normGyro_2);
	mu = fabsf(normAcc - ge);
#endif
	
	if ((mu < BETA_MU)  && (normGyro < BETA_GYRO)){
		R[78] = R_ACC_BAD;
		R[91] = R_ACC_BAD;
		R[104] = R_ACC_BAD;
	} else {
		R[78] = R_ACC_X;
		R[91] = R_ACC_Y;
		R[104] = R_ACC_Z;
	}
}
