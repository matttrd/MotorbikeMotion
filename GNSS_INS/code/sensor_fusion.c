//
//  init_sensor_fusion.c
//  GNSS_INS
//

#include "sensor_fusion.h"
#include <math.h>
#include “../../EKF/EKFstruct.h"
#include "driver_functions.h"

EKF* ekf;
GNSS_INS_t *ins_master_struct;

void init_sensor_fusion(){
    
    precision x0[13] = { 0.0 };
    precision Q[N_STATES_2] = { 0.0 };
    precision R[N_OBS_2] = { 0.0 };
    precision P[N_STATES_2] = { 0.0 };
    
    ins_master_struct = malloc(sizeof(GNSS_INS_t));
    
    ins_master_struct->hdop = 0;
    ins_master_struct->flag_GPS_pending = false;
    ins_master_struct->flag_IMU_pending = false;
    
    ekf = init_ekf(13, 12, 15, &x0[0], &Q[0], &R[0], &P[0]);
};


void handle_new_data_event(double * imu, double * gps){
    
    
    if (imu){
        
        // ENU TO NED
        
        // swap x and y axes (gyro,acc,mag in this order)
        double tmp;
        tmp = imu[1];
        imu[1] = imu[2];
        imu[2] = tmp;
        
        tmp = imu[4];
        imu[4] = imu[5];
        imu[5] = tmp;
        
        tmp = imu[7];
        imu[7] = imu[8];
        imu[8] = tmp;
        
        // - z for gyro and acc
        
        imu[3] = -imu[3];
        imu[6] = -imu[6];
        
        // multiply for g = 9.81
        imu[4]*=9.81;
        imu[5]*=9.81;
        imu[6]*=9.81;
        
        
        // gyro dps -> rps
        double scale_factor = M_PI/180.0;
        imu[1] *= scale_factor;
        imu[2] *= scale_factor;
        imu[3] *= scale_factor;
        
        if (ins_master_struct->flag_GPS_pending == true)
        {
            
            ins_master_struct->flag_GPS_pending = false;
            double timestamp_ID = imu[0];
            double current_Ts = (timestamp_ID - ins_master_struct->timestamp_ID)*ins_master_struct->IMU_Ts;
            ins_master_struct->timestamp_ID = timestamp_ID;
            
            
            for (int i = 0; i < 9; i++)
                ins_master_struct->processed_data[i+6] = imu[i+1];
            
            ins_master_struct->time += current_Ts;
            updateGNSS_INS(&ins_master_struct->processed_data[0], ekf, current_Ts,ins_master_struct);
        }
        else{
            ins_master_struct->processed_data[0] = 1717;
            ins_master_struct->processed_data[1] = 1717;
            
            double timestamp_ID = imu[0];
            double current_Ts = (timestamp_ID - ins_master_struct->timestamp_ID)*ins_master_struct->IMU_Ts;
            ins_master_struct->time += current_Ts;
            
            updateGNSS_INS(&ins_master_struct->processed_data[0], ekf, current_Ts,ins_master_struct);
            
        }
    }
    else // then GPS
        
    {
        preprocess_GPS_data(gps,ins_master_struct);
        ins_master_struct->flag_GPS_pending = true;
    }
}

void free_sensor_fusion(){
    
    free_ekf(ekf);
    free(ins_master_struct);
}

void get_ouput(out_data* out_struct){
    
    /*out_struct->time = ins_master_struct->time;
     out_struct->lat = ekf->x[0];
     out_struct->lon = ekf->x[1];
     out_struct->vn = ekf->x[3];
     out_struct->ve = ekf->x[4];
     out_struct->q0 = ekf->x[6];
     out_struct->q1 = ekf->x[7];
     out_struct->q2 = ekf->x[8];
     out_struct->q3 = ekf->x[9];
     out_struct->gbx = ekf->x[10];
     out_struct->gby = ekf->x[11];
     out_struct->gbz = ekf->x[12];*/
    
    double v_n = ekf->x[3];
    double v_e = ekf->x[4];
    double speed = sqrt(v_n*v_n + v_e*v_e);
    
    quaternion_t *acc = quaternion_new_set(0, ekf->z[0],ekf->z[1],ekf->z[2]);
    quaternion_t *q_star = quaternion_new_set(ekf->x[6], -ekf->x[7],-ekf->x[8],-ekf->x[9]);
    quaternion_t *acc_nav;
    quaternProd(acc_nav, q_star,acc);
    quaternion_t *q = quaternion_new_set(ekf->x[6], ekf->x[7],ekf->x[8],ekf->x[9]);
    quaternProd(acc_nav, acc_nav,q);
    
    // scalar product
    
    double acc_nav_n = acc_nav->q[1];
    double acc_nav_e = acc_nav->q[2];
    
    double ps = v_n*acc_nav_n + v_e*acc_nav_e;
    
    double acc_motion = ps/speed;
    double acc_lat = sqrt(acc_nav_n*acc_nav_n + acc_nav_e*acc_nav_e - acc_motion*acc_motion);
    
    
    // from NED to SceneKit convention
    double ret[4];
    double q0 = ekf->x[6];
    if (0 == q0){
        ret[3] = 1;
        ret[0] = 0;
        ret[1] = 0;
        ret[2] = 0;
    }
    else{
        ret[3] = 2 * acos(q0);
        // WARNING: changed formulas to take into account of the scenekit reference frame
        double q0q0 = q0*q0;
        ret[0] = - ekf->x[8] * invSqrt(1 - q0q0);
        ret[1] =   ekf->x[9] * invSqrt(1 - q0q0);
        ret[2] = - ekf->x[7] * invSqrt(1 - q0q0);
    }
    
}
