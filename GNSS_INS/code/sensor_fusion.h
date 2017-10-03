//
//  init_sensor_fusion.h
//  GNSS_INS

#ifndef sensor_fusion_h
#define sensor_fusion_h

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define KNOTS_TO_METERS  0.51444444444444444444444444444444444
#define INIT_STOP_COUNTER 99


typedef struct GNSS_INS{
    
    double hdop;
    double time;
    //double timestamp_IMU;
    double timestamp_GPS;
    double timestamp_ID;
    double IMU_Ts;
    bool flag_IMU_pending; //waiting for IMU
    bool flag_GPS_pending; //waiting for GPS
    double mag_variation; // for declination
    double processed_data[15];
    
}GNSS_INS_t;


typedef struct out_data{

    double time;
    double lat,lon;
    double vn, ve,speed;
    double lateral_acc, deacceleration,acceleration;
    double axis[3];
    double angle;
    double q0,q1,q2,q3;
    double gbx,gby,gbz;
    
}out_data;


void init_sensor_fusion();
void handle_new_data_event(double * imu, double * gps);
void free_sensor_fusion();
#endif /* init_sensor_fusion_h */
