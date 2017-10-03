//
//  main.c
//  GNSS_INS


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Accelerate/Accelerate.h>

#include "driver_functions.h"
#include “../EKF/EKF.h"
#include “./ekf_functions/INS_matrix_fun.h"
#include “./init/initialization.h"

#include <time.h>



int FP = 0;

const char *getfield(char* line, int num)
{
	const char* tok;
	for (tok = strtok(line, ",");
		 tok && *tok;
		 tok = strtok(NULL, ",\n"))
	{
		if (!--num)
			return tok;
	}
	return NULL;
}


int main(){

	precision dt = 0.1;
	
	precision x0[13] = { 0.0 };
	precision Q[N_STATES_2] = { 0.0 };
	precision R[N_OBS_2] = { 0.0 };
	precision P[N_STATES_2] = { 0.0 };
	
    //GNSS_INS_t* ins_struct = init_GNSS_INS_struct();
	EKF *ekf = init_ekf(N_STATES, N_OBS, N_INPUTS, &x0[0], &Q[0], &R[0], &P[0]);
	
	ekf->compute_f_F = &compute_f_and_F;
	ekf->compute_h_H = &compute_INS_h_and_H;
	ekf->update_R = &update_R;
	
	ekf->x[6] = 1;
	
	FILE* stream = fopen("data_INS.txt", "r");

	char line[1024];
	double data[15];
	
	FILE *fp;
	
    fp = fopen("output.txt”, "w");
	
    fprintf(fp, "lat, lon, h, v_n, v_e, v_d, q0, q1, q2, q3, xg_1, xg_2, xg_3\n");
	//=========================================== TIMING====================================================
	//clock_t begin, end;
	//double time_spent;
	//begin = clock();
	//======================================================================================================
	while (fgets(line, 1024, stream)){
		if (FP == 1){
			int i;
			for (i = 0; i < 15; i++){
				char* tmp = strdup(line);
				data[i] = atof(getfield(tmp, i + 1));
				free(tmp);
			}
			updateGNSS_INS(&data[0], ekf, dt, ins_struct);
			fprintf(fp, "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", ekf->x[0], ekf->x[1], ekf->x[2], ekf->x[3], ekf->x[4], ekf->x[5], ekf->x[6], ekf->x[7], ekf->x[8], ekf->x[9], ekf->x[10], ekf->x[11], ekf->x[12]);
		} else {
			FP = 1;
		}
	}
	fclose(fp);
	free_ekf(ekf);
	//=========================================== TIMING====================================================
	//end = clock();
	//time_spent = (double)(end - begin)/ CLOCKS_PER_SEC;
	//printf("time: %f\n", time_spent);
	//======================================================================================================
	return 0;
}
