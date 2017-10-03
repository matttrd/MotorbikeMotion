//
//  main.c
//  EKF
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../EKF/EKFstruct.h"
#include "../../EKF/EKF.h"
#include "roll_angle_formula.h"



const char *getfield(char* line, int num)
{
    const char* tok;
    for (tok = strtok(line, ",");
         tok && *tok;
         tok = strtok(NULL, ",\r\n"))
    {
        if (!--num)
            return tok;
    }
    return NULL;
}
	


int main(int argc, const char * argv[]) {

	double x0[2] = {0, 0};
	double Q0[4] = {0.005, 0.0,
					0.0, 0.0001};
	double R0[1] = {0.2};
	
	double P0[4] = {1, 0,
					0, 1};
	
	EKF* roll_EKF = init_ekf(2, 1, 2, &x0[0], &Q0[0], &R0[0], &P0[0]);
	
	
	roll_EKF->compute_f_F = &create_f;
	roll_EKF->compute_h_H = &create_h;
	
	roll_EKF->parameters = init_moto_default();

    FILE* stream = fopen("/Users/matt/Dropbox/B/roll_angle.csv", "r");
    FILE* fp;
    fp = fopen("/Users/matt/Dropbox/B/roll_angle_test.txt", "w");
    double data[7];
    char line[1024];
    double dt, t_prev = 0.0;
    
    int i, row = 0;

    while (fgets(line, 1024, stream)){
        if (!row){ row = 1; continue; }
        for (i = 0; i < 7; i++){
            char* tmp = strdup(line);
            data[i] = atof(getfield(tmp, i + 1));
            free(tmp);
        }
        
        dt = data[3] - t_prev;
        t_prev = data[3];
        
        roll_EKF->u[0] = data[4];
        roll_EKF->u[1] = data[0];
        
        ekf_prediction(roll_EKF, dt);
        
        roll_EKF->z[0] = data[2];
        ekf_update(roll_EKF);
        
        fprintf(fp, "%f\n", roll_EKF->x[0]);
    }
    fclose(fp);
    fclose(stream);
	free_ekf(roll_EKF);
	
    return 0;
}
