//
//  main.c
//  AHRS
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../cblas/include/cblas.h"

#include "AHRS.h"




//for testing
#include "AHRSstruct.h"
#include "mechanization.h"
#include "linalg.h"

#define type	double

//int FP = 0;
//
//const char *getfield(char* line, int num)
//{
//	const char* tok;
//	for (tok = strtok(line, ",");
//		 tok && *tok;
//		 tok = strtok(NULL, ",\n"))
//	{
//		if (!--num)
//			return tok;
//	}
//	return NULL;
//}

//int main(){
//	float data[9];
//	
//	FILE* stream = fopen("/Users/michele/Desktop/c_code/AHRS/AHRS/install/test.txt", "r");
//	
//	char line[1024];
//	
//	while (fgets(line, 1024, stream)){
//
//		if (FP == 1){
//			int i = 0;
//			for (i = 0; i < 9; i++){
//				char* tmp = strdup(line);
//				data[i] = atof(getfield(tmp, i + 1));
//				free(tmp);
//			}
//			
//			EKF_AHRS_update(data[0], data[1], data[2], 9.81 * data[3], 9.81 * data[4], 9.81 * data[5], data[6], data[7], data[8]);
//			
//		} else {FP = 1;}
//
//		
//
//	}
//}


int main(int argc, const char * argv[]) {
	// insert code here...
	
	float *matrix1 = malloc(sizeof(matrix1) * 4);
	float *matrix2 = malloc(sizeof(matrix2) * 4);
	float *matrix3 = malloc(sizeof(matrix3) * 4);
	
	matrix1[0] = 1.0; matrix1[1] = 2.0;
	matrix1[2] = 3.0; matrix1[3] = 4.0;
	
	matrix2[0] = 5.0; matrix2[1] = 6.0;
	matrix2[2] = 7.0; matrix2[3] = 8.0;
	
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 2, 2, 1.0, matrix1, 2, matrix2, 2, 0.0, matrix3, 2);
	
	printf("%f %f\n%f %f\n", matrix3[0], matrix3[1], matrix3[2], matrix3[3]);

	free(matrix1);
	free(matrix2);
	free(matrix3);
	
    return 0;
}
