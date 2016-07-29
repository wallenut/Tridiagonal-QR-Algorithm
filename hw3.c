/*
	CS440 Numerical Computation

	Homework 3: QR Tridiagonal Algorithm

	Allen Wang

	4/28/16
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

void printm(double *a, int m, int n)  //prints out a matrix to test at the end of main.
{ 
  int i, k;
  for (i = 0; i < m; i++) {
    for (k = 0; k < n; k++) {
      printf("\t%e\t", a[i * m + k]);
    }
    printf("\n");
  }
}

void mushift(int n, double mu, double shiftedb[n][n], int add){ //given a 2d matrix, adds or subtracts mu from diagonals
	int i, j;
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			if (i == j){
				if (add == 1){
					shiftedb[i][j] += mu;
				}
				else{
					shiftedb[i][j] -= mu;
				}
			}
		}
	}

}

void twodimension(int n, double* b, double result[n][n], int direction){ //converts an array into a 2d matrix
	int i, j, index;
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			if (direction == 1){
				result[i][j] = b[i*n + j];
			}
			else{
				b[i*n + j] = result[i][j];
			}
		}
	}
}

void qr_symmetric(double *a, int n, double *b)
{
    
    int rows, cols, k; // counters for loops of hte givens rotatiosn and checking for zeroed out terms and eigenvalues


    double b_rotation[n][n];//2d arrays that store the results of givens rotations, mu shifts, and transposes
    double thetaarray[n];
    double finalB[n][n];
    double bshifted[n][n];


    double b_val, d_val;
    double precision = 1e-16; //precision that Kirill says to go for (I think so)
    

    int i, j, b_d_colrow;

    twodimension(n, a, finalB, 1);

    rows = 0;

    double theta;
    
    while( rows < n)
    {
        b_d_colrow = rows + 1;
        while(fabs(finalB[rows][b_d_colrow]) > precision) //term is not yet zeroed out next to the diagonal
        {
            double mu;
            mu = finalB[rows][rows];

            for (i = 0; i < n; i++){
                memcpy(bshifted[i], finalB[i], sizeof(double) * n);
            }

            mushift(n, mu, bshifted, 0);

            for (i = 0; i < n; i++){
                memcpy(b_rotation[i], bshifted[i], sizeof(double) * n); //creating copies to re add the mu shift
            }

            //multiplying the Q matrices
            for (cols = n - 2; cols >= rows; cols -= 1)
            {
                b_d_colrow = cols + 1;
                b_val = b_rotation[cols][b_d_colrow]; // calculate theta based on B and D of an A B C D matrix being rotated
                d_val = b_rotation[b_d_colrow][b_d_colrow] + precision; 
                theta = atan(b_val / d_val);
                thetaarray[cols] = theta; // store into array
                for (k = rows; k < n; k+=1)
                {
                    double givens1; 
                    double givens2; 
                    givens1 = (cos(theta) * b_rotation[cols][k])  - (sin(theta) * b_rotation[b_d_colrow][k]); //the sum of each B, D pair across 2 whole row
                    givens2 = (sin(theta) * b_rotation[cols][k])  + (cos(theta)* b_rotation[b_d_colrow][k]);
                    b_rotation[cols][k] = givens1;
                    b_rotation[b_d_colrow][k] = givens2;
                    
                }
            }

            //multiplying Qtranspose from thetaarray

            for (cols = n - 2; cols >= rows; cols -= 1)
            { 
                b_d_colrow = cols + 1;
                theta = thetaarray[cols];
                for (k = rows; k < n; k += 1) //same as above but now vertically with swapped rows/columns and indices
                { 
                    double givens1; 
                    double givens2; 
                    givens1 = (cos(theta) * b_rotation[k][cols]) - (sin(theta) * b_rotation[k][b_d_colrow]);
                    givens2 = (sin(theta) * b_rotation[k][cols]) + (cos(theta) * b_rotation[k][b_d_colrow]);
                    b_rotation[k][cols] = givens1;
                    b_rotation[k][b_d_colrow] = givens2;
                    
                }
            }


            for (i = 0; i < n; i++){
                memcpy(finalB[i], b_rotation[i], sizeof(double) * n);

            }

            mushift(n, mu, finalB, 1);

        }
        rows += 1;
    }

    twodimension(n, b, finalB, 0);
    return;
}



int main (int argc, char* argv[]){ //multiple tests

	int n = 4;
	double * a = malloc(sizeof(double)* n * n);
	//double* result = malloc(sizeof(double) * n * n);
	a[0] = 1;
	a[1] = 4;
	a[2] = 0;
	a[3] = 0;
	a[4] = 3;
	a[5] = 4;
	a[6] = 1;
	a[7] = 0;
	a[8] = 0;
	a[9] = 2;
	a[10] = 3;
	a[11] = 4;
	a[12] = 0;
	a[13] = 0;
	a[14] = 1;
	a[15] = 3;

	double * result = malloc(sizeof(double) * n * n);
	printm(a, n, n);
	printf("\n");


	qr_symmetric(a, n, result);

	printm(a, n, n);
	printf("\n");
	printm(result, n, n);	
	printf("\n");

	double * test1 = malloc(sizeof(double) * 3 * 3);
	double * test2 = malloc(sizeof(double) * 3 * 3);

	test1[0] = 1;
	test1[1] = 4;
	test1[2] = 0;
	test1[3] = 3;
	test1[4] = 3;
	test1[5] = 4;
	test1[6] = 0;
	test1[7] = 2;
	test1[8] = 4;
	
	printm(test1, 3, 3);
	printf("\n");


	qr_symmetric(test1, 3, test2);

	printm(test1, 3, 3);
	printf("\n");
	printm(test2, 3, 3);
	printf("\n");

	



}
