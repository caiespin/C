// CMPE110 HA5: Locality Detective
// Author:  CMPE110 Spring18 staff
// Modified by Carlos Espinosa
// 5/20/2018
// The purpose of this code is to analyze different methods of optimization for a matrix multiplication kernel
// by means of leveraging spatial and temporal memory locality. The analyzed methods are:
// - Matrix Transpose
// - Matrix Block Tiling

#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#define SIZE 1024

typedef int bool;
#define TRUE  1
#define FALSE 0

__uint64_t A[SIZE][SIZE]; //Random filled Matrix A
__uint64_t B[SIZE][SIZE]; //Random filled Matrix B
__uint64_t T[SIZE][SIZE]; //Result of Transpose Calculation
__uint64_t C[SIZE][SIZE]; //Result of Base Naive Multiplication approach
__uint64_t D[SIZE][SIZE]; //Result of Transpose Multiplication approach
__uint64_t E[SIZE][SIZE]; //Result of Tiling Multiplication approach
__uint64_t F[SIZE][SIZE]; //Result of Combining Tiling and Transpose Multiplication approach

void init(__uint64_t A[][SIZE], __uint64_t B[][SIZE])
{
	int r, c;

	for (c = 0; c < SIZE; c++) {
		for (r = 0; r < SIZE; r++) {
			A[r][c] = rand();
			B[r][c] = rand();
		}
	}
}

void transpose(__uint64_t O[][SIZE], __uint64_t T[][SIZE])
{
    int rowO, colT;
    for (rowO = 0; rowO < SIZE; rowO++){
        for(colT = 0; colT < SIZE; colT++){
            T[colT][rowO]=O[rowO][colT];
        }
    }
}

void verify(__uint64_t C[][SIZE], __uint64_t D[][SIZE])
{
	int r, c, ok=1;
    bool stop = FALSE;

	for (c = 0; (c < SIZE) && !stop; c++) {
		for (r = 0; (r < SIZE) && !stop; r++) {
			if (C[r][c] != D [r][c]) {
				ok = 0;
				stop = TRUE;
				printf("Matrices do not match.\n");
			}

		}
	}
	if (ok==1)
        printf("Matrices are equal.\n");
}

void matmulBase(__uint64_t A[][SIZE], __uint64_t B[][SIZE])
{
	int rowA, colB, idx;

	for (rowA = 0; rowA < SIZE; rowA++) {
		for (colB = 0; colB < SIZE; colB++) {
			for (idx = 0; idx < SIZE; idx++) {
				C[rowA][colB] += A[rowA][idx] * B[idx][colB];
			}
		}
	}
}

void matmulTrans(__uint64_t A[][SIZE], __uint64_t B[][SIZE])
{
	int rowA, colB, idx;

	for (rowA = 0; rowA < SIZE; rowA++) {
		for (colB = 0; colB < SIZE; colB++) {
			for (idx = 0; idx < SIZE; idx++) {
				D[rowA][colB] += A[rowA][idx] * B[colB][idx];
			}
		}
	}
}

void matmulTiling(__uint64_t A[][SIZE], __uint64_t B[][SIZE], int TILE)
{
    int i, j, k;
	int rowA, colB, idx;
    //This initial loops will move the multiplication over each Tile
	for (i = 0; i < SIZE; i+=TILE) {
		for (j = 0; j < SIZE; j+=TILE) {
			for (k = 0; k < SIZE; k+=TILE) {
			    //Regular matrix multiplication in each Tile
                for (rowA = i; rowA < i+TILE; rowA++) {
                    for (colB = j; colB < j+TILE; colB++) {
                        for (idx = k; idx < k+TILE; idx++) {
                            E[rowA][colB] += A[rowA][idx] * B[idx][colB];
                        }
                    }
                }
			}
		}
	}
}

void matmulTileTrans(__uint64_t A[][SIZE], __uint64_t B[][SIZE], int TILE)
{
    int i, j, k;
	int rowA, colB, idx;
    //This initial loops will move the multiplication over each Tile
	for (i = 0; i < SIZE; i+=TILE) {
		for (j = 0; j < SIZE; j+=TILE) {
			for (k = 0; k < SIZE; k+=TILE) {
			    //Matrix multiplication by transpose in each Tile
                for (rowA = i; rowA < i+TILE; rowA++) {
                    for (colB = j; colB < j+TILE; colB++) {
                        for (idx = k; idx < k+TILE; idx++) {
                            F[rowA][colB] += A[rowA][idx] * B[colB][idx];
                        }
                    }
                }
			}
		}
	}
}

void clear(__uint64_t C[][SIZE])
{
	int r, c;

	for (c = 0; c < SIZE; c++) {
		for (r = 0; r < SIZE; r++) {
			C[r][c] = 0;
		}
	}
}

void printM(__uint64_t M[][SIZE])
{
    int r, c;

	for (c = 0; c < SIZE; c++) {
		for (r = 0; r < SIZE; r++) {
            printf("%llu  ", M[r][c]);
		}
		printf("\n");
	}
}

int main(int argc, char *argv)
{
	clock_t t;
	double time_takenBase, time_takenTrans, time_takenTiling, time_takenTileTrans;
	double x;
	int Tile;

	init(A, B);
	memset(C, 0, sizeof(__uint64_t) * SIZE * SIZE);
    transpose(B, T);
	memset(D, 0, sizeof(__uint64_t) * SIZE * SIZE);
	memset(E, 0, sizeof(__uint64_t) * SIZE * SIZE);
	memset(F, 0, sizeof(__uint64_t) * SIZE * SIZE);

	t = clock();
	matmulBase(A, B);
	t = clock() - t;
	time_takenBase = ((double)t)/CLOCKS_PER_SEC; // in seconds
	printf("--------------------------------------------------------------------------------------------------\n");
    printf("The Execution time for the Matrix Multiplication Base Approach is: %f\n", time_takenBase);

	t = clock();
	matmulTrans(A, T);
	t = clock() - t;
	time_takenTrans = ((double)t)/CLOCKS_PER_SEC; // in seconds
	printf("\n--------------------------------------------------------------------------------------------------\n");
	printf("The Execution time for the Matrix Multiplication Transpose Approach is: %f\n", time_takenTrans);
    printf("Results verified, ");
    verify(C, D);
    printf("\n--------------------------------------------------------------------------------------------------\n");
    for (x = 0; x <= 10; x++)
    {
        Tile=(int) pow(2, x);
        t = clock();
        matmulTiling(A, B, Tile);
        t = clock() - t;
        time_takenTiling = ((double)t)/CLOCKS_PER_SEC;
        printf("The Execution time for the Matrix Multiplication Tiling Approach, with a %dx%d Tile is: %f\n", Tile, Tile, time_takenTiling);
        printf("Results verified, ");
        verify(C, E);
        clear(E);
    }
    printf("\n--------------------------------------------------------------------------------------------------\n");
    for (x = 0; x <= 10; x++)
    {
        Tile=(int) pow(2, x);
        t = clock();
        matmulTileTrans(A, T, Tile);
        t = clock() - t;
        time_takenTileTrans = ((double)t)/CLOCKS_PER_SEC;
        printf("The Execution time for the Matrix Multiplication Tile-Trans Approach, with a %dx%d Tile is: %f\n", Tile, Tile, time_takenTileTrans);
        printf("Results verified, ");
        verify(C, F);
        clear(F);
    }
}
