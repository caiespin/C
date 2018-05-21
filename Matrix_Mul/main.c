// CMPE110 HA5: Locality Detective
// Author:  CMPE110 Spring18 staff
// 5/20/2018
// skeleton code of the matrix multiplication
//
//
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>

#define SIZE 1024

__uint64_t A[SIZE][SIZE];
__uint64_t B[SIZE][SIZE];
__uint64_t C[SIZE][SIZE];

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

int verify(__uint64_t C[][SIZE], __uint64_t D[][SIZE])
{
	int r, c;

	for (c = 0; c < SIZE; c++) {
		for (r = 0; r < SIZE; r++) {
			if (C[r][c] != D [r][c]) {
				printf("error!\n");
				goto out;
			}

		}
	}
	return 0;

out:
	return -1;
}

void matmul(__uint64_t A[][SIZE], __uint64_t B[][SIZE])
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

int main(int argc, char *argv)
{
	clock_t t;
	double time_taken;

	init(A, B);
	memset(C, 0, sizeof(__uint64_t) * SIZE * SIZE);

	t = clock();
	matmul(A, B);
	t = clock() - t;
	time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds

	printf("Matmul took %f seconds to execute \n", time_taken);
}
