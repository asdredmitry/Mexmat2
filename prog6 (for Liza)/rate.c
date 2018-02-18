#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <sys/resource.h>
#include <sys/time.h>

#include "matrix.h"

int max(int q, int p);

long int get_full_time(void){
	struct timeval buf;
	gettimeofday(&buf, 0);
	return buf.tv_sec * 100 + buf.tv_usec/10000;
}

int max(int q, int p){
	if (q < p){
		return p;
	}
	return q;
}

double residualrate(matrix * M, double * X, double * Y){
	double ans = 0;
	double d;
	for (int i = 0; i < (M -> n); ++i){
		d = X[i];
		for (int j = 0; j < (M -> m); ++j){
			d -= (M -> A)[i * (M -> m) + j] * Y[j];
		}
		ans += d * d;
	}
	return sqrt(ans);
}

int funcmatrix(double * A, double * X, int N, int MyRank, int TotalNumberProcess){
	for (int i = MyRank; i < N; i += TotalNumberProcess){
		X[i / TotalNumberProcess] = 0;
		for (int j = 0; j < N; ++j){
			A[(i / TotalNumberProcess) * N + j] = 1. / (double) (max(i, j) + 2);
			X[(i / TotalNumberProcess)] += (1 - j % 2) * A[(i / TotalNumberProcess) * N + j];
		}
	}
	return 0;
}
