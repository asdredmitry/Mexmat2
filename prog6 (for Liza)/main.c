#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <string.h>

#include "matrix.h"

#define EPS 1e-20

int main(int argc, char ** argv){
	int N;
	double * A;
	double * X;
	double * Y;
	int * N_;
	double coeff;
	long int st, en;
	int TotalNumberProcess, MyRank;
	MPI_Status status;

	int z, main_i, x = 0, y = 0;
	int tmpn, NP_TO, NP_FROM;
	int * tmpN;
	double d, tmp, tmp1, tmp2;
	double * BUFER_TO;
	double * BUFER_FROM;
	double rate;
	int start_rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &TotalNumberProcess);
	MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
	if (MyRank == 0){
		printf("TotalNumberProcess: %d, MyRank: %d\n", TotalNumberProcess, MyRank);
	}
	if (argc > 2){
		printf("bad nomber of arguments\n");
	}
	N = atoi(argv[1]);
	if (MyRank == 0){
		printf("Martix size: %d\n", N);
	}
	BUFER_TO = (double *) malloc(N * sizeof(double));
	BUFER_FROM = (double *) malloc(N * sizeof(double));
	tmpN = (int *) malloc(N * sizeof(int));
	if (MyRank == 0){
		X = (double *) malloc(N * sizeof(double));
	}
	if (MyRank < N % TotalNumberProcess){
		tmpn = N / TotalNumberProcess + 1;
	}
	else{
		tmpn = N / TotalNumberProcess;
	}
	Y = (double *) malloc(tmpn * sizeof(double));
	A = (double *) malloc(N * tmpn * sizeof(double));
	funcmatrix(A, Y, N, MyRank, TotalNumberProcess);
//start print
	x = 0;
	z = 0;
	st = get_full_time();
	while (x < N){
		if (x % TotalNumberProcess == MyRank){
			for (y = 0; y < N; ++y){
				printf("%lf ", A[z * N + y]);
				if (y == 10){
					break;
				}
			}
			printf("| x = %d, R = %d | Y[%d] = %lf", x, MyRank, x, Y[z]);
			++z;
			printf("\n");
		}
		++x;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		if (x == 10){
			break;
		}
	}
	if (MyRank == 0){
		printf("\n\n\n");
	}
//end print

	x = 0;
	y = 0;

	while (x < N){
		main_i = x + ((TotalNumberProcess - ((x - MyRank) % TotalNumberProcess)) % TotalNumberProcess);
		for (int i = x + ((TotalNumberProcess - ((x - MyRank) % TotalNumberProcess)) % TotalNumberProcess); i < N; i += TotalNumberProcess){
			if (fabs(A[(i / TotalNumberProcess) * N + y]) > fabs(A[(main_i / TotalNumberProcess) * N + y])){
				main_i = i;
			}
		}
		if (main_i >= N){
			main_i = -1;
		}
		NP_TO = x % TotalNumberProcess;
		if (MyRank == 0){
			NP_FROM = 0;
			d = fabs(A[(main_i / TotalNumberProcess) * N + y]);
			for (int i = 1; i < TotalNumberProcess; ++i){
				MPI_Recv(&tmpn, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
				MPI_Recv(&tmp, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &status);
				if (tmpn >= 0){
					if (main_i == -1){
						main_i = tmpn;
						NP_FROM = i;
						d = tmp;
					}
					else{
						if (tmp > d){
							main_i = tmpn;
							NP_FROM = i;
							d = tmp;
						}
						else{
							if (fabs(tmp - d) < EPS){
								if (main_i > tmpn){
									main_i = tmpn;
									NP_FROM = i;
								}
							}
						}
					}
				}
			}
		}
		else{
			if (main_i >= 0){
				tmpn = main_i;
				tmp = fabs(A[(main_i / TotalNumberProcess) * N + y]);
				MPI_Send(&tmpn, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
				MPI_Send(&tmp, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
			}
			else{
				tmpn = -1;
				tmp = -1.;
				MPI_Send(&tmpn, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
				MPI_Send(&tmp, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
			}
		}
		if (MyRank == 0){
			for (int i = 1; i < TotalNumberProcess; ++i){
				MPI_Send(&NP_FROM, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
			}
		}
		else{
			MPI_Recv(&NP_FROM, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);
		}
		if (MyRank == 0){
			for (int i = 1; i < TotalNumberProcess; ++i){
				MPI_Send(&main_i, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
			}
		}
		else{
			MPI_Recv(&main_i, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		}
		if (MyRank == NP_FROM){
			d = A[(main_i / TotalNumberProcess) * N + y];
			for (int i = y; i < N; ++i){
				A[(main_i / TotalNumberProcess) * N + i] /= d;
			}
			Y[x / TotalNumberProcess] /= d;
		}
		if (MyRank == NP_FROM){
			if (MyRank == NP_TO){
				for (int j = y; j < N; ++j){
					BUFER_FROM[j] = A[(main_i / TotalNumberProcess) * N + j];
				}
				for (int i = y; i < N; ++i){
					tmp = A[(x / TotalNumberProcess) * N + i];
					A[(x / TotalNumberProcess) * N + i] = A[(main_i / TotalNumberProcess) * N + i];
					A[(main_i / TotalNumberProcess) * N + i] = tmp;
				}
				tmp = Y[x / TotalNumberProcess];
				Y[x / TotalNumberProcess] = Y[main_i / TotalNumberProcess];
				Y[main_i / TotalNumberProcess] = tmp;
			}
			else{
				for (int j = y; j < N; ++j){
					BUFER_FROM[j] = A[(main_i / TotalNumberProcess) * N + j];
				}
				tmp1 = Y[main_i / TotalNumberProcess];
				MPI_Send(BUFER_FROM + y, N - y, MPI_DOUBLE, NP_TO, x * TotalNumberProcess * 2 + 1, MPI_COMM_WORLD);
				MPI_Send(&tmp1, 1, MPI_DOUBLE, NP_TO, 0, MPI_COMM_WORLD);
				//get "buffer to"
				MPI_Recv(BUFER_TO + y, N - y, MPI_DOUBLE, NP_TO, x * TotalNumberProcess * 2 + 2, MPI_COMM_WORLD, &status);
				MPI_Recv(&tmp2, 1, MPI_DOUBLE, NP_TO, 1, MPI_COMM_WORLD, &status);
				for (int j = y; j < N; ++j){
					A[(main_i / TotalNumberProcess) * N + j] = BUFER_TO[j];
				}
				Y[main_i / TotalNumberProcess] = tmp2;
			}
		}
		if (MyRank == NP_TO){
			if (MyRank != NP_FROM){
				for (int j = y; j < N; ++j){
					BUFER_TO[j] = A[(x / TotalNumberProcess) * N + j];
				}
				tmp2 = Y[x / TotalNumberProcess];
				//get "buffer from"
				MPI_Recv(BUFER_FROM + y, N - y, MPI_DOUBLE, NP_FROM, x * TotalNumberProcess * 2 + 1, MPI_COMM_WORLD, &status);
				MPI_Recv(&tmp1, 1, MPI_DOUBLE, NP_FROM, 0, MPI_COMM_WORLD, &status);
				//send "buffer to"
				MPI_Send(BUFER_TO + y, N - y, MPI_DOUBLE, NP_FROM, x * TotalNumberProcess * 2 + 2, MPI_COMM_WORLD);
				MPI_Send(&tmp2, 1, MPI_DOUBLE, NP_FROM, 1, MPI_COMM_WORLD);
				for (int j = y; j < N; ++j){
					A[(x / TotalNumberProcess) * N + j] = BUFER_FROM[j];
				}
				Y[x / TotalNumberProcess] = tmp1;	
			}
		}
		if (MyRank == NP_TO){
			d = Y[x / TotalNumberProcess];
		}
		MPI_Bcast(BUFER_FROM + y, N - y, MPI_DOUBLE, NP_TO, MPI_COMM_WORLD);
		MPI_Bcast(&d, 1, MPI_DOUBLE, NP_TO, MPI_COMM_WORLD);
		for (int i = MyRank; i < N; i += TotalNumberProcess){
			if (i != y){
				coeff = A[(i / TotalNumberProcess) * N + y] / BUFER_FROM[y];
				for (int j = y; j < N; ++j){
					A[(i / TotalNumberProcess) * N + j] -= coeff * BUFER_FROM[j];
				}
				Y[i / TotalNumberProcess] -= coeff * d;
			}
		}
		++x;
		++y;
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);

//start print
	x = 0;
	z = 0;
	while (x < N){
		if (x % TotalNumberProcess == MyRank){
			for (y = 0; y < N; ++y){
				printf("%lf ", A[z * N + y]);
				if (y == 10){
					break;
				}
			}
			printf("| x = %d, R = %d | Y[%d] = %lf", x, MyRank, x, Y[z]);
			++z;
			printf("\n");
		}
		++x;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		if (x == 10){
			break;
		}
	}
	if (MyRank == 0){
		printf("\n\n\n");
	}
//end print

	x = 0;
	z = 0;
	while (x < N){
		if (x % TotalNumberProcess == MyRank){
			if (MyRank != 0){
				MPI_Send(Y + z, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			}
			else{
				X[x] = Y[z];
			}
			++z;
		}
		else{
			if (MyRank == 0){
				MPI_Recv(X + x, 1, MPI_DOUBLE, x % TotalNumberProcess, 0, MPI_COMM_WORLD, &status);
			}
		}
		++x;
		MPI_Barrier(MPI_COMM_WORLD);
	}

	en = get_full_time();

	for (int i = 0; i < TotalNumberProcess; ++i){
		if (i == MyRank){
			printf("time of %d: %ld\n", MyRank, en - st);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (MyRank == 0){
		rate = 0;
		printf("\n\nanswer: ");
		for (int i = 0; i < N; ++i){
			printf("%lf | ", X[i]);
			rate += (((double) (1 - i % 2)) - X[i]) * (((double) (1 - i % 2)) - X[i]);
			if (i == 10){
				break;
			}
		}
		printf("\n\n");
		rate = sqrt(rate);
		printf("residualrate: %e\n", rate);
	}

	if (MyRank == 0){
		//free(N_);
		free(X);
	}
	free(tmpN);
	free(BUFER_FROM);
	free(BUFER_TO);
	free(Y);
	//free(TotalY);
	//delete(&M);
	MPI_Finalize();
	return 0;
}
