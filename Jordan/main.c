#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "help.h"
#include "sol.h"

int main(int argc, char **argv)
{
	int n;
	double *a;
	double *b;
	double *x;
	double t;
	int rows;
	int my_rank, p;
	int err1, err2;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if (argc == 2) n = atoi(argv[1]);
	else
	{
		if (my_rank == 0) printf("Not correct usage.\n");

		MPI_Finalize();

		return 1;
	}

	if (my_rank + 1 > n%p) rows = n/p;
	else rows = n/p + 1;

	a = AllocVector(rows * n);
	b = AllocVector(rows);
	x = AllocVector(n + 1);

	err1 = 0;
	err2 = 0;
	if (!(a && b && x)) err1 = 1;

	MPI_Allreduce(&err1, &err2, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	if (err2)
	{
		if (my_rank == 0) printf("Not enough memory.\n");

		if (a) free(a);
		if (b) free(b);
		if (x) free(x);

		MPI_Finalize();

		return 1;
	}

	InputMatrix(n, a, b, my_rank, p);

	if (my_rank == 0) printf("Matrix A:\n\n");
	OutputMatrix(n, a, b, x, my_rank, p);

	MPI_Barrier(MPI_COMM_WORLD);
	t = MPI_Wtime();

	SolveSystem(n, a, b, x, my_rank, p);

	MPI_Barrier(MPI_COMM_WORLD);
	t = MPI_Wtime() - t;

	if (my_rank == 0) printf("\nSolution:\n");
	OutputVector(n, b, x, my_rank, p);

	if (my_rank == 0) printf("\n\nSolution time = %e\n", t);

	free(a);
	free(b);
	free(x);

	MPI_Finalize();

	return 0;
}
