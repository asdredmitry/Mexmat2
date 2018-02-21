#include <math.h>
#include <mpi.h>

#include "sol.h"

void SolveSystem(int n, double *a, double *b, double *x, int my_rank, int p)
{
	int i, j, k;
	int rows;
	double tmp;

	if (my_rank + 1 > n%p) rows = n/p;
	else rows = n/p + 1;

	for (i = 0; i < n; i++)
	{
		if (my_rank == i%p)
		{
			tmp = 1.0/a[i/p * n + i];
			for (j = i; j < n; j++)
				a[i/p * n + j] *= tmp;
			b[i/p] *= tmp;

			for (j = i; j < n; j++)
				x[j] = a[i/p * n + j];
			x[n] = b[i/p];

			MPI_Bcast(x, n + 1, MPI_DOUBLE, i%p, MPI_COMM_WORLD);
			for (j = i/p + 1; j < rows; j++)
			{
				tmp = a[j * n + i];
				for (k = i; k < n; k++)
					a[j * n + k] -= tmp * a[i/p * n + k];
				b[j] -= tmp * b[i/p];
			}
			for (j = 0; j < i/p; j++)
			{
				tmp = a[j * n + i];
				for (k = i; k < n; k++)
					a[j * n + k] -= tmp * a[i/p * n + k];
				b[j] -= tmp * b[i/p];
			}
		}
		else
		{
			MPI_Bcast(x, n + 1, MPI_DOUBLE, i%p, MPI_COMM_WORLD);
			for (j = 0; j < rows; j++)
			{
				tmp = a[j * n + i];
				for (k = i; k < n; k++)
					a[j * n + k] -= tmp * x[k];
				b[j] -= tmp * x[n];
			}
		}
	}
}
