#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "help.h"

#define MAX_OUTPUT_SIZE 5

double f(int i, int j)
{
	return i - j;
}

double *AllocVector(int size)
{
	return (double*)malloc(size * sizeof(double));
}

void InputMatrix(int n, double *a, double *b, int my_rank, int p)
{
	int i, j;
	int rows;
	double x;

	if (my_rank + 1 > n%p) rows = n/p;
	else rows = n/p + 1;

	for (i = 0; i < rows; i++)
	{
		for (x = 0.0, j = 0; j < n; j++)
		{
			a[i * n + j] = f(my_rank + p * i, j);
			if (!(j % 2)) x += a[i * n + j];
		}
		b[i] = x;
	}
}

void OutputMatrix(int n, double *a, double *b, double *x, int my_rank, int p)
{
	int i, j, m;
	MPI_Status status;

	m = (n < MAX_OUTPUT_SIZE) ? n : MAX_OUTPUT_SIZE;

	for (i = 0; i < m; i++)
	{
		if (my_rank == 0)
		{
			if (my_rank == i%p)
			{
				printf("| ");
				for (j = 0; j < m; j++)
					printf("%10.3g ", a[i/p * n + j]);
				printf("|   %10.3g\n", b[i/p]);
			}
			else
			{
				MPI_Recv(x, m + 1, MPI_DOUBLE, i%p, 0, MPI_COMM_WORLD, &status);
				printf("| ");
				for (j = 0; j < m; j++)
					printf("%10.3g ", x[j]);
				printf("|   %10.3g\n", x[m]);
			}
		}
		else if (my_rank == i%p)
		{
			for (j = 0; j < m; j++)
				x[j] = a[i/p * n + j];
			x[m] = b[i/p];
			MPI_Send(x, m + 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	}
}

void OutputVector(int n, double *b, double *x, double * y, int my_rank, int p)
{
	int i, m;
	MPI_Status status;
        m = n;
        //m = (n < MAX_OUTPUT_SIZE) ? n : MAX_OUTPUT_SIZE;
	for (i = 0; i < m; i++)
	{
		if (my_rank == 0)
		{
			if (my_rank == i%p)
                        {
				printf("%10.3g ", b[i/p]);
                                y[i] = b[i/p];
                        }
			else
			{
				MPI_Recv(x, 1, MPI_DOUBLE, i%p, 0, MPI_COMM_WORLD, &status);
				printf("%10.3g ", x[0]);
                                y[i] = x[0];
			}
		}
		else if (my_rank == i%p)
		{
			x[0] = b[i/p];
			MPI_Send(x, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	}
}
double SolutionAccuracy(int n, double* x)
{
	int i;
	double tmp;
	double rezult;

	rezult = 0.0;
	for (i = 0; i < n; ++i)
	{
		tmp = x[i];

		if (i % 2 == 0)
			tmp -= 1.0;

		rezult += tmp * tmp;
	}

	return sqrt(rezult);
}
