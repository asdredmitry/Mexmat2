#ifndef __HELP_H__
#define __HELP_H__

double f(int i, int j);

double *AllocVector(int size);

void InputMatrix(int n, double *a, double *b, int my_rank, int p);

void OutputMatrix(int n, double *a, double *b, double *x, int my_rank, int p);

void OutputVector(int n, double *b, double *x, int my_rank, int p);

#endif
