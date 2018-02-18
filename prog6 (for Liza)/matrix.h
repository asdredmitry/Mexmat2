typedef struct{
	int n, m;
	double * A;
} matrix;

int funcmatrix(double * A, double * X, int N, int MyRank, int TotalNumberProcess);
long int get_full_time(void);
