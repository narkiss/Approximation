#include<math.h>

double f(double x)
{
	return exp(-x) * cos(x);
}

double P(double x, double* C, int n)
{
	double res = 0;

	for (int i = 0; i <= n; i++)
	{
		res += C[i] * pow(x, i);
	}

	return res;
}

void sort(double** A, double* B, int col, int n)
{
	int swap = col;
	double A1;
	
	for (int i = col; i < n; i++)
	{
		if (abs(A[swap][col]) < abs(A[i + 1][col])) swap = i + 1;
	}

	for (int i = 0; i <= n; i++)
	{
		A1 = A[col][i];
		A[col][i] = A[swap][i];
		A[swap][i] = A1;
	}

	A1 = B[col];
	B[col] = B[swap];
	B[swap] = A1;
}

void Gauss(double** A, double* C, double* B, int n, int m)
{
	int n1 = n;
	double A1;

	for (int k = 0; k < n1; k++)
	{
		sort(A, B, k, n1);

		for (int i = n1 - n; i < n1; i++)
		{
			if (A[i + 1][k] != 0)
			{
				A1 = A[i + 1][k];

				for (int j = n1 - n; j <= n1; j++)
				{
					A[i + 1][j] -= A[k][j] * A1 / A[k][k];
				}
				B[i + 1] -= B[k] * A1 / A[k][k];
			}
		}
		n--;
	}

	for (int i = n1; i >= 0; i--)
	{
		C[i] = B[i];

		for (int j = i + 1; j <= n1; j++)
		{
			C[i] -= A[i][j] * C[j];
		}

		C[i] /= A[i][i];
	}
}

double find_s(int n, int m, double* x, double* y, double* C)
{
	double S = 0;

	for (int i = 0; i <= m; i++)
	{
		S += (P(x[i], C, n) - y[i]) * (P(x[i], C, n) - y[i]);
	}

	return S;
}