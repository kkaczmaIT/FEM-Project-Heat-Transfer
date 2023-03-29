#include"Gauss.h"
/*
	Functions to execute Gauss-Crout Elimination Method Start
*/

bool firstLevel(double** arr, const int n)
{
	double counter;
	for (int z = 0; z < n; z++)
	{
		if (arr[z][z] == 0)
		{
			std::cout << "Wystapilo zero na przekatnej\n";
			return 0;
		}
		for (int i = z + 1; i < n; i++)
		{
			counter = arr[i][z] / arr[z][z];
			for (int j = 0; j < n + 1; j++)
			{
				arr[i][j] = arr[i][j] - (counter * arr[z][j]);
			}
		}
	}
	return 1;
}

double* secondLevel(double** arr, const int n)
{
	double* xArr = new double[n];
	double* solArr = new double[n];
	for (int z = 0; z < n; z++)
	{
		xArr[z] = arr[z][n] / arr[z][z];
	}
	for (int i = n - 1; i >= 0; i--)
	{
		solArr[i] = arr[i][n];
		for (int k = i + 1; k < n; k++)
		{
			solArr[i] -= (arr[i][k] * xArr[k]);
		}
		solArr[i] /= arr[i][i];
		xArr[i] = solArr[i];
	}
	delete[] xArr;
	return solArr;
	delete[] solArr;
}

/*
	Functions to execute Gauss-Crout Elimination Method End
*/