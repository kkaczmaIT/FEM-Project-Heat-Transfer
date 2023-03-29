#include<iostream>
#include<cmath>

using namespace std;

double fn1(double x) 
{
	return -5 * pow(x, 2) + 2 * x - 8;
}

double fn2(double x, double y)
{
	return 3 * pow(x, 2) * y + 2 * x * pow(y, 2) + 2;
}

double fn3(double x) 
{
	return pow(x, 2) - 2 * x + 2;
}

double gaussLegendre1D(double a, double b, int n, double (fn)(double))
{
	double result = 0.0;
	if (n == 2)
	{
		double x2[2][2];
		x2[0][0] = -0.57735;
		x2[0][1] = 1;
		x2[1][0] = 0.57735;
		x2[1][1] = 1;
		double t[2];
		for (int i = 0; i < n; i++)
		{
			t[i] = (a + b) / 2 + (b - a) / 2 * x2[i][0];
			result += x2[i][1] * fn(t[i]);
		}
		result *= ((b - a) / 2);

	}
	else if (n == 3)
	{
		double x2[3][2];
		x2[0][0] = -pow(3.0 / 5.0, 0.5); 
		x2[0][1] = 5.0 / 9.0;
		x2[1][0] = 0;
		x2[1][1] = 8.0 / 9.0;
		x2[2][0] = pow(3.0 / 5.0, 0.5);
		x2[2][1] = 5.0 / 9.0;
		double t[3];
		for (int i = 0; i < n; i++)
		{
			t[i] = (a + b) / 2 + (b - a) / 2 * x2[i][0];
			result += x2[i][1] * fn(t[i]);
		}
		result *= ((b - a) / 2);
	}
	else if (n == 4)
	{
		double x2[4][2];
		x2[0][0] = -0.861136;
		x2[0][1] = 0.347855;
		x2[1][0] = -0.339981;
		x2[1][1] = 0.652145;
		x2[2][0] = 0.339981;
		x2[2][1] = 0.652145;
		x2[3][0] = 0.861136;
		x2[3][1] = 0.347855;
		double t[4];
		for (int i = 0; i < n; i++)
		{
			t[i] = (a + b) / 2 + (b - a) / 2 * x2[i][0];
			result += x2[i][1] * fn(t[i]);
		}
		result *= ((b - a) / 2);

	}
	return result;
}

double gaussLegendre2D(double a, double b, int n, double (fn)(double, double))
{
	double result = 0.0;
	if (n == 2)
	{
		double x2[2][2];
		x2[0][0] = -0.57735;
		x2[0][1] = 1;
		x2[1][0] = 0.57735;
		x2[1][1] = 1;
		double t[2][2]; 
		for (int i = 0; i < n; i++)
		{
			t[i][0] = (a + b) / 2 + (b - a) / 2 * x2[i][0];
			for (int j = 0; j < n; j++)
			{
				t[j][1] = (a + b) / 2 + (b - a) / 2 * x2[j][0];
				result += x2[i][1] * x2[j][1] * fn(t[i][0], t[j][1]);
			}
		}
		result *= ((b - a) / 2);

	}
	else if (n == 3)
	{
		double x2[3][2];
		x2[0][0] = -pow(3.0 / 5.0, 0.5);
		x2[0][1] = 5.0 / 9.0;
		x2[1][0] = 0;
		x2[1][1] = 8.0 / 9.0;
		x2[2][0] = pow(3.0 / 5.0, 0.5);
		x2[2][1] = 5.0 / 9.0;
		double t[3][2];
		for (int i = 0; i < n; i++)
		{
			t[i][0] = (a + b) / 2 + (b - a) / 2 * x2[i][0];
			for (int j = 0; j < n; j++)
			{
				t[j][1] = (a + b) / 2 + (b - a) / 2 * x2[j][0];
				result += x2[i][1] * x2[j][1] * fn(t[i][0], t[j][1]);
			}
		}
		result *= ((b - a) / 2);
	}
	/*else if (n == 4)
	{
		double x2[4][2];
		x2[0][0] = -0.861136;
		x2[0][1] = 0.347855;
		x2[1][0] = -0.339981;
		x2[1][1] = 0.652145;
		x2[2][0] = 0.339981;
		x2[2][1] = 0.652145;
		x2[3][0] = 0.861136;
		x2[3][1] = 0.347855;
		double t[4];
		for (int i = 0; i < n; i++)
		{
			t[i] = (a + b) / 2 + (b - a) / 2 * x2[i][0];
			result += x2[i][1] * fn(t[i]);
		}
		result *= ((b - a) / 2);

	}*/
	return result;
}

int main()
{
	cout << "Funkcja jednej zmiennej dla 2 punktow: " << gaussLegendre1D(-1, 1, 2, fn1) << endl;
	cout << "Funkcja jednej zmiennej dla 3 punktow: " << gaussLegendre1D(-1, 1, 3, fn1) << endl;
	cout << "Funkcja dwoch zmiennych dla 2 punktow: " << gaussLegendre2D(-1, 1, 2, fn2) << endl;
	cout << "Funkcja dwoch zmiennych dla 3 punktow: " << gaussLegendre2D(-1, 1, 3, fn2) << endl;

	cout << "21 pazdziernika lekcja: \n";
	cout << "Funkcja jednej zmiennej dla 2 punktow: " << gaussLegendre1D(-2, 6, 2, fn3) << endl;
	cout << "Funkcja jednej zmiennej dla 3 punktow: " << gaussLegendre1D(-2, 6, 3, fn3) << endl;
	return 0;
}