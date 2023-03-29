#include<iostream>
#include<cmath>
using namespace std;

//Prototype
void displayMatrix(double** matrix, int row, int col);
double** transpositionMatrix(double** matrix);
double** multplicationMatrixByScalar(double** matrix, double scalar, int row, int col);
double** sumMatrix(double** matrix1, double** matrix2, int row, int col);
double** transpositionMatrix(double** matrix);
double** multiplicationMatrix(double* verMatrix, double* horMatrix);

double N1(double ksi, double eta)
{
	return 0.25 * (1 - ksi) * (1 - eta);
}

double N2(double ksi, double eta)
{
	return 0.25 * (1 + ksi) * (1 - eta);
}

double N3(double ksi, double eta)
{
	return 0.25 * (1 + ksi) * (1 + eta);
}

double N4(double ksi, double eta)
{
	return 0.25 * (1 - ksi) * (1 + eta);
}

double N1ksi(double eta)
{
	return -0.25 * (1 - eta);
}

double N2ksi(double eta)
{
	return 0.25 * (1 - eta);
}

double N3ksi(double eta)
{
	return 0.25 * (1 + eta);
}

double N4ksi(double eta)
{
	return -0.25 * (1 + eta);
}

double N1eta(double ksi)
{
	return -0.25 * (1 - ksi);
}

double N2eta(double ksi)
{
	return -0.25 * (1 + ksi);
}

double N3eta(double ksi)
{
	return 0.25 * (1 + ksi);
}

double N4eta(double ksi)
{
	return 0.25 * (1 - ksi);
}




class Element4
{
private:
	int numberNodes;
	int numberPc;
	double* weightsEta;
	double* weightsKsi;
	double* detJ;
	double* eta;
	double* ksi;
	double** dNdksi;
	double** dNdeta;
	double** dNdx;
	double** dNdy;
	double*** matrixHPc;
	double** matrixH;
public:
	friend void displayMatrix(double** matrix, int row, int col);
	friend double** transpositionMatrix(double** matrix);
	friend double** multplicationMatrixByScalar(double** matrix, double scalar, int row, int col);
	friend double** sumMatrix(double** matrix1, double** matrix2, int row, int col);
	friend double** transpositionMatrix(double** matrix);
	friend double** multiplicationMatrix(double* verMatrix, double* horMatrix);
	Element4(int numberNodes)
	{
		this->numberNodes = numberNodes;
		this->numberPc = this->numberNodes * this->numberNodes;
		detJ = new double[this->numberPc];
		weightsEta = new double[this->numberPc];
		weightsKsi = new double[this->numberPc];
		eta = new double[this->numberPc];
		ksi = new double[this->numberPc];
		dNdksi = new double* [this->numberPc];
		dNdeta = new double* [this->numberPc];
		dNdx = new double* [this->numberPc];
		dNdy = new double* [this->numberPc];
		matrixHPc = new double** [this->numberPc];
		matrixH = new double* [this->numberPc];
		for (int i = 0; i < this->numberPc; i++)
		{
			dNdksi[i] = new double[4];
			dNdeta[i] = new double[4];
			dNdy[i] = new double[4];
			dNdx[i] = new double[4];
			matrixHPc[i] = new double* [4];
			matrixH[i] = new double[4];
			for (int j = 0; j < this->numberPc; j++)
			{
				matrixHPc[i][j] = new double[4];
			}
			
		}
		for (int row = 0; row < 4; row++)
			for (int col = 0; col < 4; col++)
				matrixH[row][col] = 0;
		if (this->numberNodes == 2)
		{
			weightsEta[0] = weightsEta[1] = weightsEta[2] = weightsEta[3] = weightsKsi[0] = weightsKsi[1] = weightsKsi[2] = weightsKsi[3] = 1;
			this->eta[0] = -1.0 / sqrt(3.0);
			this->eta[1] = -1.0 / sqrt(3.0);
			this->eta[2] = 1.0 / sqrt(3.0);
			this->eta[3] = 1.0 / sqrt(3.0);

			this->ksi[0] = -1.0 / sqrt(3.0);
			this->ksi[1] = 1.0 / sqrt(3.0);
			this->ksi[2] = -1.0 / sqrt(3.0);
			this->ksi[3] = 1.0 / sqrt(3.0);
		}
		else if (this->numberNodes == 3)
		{
			this->eta[0] = -sqrt(3.0 / 5.0);
			this->eta[1] = -sqrt(3.0 / 5.0);
			this->eta[2] = -sqrt(3.0 / 5.0);
			this->eta[3] = 0;
			this->eta[4] = 0;
			this->eta[5] = 0;
			this->eta[6] = sqrt(3.0 / 5.0);
			this->eta[7] = sqrt(3.0 / 5.0);
			this->eta[8] = sqrt(3.0 / 5.0);



			this->ksi[0] = -sqrt(3.0 / 5.0);
			this->ksi[1] = 0;
			this->ksi[2] = sqrt(3.0 / 5.0);
			this->ksi[3] = -sqrt(3.0 / 5.0);
			this->ksi[4] = 0;
			this->ksi[5] = sqrt(3.0 / 5.0);
			this->ksi[6] = -sqrt(3.0 / 5.0);
			this->ksi[7] = 0;
			this->ksi[8] = sqrt(3.0 / 5.0);

			weightsEta[0] = weightsEta[1] = weightsEta[2] = weightsEta[6] = weightsEta[7] = weightsEta[8] = 5.0 / 9.0;

			weightsEta[3] = weightsEta[4] = weightsEta[5] = 8.0 / 9.0;

			weightsKsi[0] = weightsKsi[3] = weightsKsi[6] = weightsKsi[2] = weightsKsi[5] = weightsKsi[8] = 5.0 / 9.0;

			weightsKsi[1] = weightsKsi[4] = weightsKsi[7] = 8.0 / 9.0;
		}
	}

	void calcdNdksideta()
	{
		for (int i = 0; i < this->numberPc; i++)
		{
			this->dNdksi[i][0] = N1ksi(this->eta[i]);
			this->dNdksi[i][1] = N2ksi(this->eta[i]);
			this->dNdksi[i][2] = N3ksi(this->eta[i]);
			this->dNdksi[i][3] = N4ksi(this->eta[i]);

			this->dNdeta[i][0] = N1eta(this->ksi[i]);
			this->dNdeta[i][1] = N2eta(this->ksi[i]);
			this->dNdeta[i][2] = N3eta(this->ksi[i]);
			this->dNdeta[i][3] = N4eta(this->ksi[i]);
		}
	}


	void displayElementdksideta()
	{
		cout << "Element\n";
		cout << "\t eta \t\t dN1/dksi \t dN2/dksi \t dN3/dksi \t dN4/dksi \n";
		for (int i = 0; i < this->numberPc; i++)
		{
			cout << "PC" << i + 1 << " \t " << this->eta[i] << " \t";
			for (int j = 0; j < 4; j++)
			{
				cout << this->dNdksi[i][j] << "\t";
			}
			cout << endl;
		}

		cout << endl << "\t ksi \t\t dN1/deta \t dN2/deta \t dN3/deta \t dN4/deta \n";
		for (int i = 0; i < this->numberPc; i++)
		{
			cout << "PC" << i + 1 << " \t " << this->ksi[i] << " \t";
			for (int j = 0; j < 4; j++)
			{
				cout << this->dNdeta[i][j] << "\t";
			}
			cout << endl;
		}
	}

	void calcNdxNdy(int nrPc, double matrix[2][2])
	{
		for (int nrdN = 0; nrdN < 4; nrdN++)
		{
			this->dNdx[nrPc][nrdN] = matrix[0][0] * this->dNdksi[nrPc][nrdN] + matrix[0][1] * this->dNdeta[nrPc][nrdN];
			this->dNdy[nrPc][nrdN] = matrix[1][0] * this->dNdksi[nrPc][nrdN] + matrix[1][1] * this->dNdeta[nrPc][nrdN];
		}
	}

	void calcJacobianMatrixcalcdNdxdNdy(double coords[2][4])
	{
		double dydeta = 0.0, dydksi = 0.0, dxdksi = 0.0, dxdeta = 0.0;
		for (int nrPc = 0; nrPc < this->numberPc; nrPc++)
		{
			dydeta = 0.0, dydksi = 0.0, dxdksi = 0.0, dxdeta = 0.0;
			for (int nrdN = 0; nrdN < 4; nrdN++)
			{

				dydeta += coords[1][nrdN] * this->dNdeta[nrPc][nrdN];
				dydksi += coords[1][nrdN] * this->dNdksi[nrPc][nrdN];
				dxdeta += coords[0][nrdN] * this->dNdeta[nrPc][nrdN];
				dxdksi += coords[0][nrdN] * this->dNdksi[nrPc][nrdN];
			}
			detJ[nrPc] = dxdksi * dydeta - (dxdeta * dydksi);
			double counter = 1.0 / detJ[nrPc];
			double matrix[2][2] = { {(dydeta * counter), (-dydksi * counter)}, {(-dxdeta * counter), (dxdksi * counter)} };
			 
				calcNdxNdy(nrPc, matrix);
			}
	}

	void displayElementdxdy()
	{
		cout << "\t\t dN1/dx \t dN2/dx \t dN3/dx \t dN4/dx \n";
		for (int i = 0; i < this->numberPc; i++)
		{
			cout << "PC" << i + 1 << " \t ";
			for (int j = 0; j < 4; j++)
			{
				cout << this->dNdx[i][j] << "\t";
			}
			cout << endl;
		}

		cout << endl << "\t\t dN1/dy \t dN2/dy \t dN3/dy \t dN4/dy \n";
		for (int i = 0; i < this->numberPc; i++)
		{
			cout << "PC" << i + 1 << " \t " << " \t";
			for (int j = 0; j < 4; j++)
			{
				cout << this->dNdy[i][j] << "\t";
			}
			cout << endl;
		}
	}
	int getNumberPc()
	{
		return this->numberPc;
	}

	void setMatrixHPcValue(int nrPc, double** matrix)
	{
		for (int row = 0; row < 4; row++)
			for (int col = 0; col < 4; col++)
				this->matrixHPc[nrPc][row][col] = matrix[row][col];
	}

	void calcMatrixHPc(int nrPc, int k)
	{
		double** Nx, **Ny, **tmp;
		Nx = new double* [4];
		Ny = new double* [4];
		tmp = new double* [4];
		for (int i = 0; i < 4; i++)
		{
			Nx[i] = new double[4];
			Ny[i] = new double[4];
			tmp[i] = new double[4];
		}
		Nx = multiplicationMatrix(dNdx[nrPc], dNdx[nrPc]);
		Nx = transpositionMatrix(Nx);
		Ny = multiplicationMatrix(dNdy[nrPc], dNdy[nrPc]);
		Ny = transpositionMatrix(Ny);
		//displayMatrix(Nx, 4, 4);
		//cout << endl;
		//displayMatrix(Ny, 4, 4);
		tmp = multplicationMatrixByScalar(sumMatrix(Nx, Ny, 4, 4), k * detJ[nrPc], 4, 4);
		setMatrixHPcValue(nrPc, tmp);
		//displayMatrix(this->matrixHPc[nrPc], 4, 4);
		
		for (int i = 0; i < 4; i++)
		{
			tmp[i] = nullptr;
			delete tmp[i];
			Nx[i] = nullptr;
			delete Nx[i];
			Ny[i] = nullptr;
			delete Ny[i];
		}
		tmp = nullptr;
		Nx = nullptr;
		Ny = nullptr;
		delete tmp;
		delete Nx;
		delete Ny;
	}

	void displayMatrixHPc(int nrPc)
	{
		displayMatrix(this->matrixHPc[nrPc], 4, 4);
	}

	void displayDetJ()
	{
		cout << "\nDetJ\n";
		for (int i = 0; i < numberPc; i++)
			cout << detJ[i] << " ";
		cout << endl;
	}

	void calcMatrixH()
	{
		for (int row = 0; row < 4; row++)
			for (int col = 0; col < 4; col++)
				matrixH[row][col] = matrixHPc[0][row][col];
		
		for (int i = 1; i < numberPc; i++)
		{
			multplicationMatrixByScalar(matrixHPc[i], weightsEta[i] * weightsKsi[i], 4, 4);
			sumMatrix(matrixH, matrixHPc[i], 4, 4);
		}
		displayMatrix(matrixH, 4, 4);
	}

	~Element4()
	{
		eta = nullptr;
		ksi = nullptr;
		for (int i = 0; i < numberPc; i++)
		{
			matrixH[i] = nullptr;
			delete matrixH[i];
			dNdksi[i] = nullptr;
			delete dNdksi[i];
			dNdeta[i] = nullptr;
			delete dNdeta[i];
			dNdx[i] = nullptr;
			delete dNdx[i];
			dNdy[i] = nullptr;
			delete dNdy[i];
			for (int j = 0; j < numberPc; j++)
			{
				matrixHPc[i][j] = nullptr;
				delete matrixHPc[i][j];
			}
			matrixHPc[i] = nullptr;
			delete matrixHPc[i];
		
		}
		detJ = nullptr;
		weightsEta = nullptr;
		delete weightsEta;
		weightsKsi = nullptr;
		delete weightsKsi;
		delete eta;
		delete ksi;
		delete detJ;
		delete dNdy;
		delete dNdx;
		delete dNdksi;
		delete dNdeta;
		delete matrixHPc;
		delete matrixH;
	}
};

struct Element2pc
{
	double eta[4] = { -1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
	double ksi[4] = { -1.0 / sqrt(3.0), 1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
	double dNdksi[4][4];
	double dNdeta[4][4];
};

struct Element3pc
{
	double eta[9] = { -sqrt(3.0/ 5.0), -sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), 0, 0, 0, sqrt(3.0 / 5.0), sqrt(3.0 / 5.0), sqrt(3.0 / 5.0) };
	double ksi[9] = { -sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0) };
	double dNdksi[9][4];
	double dNdeta[9][4];
};

Element2pc calcElement2pc()
{
	
	Element2pc element;
	for (int i = 0; i < 4; i++)
	{
		element.dNdksi[i][0] = N1ksi(element.eta[i]);
		element.dNdksi[i][1] = N2ksi(element.eta[i]);
		element.dNdksi[i][2] = N3ksi(element.eta[i]);
		element.dNdksi[i][3] = N4ksi(element.eta[i]);

		element.dNdeta[i][0] = N1eta(element.ksi[i]);
		element.dNdeta[i][1] = N2eta(element.ksi[i]);
		element.dNdeta[i][2] = N3eta(element.ksi[i]);
		element.dNdeta[i][3] = N4eta(element.ksi[i]);
	}
	return element;
}

Element3pc calcElement3pc()
{

	Element3pc element;
	for (int i = 0; i < 9; i++)
	{
		element.dNdksi[i][0] = N1ksi(element.eta[i]);
		element.dNdksi[i][1] = N2ksi(element.eta[i]);
		element.dNdksi[i][2] = N3ksi(element.eta[i]);
		element.dNdksi[i][3] = N4ksi(element.eta[i]);

		element.dNdeta[i][0] = N1eta(element.ksi[i]);
		element.dNdeta[i][1] = N2eta(element.ksi[i]);
		element.dNdeta[i][2] = N3eta(element.ksi[i]);
		element.dNdeta[i][3] = N4eta(element.ksi[i]);
		
	}
	return element;
}



void displayElement2pc(Element2pc element)
{
	cout << "Element2pc\n";
	cout << "\t eta \t\t dN1/dksi \t dN2/dksi \t dN3/dksi \t dN4/dksi \n";
	for (int i = 0; i < 4; i++)
	{
		cout << "PC" << i + 1 << " \t " << element.eta[i] << " \t";
		for (int j = 0; j < 4; j++)
		{
			cout << element.dNdksi[i][j] << "\t";
		}
		cout << endl;
	}

	cout << endl << "\t ksi \t\t dN1/deta \t dN2/deta \t dN3/deta \t dN4/deta \n";
	for (int i = 0; i < 4; i++)
	{
		cout << "PC" << i + 1 << " \t " << element.ksi[i] << " \t";
		for (int j = 0; j < 4; j++)
		{
			cout << element.dNdeta[i][j] << "\t";
		}
		cout << endl;
	}
}

void displayElement3pc(Element3pc element)
{
	cout << "Element3pc\n";
	cout << "\t eta \t\t dN1/dksi \t dN2/dksi \t dN3/dksi \t dN4/dksi \n";
	for (int i = 0; i < 9; i++)
	{
		cout << "PC" << i + 1 << " \t " << element.eta[i] << " \t";
		for (int j = 0; j < 4; j++)
		{
			cout << element.dNdksi[i][j] << "\t";
		}
		cout << endl;
	}

	cout << endl << "\t ksi \t\t dN1/deta \t dN2/deta \t dN3/deta \t dN4/deta \n";
	for (int i = 0; i < 9; i++)
	{
		cout << "PC" << i + 1 << " \t " << element.ksi[i] << " \t";
		for (int j = 0; j < 4; j++)
		{
			cout << element.dNdeta[i][j] << "\t";
		}
		cout << endl;
	}
}

double** calcSingleJacobianMatrix2pc(double coords[2][4], Element2pc element, int pc)
{
	double dydeta = 0.0, dydksi = 0.0, dxdksi = 0.0, dxdeta = 0.0, detJ;
	for (int i = 0; i < 4; i++)
	{
		
		dydeta += coords[1][i] * element.dNdeta[pc][i];
		dydksi += coords[1][i] * element.dNdksi[pc][i];
		dxdeta += coords[0][i] * element.dNdeta[pc][i];
		dxdksi += coords[0][i] * element.dNdksi[pc][i];
	}
	detJ = dxdksi * dydeta - (dxdeta * dydksi);
	double counter = 1.0 / detJ;
	double matrix[2][2] = { {(dydeta * counter), (-dydksi * counter)}, {(-dxdeta * counter), (dxdksi * counter)}};
	cout << matrix[0][0] << " " << matrix[0][1] << endl << matrix[1][0] << " " << matrix[1][1] << endl;
	double** resPtr = new double*[2];
	resPtr[0] = new double[4];
	resPtr[1] = new double[4];
	for (int i = 0; i < 4; i++)
	{
		resPtr[0][i] = matrix[0][0] * element.dNdksi[pc][i] + matrix[0][1] * element.dNdeta[pc][i];
		resPtr[1][i] = matrix[1][0] * element.dNdksi[pc][i] + matrix[1][1] * element.dNdeta[pc][i];
	}
	
	
	return resPtr;
}

//int main()
//{
//	double coords[2][4] = { {0, 0.025, 0.025, 0}, {0, 0, 0.025, 0.025} };
//	//double coords[2][4] = { {0.0, 9.0, 8.0, 0.0}, {0.0, 0.0, 4.0, 4.0} };
//	Element4 init(3);
//	init.calcdNdksideta();
//	init.displayElementdksideta();
//	init.calcJacobianMatrixcalcdNdxdNdy(coords);
//	init.displayElementdxdy();
//	for (int nrPc = 0; nrPc < init.getNumberPc(); nrPc++)
//	{
//		cout << "Matrix H nr " << nrPc + 1 << endl;
//		init.calcMatrixHPc(nrPc, 30.0);
//		
//		init.displayMatrixHPc(nrPc);
//	}
//	cout << "\nMatrix H\n";
//	init.calcMatrixH();
//	init.displayDetJ();
//	/*Element2pc elem2pc = calcElement2pc();
//	displayElement2pc(elem2pc);
//	Element3pc elem3pc = calcElement3pc();
//	displayElement3pc(elem3pc);
//	double coords[2][4] = { {0, 0.025, 0.025, 0}, {0, 0, 0.025, 0.025} };
//	double resPc2dNdx[4][4];
//	double resPc2dNdy[4][4];
//	double** dNdxdNdy;
//	for (int i = 0; i < 4; i++)
//	{
//		dNdxdNdy = calcSingleJacobianMatrix2pc(coords, elem2pc, i);
//		for (int j = 0; j < 4; j++)
//		{
//			resPc2dNdx[i][j] = dNdxdNdy[0][j];
//			resPc2dNdy[i][j] = dNdxdNdy[1][j];
//		}
//		
//	}
//
//	cout << "Element2pc\n";
//	cout << "\t\t dN1/dx \t dN2/dx \t dN3/dx \t dN4/dx \n";
//	for (int i = 0; i < 4; i++)
//	{
//		cout << "PC" << i + 1 << " \t ";
//		for (int j = 0; j < 4; j++)
//		{
//			cout << resPc2dNdx[i][j] << "\t";
//		}
//		cout << endl;
//	}
//
//	cout << endl << "\t\t dN1/dy \t dN2/dy \t dN3/dy \t dN4/dy \n";
//	for (int i = 0; i < 4; i++)
//	{
//		cout << "PC" << i + 1 << " \t " << " \t";
//		for (int j = 0; j < 4; j++)
//		{
//			cout << resPc2dNdy[i][j] << "\t";
//		}
//		cout << endl;
//	}
//	*/
//	return 0;
//}

double** multiplicationMatrix(double* verMatrix, double* horMatrix)
{
	double** resMatrix = new double* [4];
	for (int i = 0; i < 4; i++)
	{
		resMatrix[i] = new double[4];
	}
	cout << endl;
	for (int row = 0; row < 4; row++)
		for (int col = 0; col < 4; col++)
			resMatrix[row][col] = verMatrix[row] * horMatrix[col];
	return resMatrix;
	resMatrix = nullptr;
	delete resMatrix;
}

double** transpositionMatrix(double** matrix)
{
	double tmp;
	for (int row = 0; row < 4; row++)
		for (int col = 0; col < 4; col++)
		{
			tmp = matrix[row][col];
			matrix[row][col] = matrix[col][row];
			matrix[col][row] = tmp;
		}
	return matrix;
}

void displayMatrix(double** matrix, int row, int col)
{
	cout << endl;
	for (int r = 0; r < row; r++)
	{
		for (int c = 0; c < col; c++)
			cout << matrix[r][c] << "\t";
		cout << endl;
	}
	cout << endl;
}

double** sumMatrix(double** matrix1, double** matrix2, int row, int col)
{
	for (int r = 0; r < row; r++)
		for (int c = 0; c < col; c++)
			matrix1[r][c] += matrix2[r][c];
	return matrix1;
}

double** multplicationMatrixByScalar(double** matrix, double scalar, int row, int col)
{
	for (int r = 0; r < row; r++)
		for (int c = 0; c < col; c++)
			matrix[r][c] *= scalar;
	return matrix;
}
