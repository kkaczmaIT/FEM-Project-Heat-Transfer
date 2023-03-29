#include "MainSystem.h"


MainSystem::MainSystem(int numberNodes)
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
	N = new double* [this->numberPc];
	matrixHPc = new double** [this->numberPc];
	matrixCPc = new double** [this->numberPc];
	matrixH = new double* [this->numberPc];
	matrixC = new double* [this->numberPc];
	this->vectorP = new double* [4];
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			vectorP[i] = new double[4];
		}
	}


	for (int i = 0; i < this->numberPc; i++)
	{
		dNdksi[i] = new double[4];
		dNdeta[i] = new double[4];
		dNdy[i] = new double[4];
		dNdx[i] = new double[4];
		matrixHPc[i] = new double* [4];
		matrixH[i] = new double[4];
		matrixC[i] = new double[4];
		N[i] = new double[4];
		matrixCPc[i] = new double* [4];

		for (int j = 0; j < this->numberPc; j++)
		{
			matrixHPc[i][j] = new double[4];
			matrixCPc[i][j] = new double[4];
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
		this->eta[3] = 0.0;
		this->eta[4] = 0.0;
		this->eta[5] = 0.0;
		this->eta[6] = sqrt(3.0 / 5.0);
		this->eta[7] = sqrt(3.0 / 5.0);
		this->eta[8] = sqrt(3.0 / 5.0);



		this->ksi[0] = -sqrt(3.0 / 5.0);
		this->ksi[1] = 0.0;
		this->ksi[2] = sqrt(3.0 / 5.0);
		this->ksi[3] = -sqrt(3.0 / 5.0);
		this->ksi[4] = 0.0;
		this->ksi[5] = sqrt(3.0 / 5.0);
		this->ksi[6] = -sqrt(3.0 / 5.0);
		this->ksi[7] = 0.0;
		this->ksi[8] = sqrt(3.0 / 5.0);

		weightsEta[0] = weightsEta[1] = weightsEta[2] = weightsEta[6] = weightsEta[7] = weightsEta[8] = 5.0 / 9.0;

		weightsEta[3] = weightsEta[4] = weightsEta[5] = 8.0 / 9.0;

		weightsKsi[0] = weightsKsi[3] = weightsKsi[6] = weightsKsi[2] = weightsKsi[5] = weightsKsi[8] = 5.0 / 9.0;

		weightsKsi[1] = weightsKsi[4] = weightsKsi[7] = 8.0 / 9.0;
	}

	else if (this->numberNodes == 4)
	{
		this->eta[0] = -0.861136;
		this->eta[1] = -0.861136;
		this->eta[2] = -0.861136;
		this->eta[3] = -0.861136;;
		this->eta[4] = -0.339981;
		this->eta[5] = -0.339981;
		this->eta[6] = -0.339981;
		this->eta[7] = -0.339981;
		this->eta[8] = 0.339981;
		this->eta[9] = 0.339981;
		this->eta[10] = 0.339981;
		this->eta[11] = 0.339981;
		this->eta[12] = 0.861136;
		this->eta[13] = 0.861136;
		this->eta[14] = 0.861136;
		this->eta[15] = 0.861136;


		this->ksi[0] = -0.861136;
		this->ksi[1] = -0.339981;
		this->ksi[2] = 0.339981;
		this->ksi[3] = 0.861136;
		this->ksi[4] = -0.861136;
		this->ksi[5] = -0.339981;
		this->ksi[6] = 0.339981;
		this->ksi[7] = 0.861136;
		this->ksi[8] = -0.861136;
		this->ksi[9] = -0.339981;
		this->ksi[10] = 0.339981;
		this->ksi[11] = 0.861136;
		this->ksi[12] = -0.861136;
		this->ksi[13] = -0.339981;
		this->ksi[14] = 0.339981;
		this->ksi[15] = 0.861136;

		weightsEta[0] = weightsEta[1] = weightsEta[2] = weightsEta[3] = weightsEta[12] = weightsEta[13] = weightsEta[14] = weightsEta[15] = 0.347855;

		weightsEta[4] = weightsEta[5] = weightsEta[6] = weightsEta[7] = weightsEta[8] = weightsEta[9] = weightsEta[10] = weightsEta[11] = 0.652145;

		weightsKsi[0] = weightsKsi[4] = weightsKsi[8] = weightsKsi[12] = weightsKsi[3] = weightsKsi[7] = weightsKsi[11] = weightsKsi[15] = 0.347855;

		weightsKsi[1] = weightsKsi[2] = weightsKsi[5] = weightsKsi[6] = weightsKsi[9] = weightsKsi[10] = weightsKsi[13] = weightsKsi[14] = 0.652145;
		}
	}


	double** MainSystem::getVectorP()
	{
	
		return this->vectorP;
	}

	double** MainSystem::getMatrixH()
	{
		return this->matrixH;
	}

	void MainSystem::calcN()
	{
		for (int i = 0; i < this->numberPc; i++)
		{
			this->N[i][0] = N1(this->ksi[i], this->eta[i]);
			this->N[i][1] = N2(this->ksi[i], this->eta[i]);
			this->N[i][2] = N3(this->ksi[i], this->eta[i]);
			this->N[i][3] = N4(this->ksi[i], this->eta[i]);
		}
	}

	void MainSystem::displayN()
	{
		std::cout << "N1 N2 N3 N4\n";
		for (int i = 0; i < this->numberPc; i++)
		{
			std::cout << this->N[i][0] << " " << N[i][1] << " " << N[i][2] << " " << N[i][3] << std::endl;
		}
	}

	void MainSystem::calcdNdksideta()
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

	double** MainSystem::calcMatrixHBc(double coords[2][4], bool isBc[4], double alpha, double tot)
	{
		double*** matrixHBcPc;
		double* etaLocal, * ksiLocal, * weightsEtaLocal, * weightsKsiLocal;
		double detJLocal[4];
		bool isFirst = true;
		matrixHBcPc = new double** [4];
		for (int i = 0; i < 4; i++)
		{
			matrixHBcPc[i] = new double* [this->numberNodes];
			for (int j = 0; j < this->numberNodes + 1; j++)
			{
				matrixHBcPc[i][j] = new double[4];
			}
		}
		etaLocal = new double[4 * this->numberNodes] {0};
		weightsEtaLocal = new double[4 * this->numberNodes] {0};
		ksiLocal = new double[4 * this->numberNodes] {0};
		weightsKsiLocal = new double[4 * this->numberNodes] {0};

		if (this->numberNodes == 2)
		{
			weightsEtaLocal[0] = weightsEtaLocal[1] = weightsEtaLocal[2] = weightsEtaLocal[3] = weightsEtaLocal[4] = weightsEtaLocal[5] = weightsEtaLocal[6] = weightsEtaLocal[7] = weightsKsiLocal[0] = weightsKsiLocal[1] = weightsKsiLocal[2] = weightsKsiLocal[3] = weightsKsiLocal[4] = weightsKsiLocal[5] = weightsKsiLocal[6] = weightsKsiLocal[7] = 1.0;
			etaLocal[0] = -1.0;
			etaLocal[1] = -1.0;
			etaLocal[2] = -1.0 / sqrt(3.0);
			etaLocal[3] = 1.0 / sqrt(3.0);
			etaLocal[4] = 1.0;
			etaLocal[5] = 1.0;
			etaLocal[6] = 1.0 / sqrt(3.0);
			etaLocal[7] = -1.0 / sqrt(3.0);

			ksiLocal[0] = -1.0 / sqrt(3.0);
			ksiLocal[1] = 1.0 / sqrt(3.0);
			ksiLocal[2] = 1.0;
			ksiLocal[3] = 1.0;
			ksiLocal[4] = 1.0 / sqrt(3.0);
			ksiLocal[5] = -1.0 / sqrt(3.0);
			ksiLocal[6] = -1.0;
			ksiLocal[7] = -1.0;
		}
		else if (this->numberNodes == 3)
		{
			etaLocal[0] = -1.0;//sqrt(3.0 / 5.0);
			etaLocal[1] = -1.0;//sqrt(3.0 / 5.0);
			etaLocal[2] = -1.0;//sqrt(3.0 / 5.0);
			etaLocal[3] = -sqrt(3.0 / 5.0);//0;
			etaLocal[4] = 0;
			etaLocal[5] = sqrt(3.0 / 5.0);//0;
			etaLocal[6] = 1.0;//sqrt(3.0 / 5.0);
			etaLocal[7] = 1.0;//sqrt(3.0 / 5.0);
			etaLocal[8] = 1.0;//sqrt(3.0 / 5.0);
			etaLocal[9] = sqrt(3.0 / 5.0);//0;
			etaLocal[10] = 0;// sqrt(3.0 / 5.0);
			etaLocal[11] = -sqrt(3.0 / 5.0);



			ksiLocal[0] = -sqrt(3.0 / 5.0);
			ksiLocal[1] = 0;
			ksiLocal[2] = sqrt(3.0 / 5.0);
			ksiLocal[3] = 1.0;// -sqrt(3.0 / 5.0);
			ksiLocal[4] = 1.0;// 0;
			ksiLocal[5] = 1.0;//sqrt(3.0 / 5.0);
			ksiLocal[6] = sqrt(3.0 / 5.0);
			ksiLocal[7] = 0;
			ksiLocal[8] = -sqrt(3.0 / 5.0);
			ksiLocal[9] = -1.0;// -sqrt(3.0 / 5.0);
			ksiLocal[10] = -1.0;// -sqrt(3.0 / 5.0);
			ksiLocal[11] = -1.0;//-sqrt(3.0 / 5.0);

			weightsEtaLocal[0] = weightsEtaLocal[1] = weightsEtaLocal[2] = weightsEtaLocal[3] = weightsEtaLocal[5] = weightsEtaLocal[6] = weightsEtaLocal[7] = weightsEtaLocal[8] = weightsEtaLocal[9] = weightsEtaLocal[11] = 5.0 / 9.0;

			weightsEtaLocal[4] = weightsEtaLocal[10] = 8.0 / 9.0;

			weightsKsiLocal[0] = weightsKsiLocal[2] = weightsKsiLocal[3] = weightsKsiLocal[4] = weightsKsiLocal[5] = weightsKsiLocal[6] = weightsKsiLocal[8] = weightsKsiLocal[9] = weightsKsiLocal[10] = weightsKsiLocal[11] = 5.0 / 9.0;

			weightsKsiLocal[1] = weightsKsiLocal[7] = 8.0 / 9.0;
		}

		else if (this->numberNodes == 4)
		{
			etaLocal[0] = -1.0;//-0.861136;
			etaLocal[1] = -1.0;//-0.861136;
			etaLocal[2] = -1.0;// -0.861136;
			etaLocal[3] = -1.0;// -0.861136;;
			etaLocal[4] = -0.861136;//-0.339981;
			etaLocal[5] = -0.339981;
			etaLocal[6] = 0.339981;//-0.339981;
			etaLocal[7] = 0.861136;//-0.339981;
			etaLocal[8] = 1.0;//0.339981;
			etaLocal[9] = 1.0;//0.339981;
			etaLocal[10] = 1.0;//0.339981;
			etaLocal[11] = 1.0;//0.339981;
			etaLocal[12] = 0.861136;
			etaLocal[13] = 0.339981;//0.861136;
			etaLocal[14] = -0.339981;//0.861136;
			etaLocal[15] = -0.861136;//0.861136;


			ksiLocal[0] = -0.861136;
			ksiLocal[1] = -0.339981;
			ksiLocal[2] = 0.339981;
			ksiLocal[3] = 0.861136;
			ksiLocal[4] = 1.0;//-0.861136;
			ksiLocal[5] = 1.0;// -0.339981;
			ksiLocal[6] = 1.0;//0.339981;
			ksiLocal[7] = 1.0;//0.861136;
			ksiLocal[8] = 0.861136;//-0.861136;
			ksiLocal[9] = 0.339981;//-0.339981;
			ksiLocal[10] = -0.339981;//0.339981;
			ksiLocal[11] = -0.861136;//0.861136;
			ksiLocal[12] = -1.0;// -0.861136;
			ksiLocal[13] = -1.0;//-0.339981;
			ksiLocal[14] = -1.0;//0.339981;
			ksiLocal[15] = -1.0;// 0.861136;



			weightsEtaLocal[0] = weightsEtaLocal[3] = weightsEtaLocal[4] = weightsEtaLocal[7] = weightsEtaLocal[8] = weightsEtaLocal[11] = weightsEtaLocal[12] = weightsEtaLocal[15] = 0.347855;

			weightsEtaLocal[1] = weightsEtaLocal[2] = weightsEtaLocal[5] = weightsEtaLocal[6] = weightsEtaLocal[9] = weightsEtaLocal[10] = weightsEtaLocal[13] = weightsEtaLocal[14] = 0.652145;

			weightsKsiLocal[0] = weightsKsiLocal[3] = weightsKsiLocal[4] = weightsKsiLocal[7] = weightsKsiLocal[8] = weightsKsiLocal[11] = weightsKsiLocal[12] = weightsKsiLocal[15] = 0.347855;

			weightsKsiLocal[1] = weightsKsiLocal[2] = weightsKsiLocal[5] = weightsKsiLocal[6] = weightsKsiLocal[9] = weightsKsiLocal[10] = weightsKsiLocal[13] = weightsKsiLocal[14] = 0.652145;

		}

		for (int side = 0; side < 4; side++)
		{
			for (int nrPc = 0, index = side * this->numberNodes; nrPc < this->numberNodes; nrPc++, index++)
			{
				matrixHBcPc[side][nrPc][0] = N1(ksiLocal[index], etaLocal[index]);
				matrixHBcPc[side][nrPc][1] = N2(ksiLocal[index], etaLocal[index]);
				matrixHBcPc[side][nrPc][2] = N3(ksiLocal[index], etaLocal[index]);
				matrixHBcPc[side][nrPc][3] = N4(ksiLocal[index], etaLocal[index]);
			}
			detJLocal[side] = sqrt(pow(coords[0][side] - coords[0][(side + 1) % 4], 2) + pow(coords[1][side] - coords[1][(side + 1) % 4], 2)) / 2;
		}
		double** matrixTemp, ** matrixHBTemp = NULL, ** matrixRes = NULL, * vectorPPart = NULL;

		for (int side = 0; side < 4; side++)
		{
			matrixTemp = NULL;
			matrixHBTemp = NULL;
			vectorPPart = NULL;
			for (int nrPc = 0; nrPc < this->numberNodes; nrPc++)
			{
				matrixTemp = multiplicationMatrix(matrixHBcPc[side][nrPc], matrixHBcPc[side][nrPc]);
				// wektor P
				if (side == 0 || side == 2)
					vectorPPart = multplicationVectorByScalar(matrixHBcPc[side][nrPc], tot * weightsKsiLocal[side * this->numberNodes + nrPc], 4);
				else if (side == 1 || side == 3)
					vectorPPart = multplicationVectorByScalar(matrixHBcPc[side][nrPc], tot * weightsEtaLocal[side * this->numberNodes + nrPc], 4);
				if (nrPc == 0)
					setVectorValuesByVector(vectorP[side], vectorPPart, 4);
				else
					setVectorValuesByVector(vectorP[side], sumVector(vectorP[side], vectorPPart, 4), 4);
				if (side == 0 || side == 2)
					matrixTemp = multplicationMatrixByScalar(matrixTemp, /*weightsEtaLocal[side * this->numberNodes + nrPc] * */ weightsKsiLocal[side * this->numberNodes + nrPc], 4, 4);
				else if (side == 1 || side == 3)
					matrixTemp = multplicationMatrixByScalar(matrixTemp, weightsEtaLocal[side * this->numberNodes + nrPc] /* * weightsKsiLocal[side * this->numberNodes + nrPc]*/, 4, 4);

				if (nrPc == 0)
					matrixHBTemp = matrixTemp;
				else if (nrPc != 0)
					matrixHBTemp = sumMatrix(matrixHBTemp, matrixTemp, 4, 4);
			}
			matrixHBTemp = multplicationMatrixByScalar(matrixHBTemp, alpha * detJLocal[side], 4, 4);
			setVectorValuesByVector(vectorP[side], multplicationVectorByScalar(vectorP[side], detJLocal[side] * alpha, 4), 4);
			if (isBc[side] && isBc[(side + 1) % 4] && isFirst)
			{
				matrixRes = matrixHBTemp;
				isFirst = false;
			}
			else if (isBc[side] && isBc[(side + 1) % 4])
				matrixRes = sumMatrix(matrixRes, matrixHBTemp, 4, 4);
		}

		for (int i = 0; i < 4; i++)
		{

			for (int j = 0; j < this->numberNodes; j++)
			{
				matrixHBcPc[i][j] = nullptr;
				delete matrixHBcPc[i][j];
			}
			matrixHBcPc[i] = nullptr;
			delete matrixHBcPc[i];
		}

		matrixHBcPc = nullptr;
		delete matrixHBcPc;

		etaLocal = nullptr;
		delete[] etaLocal;
		ksiLocal = nullptr;
		delete[] ksiLocal;
		weightsKsiLocal = nullptr;
		delete[] weightsKsiLocal;
		weightsEtaLocal = nullptr;
		delete[] weightsEtaLocal;
		vectorPPart = nullptr;
		delete vectorPPart;

		return matrixRes;
	}

	void MainSystem::calcMatrixCPc(int nrPc, double c, double p)
	{
		double** matrixCTemp = NULL;
		matrixCTemp = multiplicationMatrix(N[nrPc], N[nrPc]);
		matrixCTemp = multplicationMatrixByScalar(matrixCTemp, c * p * weightsEta[nrPc] * weightsKsi[nrPc] * detJ[nrPc], 4, 4);
		setMatrixCPcValue(nrPc, matrixCTemp);
	}

	void MainSystem::setMatrixCPcValue(int nrPc, double** matrixCTemp)
	{
		for (int row = 0; row < 4; row++)
		{
			for (int column = 0; column < 4; column++)
				matrixCPc[nrPc][row][column] = matrixCTemp[row][column];
		}
	}

	void MainSystem::calcMatrixC()
	{
		for (int row = 0; row < 4; row++)
		{
			for (int column = 0; column < 4; column++)
			{
				matrixC[row][column] = matrixCPc[0][row][column];
			}
		}

		for (int nrPc = 1; nrPc < this->numberPc; nrPc++)
		{

			sumMatrix(matrixC, matrixCPc[nrPc], 4, 4);
		}


	}

	double** MainSystem::getMatrixC()
	{
		return this->matrixC;
	}


	/*Test function weights*/
	void MainSystem::displayMultplicationWeight()
	{
		std::cout << "Iloczyn wag\n";
		for (int i = 0; i < this->numberPc; i++)
		{
			std::cout << this->weightsEta[i] * this->weightsKsi[i] << " ";
		}
		std::cout << std::endl;
	}
	/*Test function weights*/

	void MainSystem::displayElementdksideta()
	{
		std::cout << "Element\n";
		std::cout << "\t eta \t\t dN1/dksi \t dN2/dksi \t dN3/dksi \t dN4/dksi \n";
		for (int i = 0; i < this->numberPc; i++)
		{
			std::cout << "PC" << i + 1 << " \t " << this->eta[i] << " \t";
			for (int j = 0; j < 4; j++)
			{
				std::cout << this->dNdksi[i][j] << "\t";
			}
			std::cout << std::endl;
		}

		std::cout << std::endl << "\t ksi \t\t dN1/deta \t dN2/deta \t dN3/deta \t dN4/deta \n";
		for (int i = 0; i < this->numberPc; i++)
		{
			std::cout << "PC" << i + 1 << " \t " << this->ksi[i] << " \t";
			for (int j = 0; j < 4; j++)
			{
				std::cout << this->dNdeta[i][j] << "\t";
			}
			std::cout << std::endl;
		}
	}

	void MainSystem::calcNdxNdy(int nrPc, double matrix[2][2])
	{
		for (int nrdN = 0; nrdN < 4; nrdN++)
		{
			this->dNdx[nrPc][nrdN] = matrix[0][0] * this->dNdksi[nrPc][nrdN] + matrix[0][1] * this->dNdeta[nrPc][nrdN];
			this->dNdy[nrPc][nrdN] = matrix[1][0] * this->dNdksi[nrPc][nrdN] + matrix[1][1] * this->dNdeta[nrPc][nrdN];
		}
	}

	void MainSystem::calcJacobianMatrixcalcdNdxdNdy(double coords[2][4])
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

	void MainSystem::displayElementdxdy()
	{
		std::cout << "\t\t dN1/dx \t dN2/dx \t dN3/dx \t dN4/dx \n";
		for (int i = 0; i < this->numberPc; i++)
		{
			std::cout << "PC" << i + 1 << " \t ";
			for (int j = 0; j < 4; j++)
			{
				std::cout << this->dNdx[i][j] << "\t";
			}
			std::cout << std::endl;
		}

		std::cout << std::endl << "\t\t dN1/dy \t dN2/dy \t dN3/dy \t dN4/dy \n";
		for (int i = 0; i < this->numberPc; i++)
		{
			std::cout << "PC" << i + 1 << " \t " << " \t";
			for (int j = 0; j < 4; j++)
			{
				std::cout << this->dNdy[i][j] << "\t";
			}
			std::cout << std::endl;
		}
	}
	int MainSystem::getNumberPc()
	{
		return this->numberPc;
	}

	void MainSystem::setMatrixHPcValue(int nrPc, double** matrix)
	{
		for (int row = 0; row < 4; row++)
			for (int col = 0; col < 4; col++)
				this->matrixHPc[nrPc][row][col] = matrix[row][col];
	}

	void MainSystem::displayPCElements()
	{
		for (int i = 0; i < this->numberPc; i++)
		{
			std::cout << "nrPc " << i + 1 << " eta = " << this->eta[i] << " ksi = " << this->ksi[i] << " weightEta = " << this->weightsEta[i] << " weightKsi = " << this->weightsKsi[i] << std::endl;
		}
	}

	void MainSystem::calcMatrixHPc(int nrPc, int k)
	{
		double** Nx, ** Ny, ** tmp;
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
		tmp = multplicationMatrixByScalar(sumMatrix(Nx, Ny, 4, 4), k * detJ[nrPc] * weightsEta[nrPc] * weightsKsi[nrPc], 4, 4);

		setMatrixHPcValue(nrPc, tmp);

		for (int i = 0; i < 4; i++)
		{
			tmp[i] = nullptr;
			delete[] tmp[i];
			Nx[i] = nullptr;
			delete[] Nx[i];
			Ny[i] = nullptr;
			delete[] Ny[i];
		}
		tmp = nullptr;
		Nx = nullptr;
		Ny = nullptr;
		delete tmp;
		delete Nx;
		delete Ny;
	}



	void MainSystem::displayMatrixHPc(int nrPc)
	{
		displayMatrix(this->matrixHPc[nrPc], 4, 4);
	}

	void MainSystem::displayDetJ()
	{
		std::cout << "\nDetJ\n";
		for (int i = 0; i < numberPc; i++)
			std::cout << detJ[i] << " ";
		std::cout << std::endl;
	}

	void MainSystem::calcMatrixH()
	{
		for (int row = 0; row < 4; row++)
			for (int col = 0; col < 4; col++)
				matrixH[row][col] = matrixHPc[0][row][col];


		for (int i = 1; i < numberPc; i++)
		{
			sumMatrix(matrixH, matrixHPc[i], 4, 4);
		}
	}

	MainSystem::~MainSystem()
	{
		eta = nullptr;
		ksi = nullptr;
		for (int i = 0; i < numberPc; i++)
		{
			N[i] = nullptr;
			delete[] N[i];
			matrixH[i] = nullptr;
			delete[] matrixH[i];
			matrixC[i] = nullptr;
			delete[] matrixC[i];
			dNdksi[i] = nullptr;
			delete[] dNdksi[i];
			dNdeta[i] = nullptr;
			delete[] dNdeta[i];
			dNdx[i] = nullptr;
			delete[] dNdx[i];
			dNdy[i] = nullptr;
			delete[] dNdy[i];
			for (int j = 0; j < numberPc; j++)
			{
				matrixCPc[i][j] = nullptr;
				matrixHPc[i][j] = nullptr;
				delete[] matrixHPc[i][j];
				delete[] matrixCPc[i][j];
			}
			matrixHPc[i] = nullptr;
			matrixCPc[i] = nullptr;
			delete[] matrixHPc[i];
			delete[] matrixCPc[i];
		}
		N = nullptr;
		delete N;
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
		matrixCPc = nullptr;
		delete matrixCPc;
		delete matrixH;
		matrixC = nullptr;
		delete matrixC;
		for (int i = 0; i < 4; i++)
		{
			vectorP[i] = nullptr;
			delete vectorP[i];
		}
		vectorP = nullptr;
		delete vectorP;
	}
