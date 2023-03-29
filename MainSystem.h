#include<iostream>
#include"utility.h"

#pragma once
class MainSystem
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
	double** N;
	// Wektor obciazen
	double** vectorP;
	double*** matrixHPc;
	double*** matrixCPc;
	double** matrixC;
	double** matrixH;
	double matrixHBc[4][4];
public:
	friend void displayMatrix(double** matrix, int row, int col);
	friend double** transpositionMatrix(double** matrix);
	friend double** multplicationMatrixByScalar(double** matrix, double scalar, int row, int col);
	friend double** sumMatrix(double** matrix1, double** matrix2, int row, int col);
	friend double* sumVector(double* vector1, double* vector2, int col);
	friend double** transpositionMatrix(double** matrix);
	friend double** multiplicationMatrix(double* verMatrix, double* horMatrix);
	friend double* multplicationVectorByScalar(double* vector, double scalar, int col);
	friend void displayVector(double* vector, int col);
	friend void setVectorValuesByVector(double* vectorReceiver, double* vectorSender, int col);
	MainSystem(int numberNodes);

	double** getVectorP();

	double** getMatrixH();

	void calcN();

	void displayN();

	void calcdNdksideta();

	double** calcMatrixHBc(double coords[2][4], bool isBc[4], double alpha, double tot);
	

	void calcMatrixCPc(int nrPc, double c, double p);

	void setMatrixCPcValue(int nrPc, double** matrixCTemp);




	void calcMatrixC();

	double** getMatrixC();


	/*Test function weights*/
	void displayMultplicationWeight();
	/*Test function weights*/

	void displayElementdksideta();

	void calcNdxNdy(int nrPc, double matrix[2][2]);

	void calcJacobianMatrixcalcdNdxdNdy(double coords[2][4]);

	void displayElementdxdy();
	
	int getNumberPc();

	void setMatrixHPcValue(int nrPc, double** matrix);

	void displayPCElements();

	void calcMatrixHPc(int nrPc, int k);

	void displayMatrixHPc(int nrPc);

	void displayDetJ();

	void calcMatrixH();

	~MainSystem();

};
