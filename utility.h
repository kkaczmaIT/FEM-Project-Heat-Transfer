#include<iostream>
#include<fstream>
#include<string>
#include "Globaldata.h"
#include "Grid.h"

//#pragma once
enum propertyCode
{
	simulationTime,
	simulationStepTime,
	conductivity,
	alfa,
	tot,
	initialTemp,
	density,
	specificHeat,
	nodesNumber,
	elementsNumber,
	node,
	element,
	BCkey
};

propertyCode translate(std::string property);

bool readData(std::string filename, GlobalData& globalInfo, Grid &grid);


//Prototype
void displayMatrix(double** matrix, int row, int col);
double** transpositionMatrix(double** matrix);
double** multplicationMatrixByScalar(double** matrix, double scalar, int row, int col);
double** sumMatrix(double** matrix1, double** matrix2, int row, int col);
double* sumVector(double* vector1, double* vector2, int col);
double** transpositionMatrix(double** matrix);
double** multiplicationMatrix(double* verMatrix, double* horMatrix);
double* multplicationVectorByScalar(double* vector, double scalar, int col);
void displayVector(double* vector, int col);
void setVectorValuesByVector(double* vectorReceiver, double* vectorSender, int col);
double** divideMatrixByScalar(double** matrix, double scalar, int row, int column);
double* multplicationVectorByMatrix(double** matrix, double* vector, double* resVector, int row, int column);
double* findMaxandMinInVector(double* vector, int row);

double N1(double ksi, double eta);

double N2(double ksi, double eta);

double N3(double ksi, double eta);

double N4(double ksi, double eta);

double N1ksi(double eta);

double N2ksi(double eta);

double N3ksi(double eta);

double N4ksi(double eta);

double N1eta(double ksi);

double N2eta(double ksi);

double N3eta(double ksi);

double N4eta(double ksi);

void loadElementNodesToCoords(Grid grid, Element el, double coords[2][4]);

double* addLocalVectorToGlobalVector(double* vectorGlobal, double* vectorLocal, int* nodesID, int col);

double** addMatrixCToMatrixCG(double** matrixCG, double** matrixC, int* nodesID, int row, int col);

double** addMatrixHToMatrixHG(double** matrixHG, double** matrixH, int* nodesID, int row, int col);

