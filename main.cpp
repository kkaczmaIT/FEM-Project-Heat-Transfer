#include<iostream>
#include<cmath>
#include"Element.h"
#include "Grid.h"
#include"Gauss.h"
#include"MainSystem.h"

using namespace std;

int main()
{
	bool isPrintfMatrix = false;
	MainSystem init(3);
	//Element4 init(4);
	GlobalData globalInfo;
	Grid grid;
	Node localNode;
	Element currentElement;
	double coords[2][4];
	bool isBc[4];
	bool isCalcHBcVectorP = false;
	double** elementMatrixH;
	double** elementMatrixC;
	int* elementsNodesID;
	//string filename = "Test1_4_4.txt";
	string filename = "Test2_4_4_MixGrid.txt";
	//string filename = "Test3_31_31_kwadrat.txt";
	readData(filename, globalInfo, grid);



	//if (filename != "Test3_31_31_kwadrat.txt")
	//{
	//	for (int i = 0; i < grid.getNumberNodes(); i++)
	//		if (grid.getSingleND(i).getY() == 0)
	//				for (int index = 0; index < grid.getNumberNodes(); index++)
	//					while (grid.getSingleND(index).getY() == 0)
	//					{
	//						readData(filename, globalInfo, grid);
	//					}
	//}
	//
	//if(filename == "Test3_31_31_kwadrat.txt")
	//	grid.setSingleND(958, 0.00666666683, -0.0949999988, 1);
	//cout << grid.getNumberNodes() << endl;
	//for (int i = 0; i < grid.getNumberNodes(); i++)
	//	std::cout << "Punkt " << i + 1 << " " << grid.getSingleND(i).getX() << " " << grid.getSingleND(i).getY() << " BC " << grid.getSingleND(i).getBC() << endl;
	//


	//std::cout << "Punkt " << 1 << " " << grid.getSingleND(13).getX() << " " << grid.getSingleND(13).getY() << endl;
	double** matrixHG;
	double** matrixHGHBc;
	double** matrixCG;
	double** matrixHGvectorPGauss;
	double** vectorPLocal;
	double* vectorCdTT0 = new double[grid.getNumberNodes()] {0.0};
	double* vectorPG = new double[grid.getNumberNodes()] {0.0};
	double* solutionVectorFromMatrixHGvectorPGauss = NULL;
	matrixHG = new double* [grid.getNumberNodes()];
	matrixHGHBc = new double* [grid.getNumberNodes()];
	matrixCG = new double* [grid.getNumberNodes()];
	for (int i = 0; i < grid.getNumberNodes(); i++)
	{
		matrixHG[i] = new double[grid.getNumberNodes()];
		matrixHGHBc[i] = new double[grid.getNumberNodes()];
		matrixCG[i] = new double[grid.getNumberNodes()];
	}

	for (int row = 0; row < grid.getNumberNodes(); row++)
		for (int column = 0; column < grid.getNumberNodes(); column++)
		{
			matrixHG[row][column] = 0.0;
			matrixHGHBc[row][column] = 0.0;
			matrixCG[row][column] = 0.0;
		}

	matrixHGvectorPGauss = new double* [grid.getNumberNodes()];
	for (int i = 0; i < grid.getNumberNodes(); i++)
		matrixHGvectorPGauss[i] = new double[grid.getNumberNodes() + 1];

	init.calcdNdksideta();
	init.calcN();
	//cout << "N dla kazdego PC\n";
	//init.displayN();
	//system("pause");

	for (int index = 0; index < grid.getNumberElements(); index++)
	{
		currentElement = grid.getSingleLE(index);
		elementsNodesID = currentElement.getAllID();

		for (int i = 0; i < 4; i++)
		{
			isBc[i] = grid.getSingleBC(currentElement.getSingleID(i) - 1);
		}

		loadElementNodesToCoords(grid, currentElement, coords);

		init.calcJacobianMatrixcalcdNdxdNdy(coords);

		if (isPrintfMatrix)
		{
			std::cout << "Current element\n";
			grid.displaySingleLE(index);

			for (int i = 0; i < 4; i++)
			{
				grid.displaySingleND(elementsNodesID[i] - 1);
				//std::cout << "x: " << coords[0][i] << "\ny: " << coords[1][i] << endl;
			}
			std::cout << "\nNodes\n";

			for (int i = 0; i < 4; i++)
			{
				std::cout << "x: " << coords[0][i] << "\ny: " << coords[1][i] << endl;
			}
		}
		/*
			Calc HBc in depends on BC value of each node
		*/

		for (int id = 0; id < 4; id++)
		{
			if (grid.getSingleND(elementsNodesID[id] - 1).getBC() == 1)
			{
				isCalcHBcVectorP = true;
				currentElement.setHBc(init.calcMatrixHBc(coords, isBc, globalInfo.getAlfa(), globalInfo.getTot()));
				grid.setSingleLE(index, elementsNodesID, currentElement.getHBc());

				if (isPrintfMatrix)
				{
					std::cout << "Element nr 1 macierz HBc: \n";
					displayMatrix(grid.getSingleLE(index).getHBc(), 4, 4);
					std::cout << endl;
				}
				addMatrixHToMatrixHG(matrixHGHBc, currentElement.getHBc(), elementsNodesID, 4, 4);
				break;
			}


		}
		if (isCalcHBcVectorP == true)
		{
			if (isPrintfMatrix)
			{
				std::cout << "vectorPLocal\n";
			}
			vectorPLocal = init.getVectorP();
			for (int indexEl = 0; indexEl < 4; indexEl++)
			{
				if (grid.getSingleBC(currentElement.getSingleID(indexEl) - 1) == 1 && grid.getSingleBC(currentElement.getSingleID((indexEl + 1) % 4) - 1) == 1/* && grid.getSingleBC(currentElement.getSingleID(indexElTmp) - 1) == 1 && (grid.getSingleND(currentElement.getSingleID(indexEl) - 1).getX() == grid.getSingleND(currentElement.getSingleID(indexElTmp) - 1).getX() || grid.getSingleND(currentElement.getSingleID(indexEl) - 1).getY() == grid.getSingleND(currentElement.getSingleID(indexElTmp) - 1).getY())*/)
				{
					vectorPG = addLocalVectorToGlobalVector(vectorPG, vectorPLocal[indexEl], elementsNodesID, 4);
				}

			}

			isCalcHBcVectorP = false;
		}

		for (int nrPc = 0; nrPc < init.getNumberPc(); nrPc++)
		{
			init.calcMatrixHPc(nrPc, globalInfo.getConductivity());
			init.calcMatrixCPc(nrPc, globalInfo.getSpecificHeat(), globalInfo.getDensity());
		}
		init.calcMatrixH();
		init.calcMatrixC();

		elementMatrixH = init.getMatrixH();
		elementMatrixC = init.getMatrixC();
		if (isPrintfMatrix)
		{
			std::cout << "MatrixH nr " << index + 1 << endl;
			displayMatrix(elementMatrixH, 4, 4);
			std::cout << "Element matrix C\n";
			displayMatrix(elementMatrixC, 4, 4);
			std::cout << "Wyznacznik\n";
			init.displayDetJ();
		}
		addMatrixHToMatrixHG(matrixHG, elementMatrixH, elementsNodesID, 4, 4);
		addMatrixCToMatrixCG(matrixCG, elementMatrixC, elementsNodesID, 4, 4);

	}
	if (isPrintfMatrix)
	{
		std::cout << "wezly ilosc: " << grid.getNumberNodes() << endl;
		std::cout << "Matrix HBc\n";
		displayMatrix(matrixHGHBc, grid.getNumberNodes(), grid.getNumberNodes());
		std::cout << "Matrix HG";
		displayMatrix(matrixHG, grid.getNumberNodes(), grid.getNumberNodes());
		std::cout << "Matrix HG po dodaniu wymiany ciepla przez sciany matrix HBc\n";
	}
	sumMatrix(matrixHGHBc, matrixHG, grid.getNumberNodes(), grid.getNumberNodes());
	if (isPrintfMatrix)
	{
		displayMatrix(matrixHGHBc, grid.getNumberNodes(), grid.getNumberNodes());
		std::cout << "Wektor obciazen vectorPG\n";
		displayVector(vectorPG, grid.getNumberNodes());
	}
	for (int row = 0; row < grid.getNumberNodes(); row++)
	{
		for (int column = 0; column < grid.getNumberNodes(); column++)
			matrixHGvectorPGauss[row][column] = matrixHGHBc[row][column];
		matrixHGvectorPGauss[row][grid.getNumberNodes()] = vectorPG[row];
	}

	if (isPrintfMatrix)
	{
		std::cout << "Macierz HG + HBc + vectorP do obliczenia ukladu rownan\n";
		displayMatrix(matrixHGvectorPGauss, grid.getNumberNodes(), grid.getNumberNodes() + 1);
	}
	if (firstLevel(matrixHGvectorPGauss, grid.getNumberNodes()))
	{
		solutionVectorFromMatrixHGvectorPGauss = secondLevel(matrixHGvectorPGauss, grid.getNumberNodes());
		if (isPrintfMatrix)
		{
			std::cout << "Rozwiazana uklad rownan po metodzie eleminacji Gaussa\n ";
			displayMatrix(matrixHGvectorPGauss, grid.getNumberNodes(), grid.getNumberNodes() + 1);
			std::cout << "Tylko rozwiazania ukladu rownan\n";
			displayVector(solutionVectorFromMatrixHGvectorPGauss, grid.getNumberNodes());
			std::cout << endl;
		}
	}

	if (isPrintfMatrix)
	{
		std::cout << " Matrix CG\n";
		displayMatrix(matrixCG, grid.getNumberNodes(), grid.getNumberNodes());

		std::cout << "Matrix H + CG/dT\n";
		std::cout << "Iteration 0\n";
		std::cout << "Matrix H + CG/dT\n";
	}
	sumMatrix(matrixHGHBc, multplicationMatrixByScalar(matrixCG, 1.0 / globalInfo.getSimulationStepTime(), grid.getNumberNodes(), grid.getNumberNodes()), grid.getNumberNodes(), grid.getNumberNodes());
	if (isPrintfMatrix)
	{
		displayMatrix(matrixHGHBc, grid.getNumberNodes(), grid.getNumberNodes());
		std::cout << " Matrix CG\n";
		displayMatrix(matrixCG, grid.getNumberNodes(), grid.getNumberNodes());
		std::cout << "Vector P with C/dT * T0\n";
	}
	double* vectorT0 = new double[grid.getNumberNodes()];
	double** matrixHGHBcCGdTTemp = new double* [grid.getNumberNodes()];
	double* vectorMaxMin;
	double* vectorPGTemp = new double[grid.getNumberNodes()];
	for (int r = 0; r < grid.getNumberNodes(); r++)
	{
		matrixHGHBcCGdTTemp[r] = new double[grid.getNumberNodes() + 1];
	}

	for (int r = 0; r < grid.getNumberNodes(); r++)
		vectorT0[r] = globalInfo.getInitialTemp();
	if (isPrintfMatrix)
	{
		std::cout << "vector0 \n";
		displayVector(vectorT0, grid.getNumberNodes());
	}
	setVectorValuesByVector(vectorPGTemp, vectorPG, grid.getNumberNodes());
	sumVector(vectorPGTemp, multplicationVectorByMatrix(matrixCG, vectorT0, vectorCdTT0, grid.getNumberNodes(), grid.getNumberNodes()), grid.getNumberNodes());
	if (isPrintfMatrix)
	{
		displayVector(vectorPG, grid.getNumberNodes());
	}
	for (int r = 0; r < grid.getNumberNodes(); r++)
	{
		for (int c = 0; c < grid.getNumberNodes(); c++)
		{
			matrixHGHBcCGdTTemp[r][c] = matrixHGHBc[r][c];
		}
		matrixHGHBcCGdTTemp[r][grid.getNumberNodes()] = vectorPGTemp[r];

	}
	if (firstLevel(matrixHGHBcCGdTTemp, grid.getNumberNodes()))
	{
		solutionVectorFromMatrixHGvectorPGauss = secondLevel(matrixHGHBcCGdTTemp, grid.getNumberNodes());
		if (isPrintfMatrix)
		{
			std::cout << "Rozwiazana uklad rownan po metodzie eleminacji Gaussa\n ";
			displayMatrix(matrixHGHBcCGdTTemp, grid.getNumberNodes(), grid.getNumberNodes() + 1);
			std::cout << "Tylko rozwiazania ukladu rownan\n";
			displayVector(solutionVectorFromMatrixHGvectorPGauss, grid.getNumberNodes());
			std::cout << endl;
		}
		vectorMaxMin = findMaxandMinInVector(solutionVectorFromMatrixHGvectorPGauss, grid.getNumberNodes());
		if (isPrintfMatrix)
		{
			std::cout << "Min temp = " << vectorMaxMin[0] << " Max temp = " << vectorMaxMin[1] << endl;
		}
		setVectorValuesByVector(vectorT0, solutionVectorFromMatrixHGvectorPGauss, grid.getNumberNodes());
	}

	if (isPrintfMatrix)
	{
		std::cout << "Iteration 1\n";
		std::cout << "Matrix H + CG/dT\n";
		displayMatrix(matrixHGHBc, grid.getNumberNodes(), grid.getNumberNodes());
		std::cout << " Matrix CG\n";
		displayMatrix(matrixCG, grid.getNumberNodes(), grid.getNumberNodes());

		std::cout << "vector0 \n";
		displayVector(vectorT0, grid.getNumberNodes());
		std::cout << "Vector P with C/dT * T0\n";
	}
	setVectorValuesByVector(vectorPGTemp, vectorPG, grid.getNumberNodes());
	sumVector(vectorPGTemp, multplicationVectorByMatrix(matrixCG, vectorT0, vectorCdTT0, grid.getNumberNodes(), grid.getNumberNodes()), grid.getNumberNodes());
	if (isPrintfMatrix)
	{
		displayVector(vectorPGTemp, grid.getNumberNodes());
	}
	std::cout << "Simulation\n";
	std::cout << "Time[s]\tMinTemp\tMaxTemp\n";
	for (int r = 0; r < grid.getNumberNodes(); r++)
		vectorT0[r] = globalInfo.getInitialTemp();
	for (int time = globalInfo.getSimulationStepTime(); time <= globalInfo.getSimulationTime(); time += globalInfo.getSimulationStepTime())
	{
		setVectorValuesByVector(vectorPGTemp, vectorPG, grid.getNumberNodes());
		sumVector(vectorPGTemp, multplicationVectorByMatrix(matrixCG, vectorT0, vectorCdTT0, grid.getNumberNodes(), grid.getNumberNodes()), grid.getNumberNodes());
		for (int r = 0; r < grid.getNumberNodes(); r++)
		{
			for (int c = 0; c < grid.getNumberNodes(); c++)
			{
				matrixHGHBcCGdTTemp[r][c] = matrixHGHBc[r][c];
			}
			matrixHGHBcCGdTTemp[r][grid.getNumberNodes()] = vectorPGTemp[r];

		}

		if (firstLevel(matrixHGHBcCGdTTemp, grid.getNumberNodes()))
		{
			solutionVectorFromMatrixHGvectorPGauss = secondLevel(matrixHGHBcCGdTTemp, grid.getNumberNodes());
			vectorMaxMin = findMaxandMinInVector(solutionVectorFromMatrixHGvectorPGauss, grid.getNumberNodes());
			setVectorValuesByVector(vectorT0, solutionVectorFromMatrixHGvectorPGauss, grid.getNumberNodes());

			std::cout << time << "\t" << vectorMaxMin[0] << "\t" << vectorMaxMin[1] << endl;
		}
	}


	vectorT0 = nullptr;
	delete vectorT0;
	vectorPGTemp = nullptr;
	delete vectorPGTemp;

	matrixHGHBcCGdTTemp = nullptr;
	delete matrixHGHBcCGdTTemp;
	for (int i = 0; i < grid.getNumberNodes(); i++)
	{
		matrixHG[i] = nullptr;
		delete matrixHG[i];
		matrixHGHBc[i] = nullptr;
		delete matrixHGHBc[i];
		matrixHGvectorPGauss[i] = nullptr;
		delete matrixHGvectorPGauss[i];
		matrixCG[i] = nullptr;
		delete matrixCG[i];
	}
	matrixHG = nullptr;
	matrixHGHBc = nullptr;
	vectorPG = nullptr;
	matrixCG = nullptr;
	delete vectorPG;
	vectorCdTT0 = nullptr;
	delete vectorCdTT0;
	matrixHGvectorPGauss = nullptr;
	delete matrixHGvectorPGauss;
	delete matrixHG;
	delete matrixHGHBc;
	delete matrixCG;
	elementMatrixH = nullptr;
	delete elementMatrixH;
	elementsNodesID = nullptr;
	delete elementsNodesID;
	return 0;
}


