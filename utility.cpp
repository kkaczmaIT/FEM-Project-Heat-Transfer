#include "utility.h"
#include<string>

using namespace std;

propertyCode translate(string property)
{
	if (property == "SimulationTime") return simulationTime;
	if (property == "SimulationStepTime") return simulationStepTime;
	if (property == "Conductivity") return conductivity;
	if (property == "Alfa") return alfa;
	if (property == "Tot") return tot;
	if (property == "InitialTemp") return initialTemp;
	if (property == "Density") return density;
	if (property == "SpecificHeat") return specificHeat;
	if (property == "Nodes") return nodesNumber;
	if (property == "Elements") return elementsNumber;
	if (property == "*Node") return node;
	if (property == "*Element") return element;
	if (property == "*BC") return BCkey;
}

bool readData(std::string filename, GlobalData& globalInfo, Grid &grid)
{
	fstream File(filename, ios::in);
	string property, input, valueInStr;
	Node* nodes = NULL;
	double value, xNode, yNode;
	double thirdFromEnd[2];
	int counter;
	int numberNode, numberElement;
	int elementNode[4] = { 0, 0, 0, 0 };
	int code = 0;
	//xNode = 0.1;
	//yNode = 0.1;
	numberNode = 0;
	numberElement = 0;
	valueInStr = "";
	//std::cout << xNode << yNode << numberNode << endl;
	if (File.good())
	{
		while (!File.eof())
		{
			File >> property;
			if (property == "Nodes")
			{
				File >> property;
				property = "Nodes";
			}
			else if (property == "Elements")
			{
				File >> property;
				property = "Elements";
			}
			else if (property == "*Element,")
			{
				File >> property;
				property = "*Element";
			}
			if (property != "*Node" && property != "*Element" && property != "*BC")
				File >> value;
			//std::cout << "Property = " << property << " value = " << value << endl;
			switch (translate(property))
			{
			case simulationTime:
				globalInfo.setSimulationTime(value);
				//std::cout << globalInfo.getSimulationTime() << endl;
				break;
			case simulationStepTime:
				globalInfo.setSimulationStepTime(value);
				//std::cout << globalInfo.getSimulationStepTime() << endl;
				break;
			case conductivity:
				globalInfo.setConductivity(value);
				//std::cout << globalInfo.getConductivity() << endl;
				break;
			case tot:
				globalInfo.setTot(value);
				//std::cout << globalInfo.getTot() << endl;
				break;
			case initialTemp:
				globalInfo.setInitialTemp(value);
				//std::cout << globalInfo.getInitialTemp() << endl;
				break;
			case density:
				globalInfo.setDensity(value);
				//std::cout << globalInfo.getDensity() << endl;
				break;
			case specificHeat:
				globalInfo.setSpecificHeat(value);
				//std::cout << globalInfo.getSpecificHeat() << endl;
				break;
			case alfa:
				globalInfo.setAlfa(value);
				//std::cout << globalInfo.getAlfa() << endl;
				break;
			case nodesNumber:
				grid.setNumberNodes(value);
				nodes = new Node[value];
				//std::cout << grid.getNumberNodes();
				break;
			case elementsNumber:
				grid.setNumberElements(value);
				//std::cout << grid.getNumberElements();
				break;
			case node:
				//std::cout << "node\n";
				valueInStr = "";
				//std::cout << grid.getNumberNodes() << " " << xNode << " " << yNode << " " << numberNode << endl;
				for (int index = 0; index < grid.getNumberNodes() + 1; index++)
				{
					code++;
					valueInStr = "";
					getline(File, input);
					counter = 0;
					//std::cout << input << endl;
					for (int pos = 0; pos < input.size(); pos++)
					{
						//std::cout << "inpuut: " << input[pos] << endl;
						if (input[pos] != ',')
						{
							//std::cout << "start: ";
							valueInStr += input[pos];
							//std::cout << valueInStr << endl;
						}
						else if (input[pos] == ',')
						{
							//std::cout << "przecinek\n";
							if (counter == 0)
							{
								//std::cout << "numer wezla\n";
								//std::cout << "liczba: " << valueInStr << endl;
								numberNode = stoi(valueInStr);
							}
							else if (counter == 1)
							{
								//std::cout << "wspolrzedna x\n";
								if (valueInStr[0] == '-')
								{
									valueInStr[0] = ' ';
									xNode = stod(valueInStr);
									xNode = -xNode;
								}
								else
									xNode = stod(valueInStr);
							}
							counter++;
							valueInStr = "";
						}
						if (pos == input.size() - 1 && counter == 2)
						{
							if (valueInStr[0] == '-')
							{
								valueInStr[0] = ' ';
								yNode = stod(valueInStr);
								yNode = -yNode;
							}
							else
								yNode = stod(valueInStr);


							//grid.displaySingleND(numberNode - 1);
							//grid.getSingleND(numberNode - 1).setY(yNode);

						}
					}
					if (counter == 2)
					{
						grid.setSingleND(numberNode - 1, xNode, yNode);
						//std::cout << "z obiektu grid Node " << numberNode << "x = " << grid.getSingleND(numberNode - 1).getX() << " y = " << grid.getSingleND(numberNode - 1).getY() << endl;
						if (numberNode - 1 + 3 == grid.getNumberNodes())
						{
							thirdFromEnd[0] = grid.getSingleND(numberNode - 1).getX();
							thirdFromEnd[1] = grid.getSingleND(numberNode - 1).getY();
						}
					}

					//std::cout << "Skladowe " << numberNode << " " << xNode << " " << yNode << endl;

				}
				//cout << "Po wczytaniu danych\n";
				for (int i = 0; i < grid.getNumberNodes(); i++)
				{
					//std::cout << "z obiektu grid Node " << i + 1 << "x = " << grid.getSingleND(i).getX() << " y = " << grid.getSingleND(i).getY() << endl;
					nodes[i].setX(grid.getSingleND(i).getX());
					nodes[i].setY(grid.getSingleND(i).getY());
				}
				//cout << "code: " << code << endl;*/
				//nodes = grid.getND();

				break;

			case element:
				//std::cout << "element class\n";
				for (int index = 0; index < grid.getNumberElements() + 1; index++)
				{
					getline(File, input);
					//std::cout << input << endl;
					valueInStr = "";
					for (int pos = 0, counter = 0; pos < input.size() + 1; pos++)
					{
						if (input[pos] != ',')
							valueInStr += input[pos];
						else if (input[pos] == ',')
						{
							if (counter == 0)
							{
								numberElement = stoi(valueInStr);
							}
							if (counter != 0)
							{
								elementNode[counter - 1] = stoi(valueInStr);
							}

							//std::cout << "elementNode: " << elementNode[counter - 1] << endl;
							valueInStr = "";
							counter++;
						}
						if (pos == input.size() - 1)
						{
							elementNode[3] = stoi(valueInStr);
							//std::cout << "elementNode: " << elementNode[3] << endl;
						}
					}

					grid.setSingleLE(numberElement - 1, elementNode);
				}
				break;

			case BCkey:
				valueInStr = "";
				getline(File, input);
				getline(File, input);
				//std::cout << input << endl;
				for (int index = 0; index < grid.getNumberNodes() + 1; index++)
				{
					valueInStr = "";
					for (int pos = 0; pos < input.size() + 1; pos++)
					{
						if (input[pos] != ',')
						{
							valueInStr += input[pos];
						}
						else if (input[pos] == ',')
						{
							numberNode = stoi(valueInStr);
							//std::cout << "numberElement: " << numberNode << endl;
							valueInStr = "";
							grid.setSingleBC(numberNode - 1, 1);
						}

						if (pos == input.size() - 1)
						{
							numberNode = stoi(valueInStr);
							grid.setSingleBC(numberNode - 1, 1);
							valueInStr = "";
						}
						//std::cout << "BC: " << valueInStr << endl;
					}
				}
				break;
			default:
				break;
			}

		}

		File.close();

		//cout << "node elements: " << grid.getNumberNodes() << endl;
		//cout << "Node element X: " << grid.getSingleND(13).getX() << " Y: " << grid.getSingleND(13).getY() << endl;
		for (int i = 0; i < grid.getNumberNodes(); i++)
		{
			if (i != grid.getNumberNodes() - 3)
			{
				grid.setSingleND(i, nodes[i].getX(), nodes[i].getY(), grid.getSingleND(i).getBC());
			}
			else
			{
				grid.setSingleND(i, thirdFromEnd[0], thirdFromEnd[1], grid.getSingleND(i).getBC());
			}

			//std::cout << "poza casem grid Node " << i + 1 << "x = " << grid.getSingleND(i).getX() << " y = " << grid.getSingleND(i).getY() << endl;
		}

		//cout << "code poza casem: " << code << endl;
	}
	else
	{
		std::cout << "Error\n";
	}

	return true;
}



double* findMaxandMinInVector(double* vector, int row)
{
	double* resVector = new double[2];
	double min = vector[0];
	double max = vector[0];
	for (int r = 1; r < row; r++)
	{
		if (min > vector[r])
			min = vector[r];
		if (max < vector[r])
			max = vector[r];
	}
	resVector[0] = min;
	resVector[1] = max;
	return resVector;
}

double** multiplicationMatrix(double* verMatrix, double* horMatrix)
{
	double** resMatrix = new double* [4];
	for (int i = 0; i < 4; i++)
	{
		resMatrix[i] = new double[4];
	}
	//std::cout << endl;
	for (int row = 0; row < 4; row++)
		for (int col = 0; col < 4; col++)
			resMatrix[row][col] = verMatrix[row] * horMatrix[col];
	return resMatrix;
	resMatrix = nullptr;
	delete resMatrix;
}

double** divideMatrixByScalar(double** matrix, double scalar, int row, int column)
{
	for (int r = 0; r < row; r++)
	{
		for (int c = 0; c < column; c++)
		{
			matrix[r][c] /= scalar;
		}
	}
	return matrix;
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
	std::cout << endl;
	for (int r = 0; r < row; r++)
	{
		for (int c = 0; c < col; c++)
			std::cout << matrix[r][c] << "\t";
		std::cout << endl;
	}
	std::cout << endl;
}

void displayVector(double* vector, int col)
{
	std::cout << endl;
	for (int c = 0; c < col; c++)
		std::cout << vector[c] << "\t";
	std::cout << endl;
}

double** sumMatrix(double** matrix1, double** matrix2, int row, int col)
{
	for (int r = 0; r < row; r++)
		for (int c = 0; c < col; c++)
			matrix1[r][c] += matrix2[r][c];
	return matrix1;
}

double* sumVector(double* vector1, double* vector2, int col)
{
	for (int c = 0; c < col; c++)
		vector1[c] += vector2[c];
	return vector1;
}

double* multplicationVectorByMatrix(double** matrix, double* vector, double* resVector, int row, int column)
{
	for (int i = 0; i < row; i++)
		resVector[i] = 0.0;
	for (int r = 0; r < row; r++)
	{
		for (int c = 0; c < column; c++)
		{
			resVector[r] += matrix[r][c] * vector[c];
		}
	}
	return resVector;
}

double** multplicationMatrixByScalar(double** matrix, double scalar, int row, int col)
{
	for (int r = 0; r < row; r++)
		for (int c = 0; c < col; c++)
			matrix[r][c] *= scalar;
	return matrix;
}

double* multplicationVectorByScalar(double* vector, double scalar, int col)
{
	for (int c = 0; c < col; c++)
		vector[c] *= scalar;
	return vector;
}


void setVectorValuesByVector(double* vectorReceiver, double* vectorSender, int col)
{
	for (int c = 0; c < col; c++)
	{
		vectorReceiver[c] = vectorSender[c];
	}
}

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

void loadElementNodesToCoords(Grid grid, Element el, double coords[2][4])
{
	Node localNode;
	int indexNode;
	for (int i = 0; i < 4; i++)
	{
		indexNode = el.getSingleID(i) - 1;
		localNode = grid.getSingleND(indexNode);
		coords[0][i] = localNode.getX();
		coords[1][i] = localNode.getY();
	}
}

// Functions to handle calculate matrix HG start


double** addMatrixHToMatrixHG(double** matrixHG, double** matrixH, int* nodesID, int row, int col)
{
	for (int r = 0; r < row; r++)
		for (int c = 0; c < col; c++)
			matrixHG[nodesID[r] - 1][nodesID[c] - 1] += matrixH[r][c];
	return matrixHG;
}

double** addMatrixCToMatrixCG(double** matrixCG, double** matrixC, int* nodesID, int row, int col)
{
	for (int r = 0; r < row; r++)
		for (int c = 0; c < col; c++)
			matrixCG[nodesID[r] - 1][nodesID[c] - 1] += matrixC[r][c];
	return matrixCG;
}

// Functions to handle calculate matrix HG end

double* addLocalVectorToGlobalVector(double* vectorGlobal, double* vectorLocal, int* nodesID, int col)
{
	for (int c = 0; c < col; c++)
		vectorGlobal[nodesID[c] - 1] += vectorLocal[c];
	return vectorGlobal;
}