#include<iostream> 
#include<fstream>
#include<string>

using namespace std;

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

class Node 
{
	private:
		double x;
		double y;
		int BC;
	public:
		Node() {}
		Node(double x, double y, int BC)
		{
			this->x = x;
			this->y = y;
			this->BC = BC;
		}

		//setters
		void setX(double x)
		{
			this->x = x;
		}

		void setY(double y)
		{
			this->y = y;
		}

		void setBC(int BC)
		{
			this->BC = BC;
		}

		// getters
		double getX()
		{
			return this->x;
		}

		double getY()
		{
			return this->y;
		}

		int getBC()
		{
			return this->BC;
		}

		~Node() {}
};

class Element
{
	private:
		int id[4];
	public:
		Element() {}
		Element(int id[4])
		{
			for (int i = 0; i < 4; i++)
			{
				this->id[i] = id[i];
			}
			
		}

		// setters
		void setAllID(int id[4])
		{
			for (int i = 0; i < 4; i++)
			{
				this->id[i] = id[i];
			}
		}

		void setSingleID(int id, int index)
		{
			if(index < 4 && index >= 0)
				this->id[index] = id;
		}


		// getters
		
		int* getAllID()
		{
			return this->id;
		}
		
		int getSingleID(int index)
		{
			if (index >= 0 && index < 4)
				return this->id[index];
			else
				return -1;
		}


		~Element() {}
};

class Grid
{
	private:
		int numberNodes;
		int numberElements;
		Node* ND;
		Element* LE;
	public:
		Grid() {}
		Grid(int numberNodes, int numberElements)
		{
			this->numberNodes = numberNodes;
			this->numberElements = numberElements;
			this->ND = new Node[this->numberNodes];
			this->LE = new Element[this->numberElements];
		}

		// Setters
		void setNumberNodes(int numberNodes) 
		{
			this->numberNodes = numberNodes;
			this->ND = new Node[this->numberNodes];
		}

		void setNumberElements(int numberElements)
		{
			this->numberElements = numberElements;
			this->LE = new Element[this->numberElements];
		}

		void setSingleND(int index, double x, double y, int BC = 0)
		{
			if (index >= 0 && index < this->numberNodes)
			{
				ND[index].setX(x);
				ND[index].setY(y);
				ND[index].setBC(BC);
			}
		}
		
		void setSingleBC(int index, int BC)
		{
			ND[index].setBC(BC);
		}

		void setND(int** nodes)
		{
			for (int index = 0; index < this->numberNodes; index++)
			{
				this->setSingleND(index, nodes[index][0], nodes[index][1]);
			}
		}

		void setSingleLE(int index, int id[4])
		{
			LE[index].setAllID(id);
			
		}


		void setLE(int** elements, int* BC)
		{
			for (int index = 0; index < this->numberElements; index++)
			{
				this->setSingleLE(index, elements[index]);
			}
		}

		// Getters
		int getNumberNodes()
		{
			return this->numberNodes;
		}

		int getNumberElements()
		{
			return this->numberElements;
		}

		Node getSingleND(int index)
		{
			if (index >= 0 && index < this->numberNodes)
			{
				return ND[index];
			}
		}

		Node* getND()
		{
			return this->ND;
			
		}

		Element getSingleLE(int index)
		{
			return LE[index];
		}

		int getSingleBC(int index)
		{
			return ND[index].getBC();
		}

		Element* getLE()
		{
			return this->LE;
			
		}
		
		//display property
		void displaySingleND(int index)
		{
			cout << "x = " << ND[index].getX() << " y = " << ND[index].getY();
			cout << " BC = " << ND[index].getBC() << endl;
		}

		void displaySingleLE(int index)
		{
			cout << "Nodes: ";
			for (int i = 0; i < 4; i++)
			{
				cout << LE[index].getSingleID(i);
				i <= 2 ? cout << ", " : cout << "\n";
			}
			
			
		}

		~Grid()
		{
			ND = nullptr;
			LE = nullptr;
			delete ND;
			delete LE;
		}
};

class GlobalData
{
	private:
		int simulationTime;
		int simulationStepTime;
		double conductivity;
		double alfa;
		double tot;
		double initialTemp;
		double density;
		double specificHeat;
	public:
		GlobalData() {}
		GlobalData(int simulationTime, int simulationStepTime, double conductivity, double alfa, double tot, double initialTemp, double density, double specificHeat)
		{
			this->simulationTime = simulationTime;
			this->simulationStepTime = simulationStepTime;
			this->conductivity = conductivity;
			this->alfa = alfa;
			this->tot = tot;
			this->initialTemp = initialTemp;
			this->density = density;
			this->specificHeat = specificHeat;
		}

		// setters
		
		void setSimulationTime(int simulationTime)
		{
			this->simulationTime = simulationTime;
		}

		void setSimulationStepTime(int simulationStepTime)
		{
			this->simulationStepTime = simulationStepTime;
		}

		void setConductivity(double conductivity)
		{
			this->conductivity = conductivity;
		}

		void setAlfa(double alfa)
		{
			this->alfa = alfa;
		}

		void setTot(double tot)
		{
			this->tot = tot;
		}

		void setInitialTemp(double initialTemp)
		{
			this->initialTemp = initialTemp;
		}

		void setDensity(double density)
		{
			this->density = density;
		}

		void setSpecificHeat(double specificHeat)
		{
			this->specificHeat = specificHeat;
		}

		// getters
		int getSimulationTime()
		{
			return this->simulationTime;
		}

		int getSimulationStepTime()
		{
			return this->simulationStepTime;
		}

		double getConductivity()
		{
			return this->conductivity;
		}

		double getAlfa()
		{
			return this->alfa;
		}

		double getTot()
		{
			return this->tot;
		}

		double getInitialTemp()
		{
			return this->initialTemp;
		}

		double getDensity()
		{
			return this->density;
		}

		double getSpecificHeat()
		{
			return this->specificHeat;
		}
		~GlobalData() {}
};

bool readData(string filename, GlobalData &globalInfo, Grid &grid)
{
	fstream File(filename , ios::in);
	string property, input, valueInStr;
	double value, xNode, yNode;
	int numberNode, numberElement;
	int elementNode[4] = {0, 0, 0, 0};
	//int counter;
	xNode = 0.1;
	yNode = 0.1;
	numberNode = 0;
	numberElement = 0;
	valueInStr = "";
	//cout << xNode << yNode << numberNode << endl;
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
			else if(property == "Elements")
			{
				File >> property;
				property = "Elements";
			}
			else if (property == "*Element,")
			{
				File >> property;
				property = "*Element";
			}
			if(property != "*Node" && property != "*Element" && property != "*BC")
				File >> value;
			//cout << "Property = " << property << " value = " << value << endl;
			switch (translate(property))
			{
				case simulationTime:
					globalInfo.setSimulationTime(value);
					//cout << globalInfo.getSimulationTime() << endl;
					break;
				case simulationStepTime:
					globalInfo.setSimulationStepTime(value);
					//cout << globalInfo.getSimulationStepTime() << endl;
					break;
				case conductivity:
					globalInfo.setConductivity(value);
					//cout << globalInfo.getConductivity() << endl;
					break;
				case tot:
					globalInfo.setTot(value);
					//cout << globalInfo.getTot() << endl;
					break;
				case initialTemp:
					globalInfo.setInitialTemp(value);
					//cout << globalInfo.getInitialTemp() << endl;
					break;
				case density:
					globalInfo.setDensity(value);
					//cout << globalInfo.getDensity() << endl;
					break;
				case specificHeat:
					globalInfo.setSpecificHeat(value);
					//cout << globalInfo.getSpecificHeat() << endl;
					break;
				case alfa:
					globalInfo.setAlfa(value);
					//cout << globalInfo.getAlfa() << endl;
					break;
				case nodesNumber:
					grid.setNumberNodes(value);
					//cout << grid.getNumberNodes();
					break;
				case elementsNumber:
					grid.setNumberElements(value);
					//cout << grid.getNumberElements();
					break;
				case node:
					//cout << "node\n";
					valueInStr = "";
					//cout << grid.getNumberNodes() << " " << xNode << " " << yNode << " " << numberNode << endl;
					for (int index = 0; index < grid.getNumberNodes() + 1; index++)
					{
						getline(File, input);
						//cout << input << endl;
						for (int counter = 0, pos = 0; pos < input.size(); pos++)
						{
							//cout << "inpuut: " << input[pos] << endl;
							if (input[pos] != ',')
							{
								//cout << "start: ";
								valueInStr += input[pos];
								//cout << valueInStr << endl;
							}
							else if (input[pos] == ',')
							{
								//cout << "przecinek\n";
								if (counter == 0)
								{
									//cout << "numer wezla\n";
									//cout << "liczba: " << valueInStr << endl;
									numberNode = stoi(valueInStr);
								}
								else if (counter == 1)
								{
									//cout << "wspolrzedna x\n";
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
							if (pos == input.size() - 1)
							{
								if (valueInStr[0] == '-')
								{
									valueInStr[0] = ' ';
									yNode = stod(valueInStr);
								}
								else
									yNode = stod(valueInStr);

								grid.setSingleND(numberNode - 1, xNode, yNode);
							}
						}
						
						valueInStr = "";
						//cout << "Skladowe " << numberNode << " " << xNode << " " << yNode << endl;
						

					}
					
					break;

				case element:
					//cout << "element class\n";
					for (int index = 0; index < grid.getNumberElements() + 1; index++)
					{
						getline(File, input);
						//cout << input << endl;
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
								
								//cout << "elementNode: " << elementNode[counter - 1] << endl;
								valueInStr = "";
								counter++;
							}
							if (pos == input.size() - 1)
							{
								elementNode[3] = stoi(valueInStr);
								//cout << "elementNode: " << elementNode[3] << endl;
							}
						}

						grid.setSingleLE(numberElement - 1, elementNode);
					}
					break;

				case BCkey:
					valueInStr = "";
					getline(File, input);
					getline(File, input);
					//cout << input << endl;
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
								//cout << "numberElement: " << numberNode << endl;
								valueInStr = "";
								grid.setSingleBC(numberNode - 1, 1);
							}

							if (pos == input.size() - 1)
							{
								numberNode = stoi(valueInStr);
								grid.setSingleBC(numberNode - 1, 1);
								valueInStr = "";
							}
							//cout << "BC: " << valueInStr << endl;
						}
					}
					break;
				default:
					break;
			}
			
		}
		
		File.close();
	}
	else
	{
		cout << "Error\n";
	}

	return true;
}

//int main()
//{
//	GlobalData globalInfo;
//	Grid grid;
//	Node tmp;
//	//readData("Test1_4_4.txt", globalInfo, grid);
//	//readData("Test2_4_4_MixGrid.txt", globalInfo, grid);
//	readData("Test3_31_31_kwadrat.txt", globalInfo, grid);
//	cout << "Grid's Data\n";
//	cout << "Simulation Time: " <<globalInfo.getSimulationTime() << endl;
//	cout << "Simulation Step Time: " << globalInfo.getSimulationStepTime() << endl;
//	cout << "Conductivity: " << globalInfo.getConductivity() << endl;
//	cout << "Alfa: " << globalInfo.getAlfa() << endl;
//	cout << "Tot: " << globalInfo.getTot() << endl;
//	cout << "Inititial Temp: " << globalInfo.getInitialTemp() << endl;
//	cout << "Density: " << globalInfo.getDensity() << endl;
//	cout << "Specific Heat: " << globalInfo.getSpecificHeat() << endl;
//
//	cout << "Nodes\n";
//	cout.precision(10);
//	for (int index = 0; index < grid.getNumberNodes(); index++)
//	{	
//		cout << index + 1 << " ";
//		grid.displaySingleND(index);
//	}
//	cout << endl << "ELemnents\n";
//	for (int index = 0; index < grid.getNumberElements(); index++)
//	{
//		cout << index + 1 << " ";
//		grid.displaySingleLE(index);
//	}
//
//	return 0;
//}