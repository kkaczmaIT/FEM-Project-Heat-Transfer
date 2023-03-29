#include"Grid.h"

Grid::Grid() {}
Grid::Grid(int numberNodes, int numberElements)
{
	this->numberNodes = numberNodes;
	this->numberElements = numberElements;
	this->ND = new Node[this->numberNodes];
	this->LE = new Element[this->numberElements];
}

// Setters
void Grid::setNumberNodes(int numberNodes)
{
	this->numberNodes = numberNodes;
	this->ND = new Node[this->numberNodes];
}

void Grid::setNumberElements(int numberElements)
{
	this->numberElements = numberElements;
	this->LE = new Element[this->numberElements];
}

void Grid::setSingleND(int index, double x, double y, int BC)
{
	if (index >= 0 && index < this->numberNodes)
	{
		ND[index].setX(x);
		ND[index].setY(y);
		ND[index].setBC(BC);
	}
}

void Grid::setSingleBC(int index, int BC)
{
	ND[index].setBC(BC);
}

void Grid::setND(int** nodes)
{
	for (int index = 0; index < this->numberNodes; index++)
	{
		this->setSingleND(index, nodes[index][0], nodes[index][1]);
	}
}

void Grid::setSingleLE(int index, int id[4])
{
	LE[index].setAllID(id);

}

void Grid::setSingleLE(int index, int id[4], double** HBc)
{
	LE[index].setAllID(id);
	LE[index].setHBc(HBc);

}

void Grid::setLE(int** elements, int* BC)
{
	for (int index = 0; index < this->numberElements; index++)
	{
		this->setSingleLE(index, elements[index]);
	}
}

// Getters
int Grid::getNumberNodes()
{
	return this->numberNodes;
}

int Grid::getNumberElements()
{
	return this->numberElements;
}

Node Grid::getSingleND(int index)
{
	if (index >= 0 && index < this->numberNodes)
	{
		return ND[index];
	}
}

Node* Grid::getND()
{
	return this->ND;

}

Element Grid::getSingleLE(int index)
{
	return LE[index];
}

int Grid::getSingleBC(int index)
{
	return ND[index].getBC();
}

Element* Grid::getLE()
{
	return this->LE;

}

//display property
void Grid::displaySingleND(int index)
{
	std::cout << "x = " << ND[index].getX() << " y = " << ND[index].getY();
	std::cout << " BC = " << ND[index].getBC() << std::endl;
}

void Grid::displaySingleLE(int index)
{
	std::cout << "Nodes: ";
	for (int i = 0; i < 4; i++)
	{
		std::cout << LE[index].getSingleID(i);
		i <= 2 ? std::cout << ", " : std::cout << "\n";
	}


}

Grid::~Grid()
{
	ND = nullptr;
	LE = nullptr;
	delete ND;
	delete LE;
}