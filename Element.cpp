#include "Element.h"

Element::Element() {}

Element::Element(int id[4])
{
	for (int row = 0; row < 4; row++)
	{
		for (int column = 0; column < 4; column++)
			HBc[row][column] = 0.0;
	}
}

// setters
void Element::setHBc(double** matrix)
{
	for (int row = 0; row < 4; row++)
		for (int column = 0; column < 4; column++)
			HBc[row][column] = matrix[row][column];
}

void Element::setAllID(int id[4])
{
	for (int i = 0; i < 4; i++)
	{
		this->id[i] = id[i];
	}
}

void Element::setSingleID(int id, int index)
{
	if (index < 4 && index >= 0)
		this->id[index] = id;
}


// getters

void Element::getHBc(double matrix[4][4])
{
	for (int row = 0; row < 4; row++)
		for (int column = 0; column < 4; column++)
			matrix[row][column] = HBc[row][column];
}

double** Element::getHBc()
{
	double** matrix;
	matrix = new double* [4];
	for (int i = 0; i < 4; i++)
		matrix[i] = new double[4];
	for (int row = 0; row < 4; row++)
		for (int column = 0; column < 4; column++)
			matrix[row][column] = HBc[row][column];

	return matrix;
	for (int i = 0; i < 4; i++)
	{
		matrix[i] = nullptr;
		delete matrix[i];
	}
	matrix = nullptr;
	delete matrix;
}

int* Element::getAllID()
{
	return this->id;
}

int Element::getSingleID(int index)
{
	if (index >= 0 && index < 4)
		return this->id[index];
	else
		return -1;
}

Element::~Element() {}