#pragma once
#include<iostream>
#include"Node.h"
#include"Element.h"
class Grid
{
private:
	int numberNodes;
	int numberElements;
	Node* ND;
	Element* LE;
public:
	Grid();
	Grid(int numberNodes, int numberElements);

	// Setters
	void setNumberNodes(int numberNodes);

	void setNumberElements(int numberElements);

	void setSingleND(int index, double x, double y, int BC = 0);

	void setSingleBC(int index, int BC);

	void setND(int** nodes);

	void setSingleLE(int index, int id[4]);

	void setSingleLE(int index, int id[4], double** HBc);

	void setLE(int** elements, int* BC);

	// Getters
	int getNumberNodes();

	int getNumberElements();

	Node getSingleND(int index);

	Node* getND();

	Element getSingleLE(int index);

	int getSingleBC(int index);

	Element* getLE();

	//display property
	void displaySingleND(int index);

	void displaySingleLE(int index);

	~Grid();
};