#pragma once

class Element
{
private:
	int id[4];
	double HBc[4][4];
public:
	Element();
	Element(int id[4]);

	// setters
	void setHBc(double** matrix);

	void setAllID(int id[4]);

	void setSingleID(int id, int index);

	// getters

	void getHBc(double matrix[4][4]);

	double** getHBc();

	int* getAllID();

	int getSingleID(int index);

	~Element();
};
