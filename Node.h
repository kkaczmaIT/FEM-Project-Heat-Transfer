#pragma once

class Node
{
private:
	double x;
	double y;
	int BC;
public:
	Node();
	Node(double x, double y, int BC);

	//setters
	void setX(double x);

	void setY(double y);

	void setBC(int BC);

	// getters
	double getX();

	double getY();

	int getBC();

	~Node();
};