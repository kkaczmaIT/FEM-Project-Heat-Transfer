#include "Node.h"

Node::Node() {}

Node::Node(double x, double y, int BC)
{
	this->x = x;
	this->y = y;
	this->BC = BC;
}

void Node::setX(double x)
{
	this->x = x;
}

void Node::setY(double y)
{
	this->y = y;
}

void Node::setBC(int BC)
{
	this->BC = BC;
}

double Node::getX()
{
	return this->x;
}

double Node::getY()
{
	return this->y;
}

int Node::getBC()
{
	return this->BC;
}

Node::~Node() {};