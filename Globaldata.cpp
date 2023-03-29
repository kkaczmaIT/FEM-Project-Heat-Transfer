#include"Globaldata.h"

GlobalData::GlobalData() {}
GlobalData::GlobalData(int simulationTime, int simulationStepTime, double conductivity, double alfa, double tot, double initialTemp, double density, double specificHeat)
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

void GlobalData::setSimulationTime(int simulationTime)
{
	this->simulationTime = simulationTime;
}

void GlobalData::setSimulationStepTime(int simulationStepTime)
{
	this->simulationStepTime = simulationStepTime;
}

void GlobalData::setConductivity(double conductivity)
{
	this->conductivity = conductivity;
}

void GlobalData::setAlfa(double alfa)
{
	this->alfa = alfa;
}

void GlobalData::setTot(double tot)
{
	this->tot = tot;
}

void GlobalData::setInitialTemp(double initialTemp)
{
	this->initialTemp = initialTemp;
}

void GlobalData::setDensity(double density)
{
	this->density = density;
}

void GlobalData::setSpecificHeat(double specificHeat)
{
	this->specificHeat = specificHeat;
}

// getters
int GlobalData::getSimulationTime()
{
	return this->simulationTime;
}

int GlobalData::getSimulationStepTime()
{
	return this->simulationStepTime;
}

double GlobalData::getConductivity()
{
	return this->conductivity;
}

double GlobalData::getAlfa()
{
	return this->alfa;
}

double GlobalData::getTot()
{
	return this->tot;
}

double GlobalData::getInitialTemp()
{
	return this->initialTemp;
}

double GlobalData::getDensity()
{
	return this->density;
}

double GlobalData::getSpecificHeat()
{
	return this->specificHeat;
}
GlobalData::~GlobalData() {}