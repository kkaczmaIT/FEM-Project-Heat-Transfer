#pragma once
#include<iostream>

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
	GlobalData();
	GlobalData(int simulationTime, int simulationStepTime, double conductivity, double alfa, double tot, double initialTemp, double density, double specificHeat);

	// setters

	void setSimulationTime(int simulationTime);

	void setSimulationStepTime(int simulationStepTime);

	void setConductivity(double conductivity);

	void setAlfa(double alfa);

	void setTot(double tot);

	void setInitialTemp(double initialTemp);

	void setDensity(double density);

	void setSpecificHeat(double specificHeat);

	// getters
	int getSimulationTime();

	int getSimulationStepTime();

	double getConductivity();

	double getAlfa();

	double getTot();

	double getInitialTemp();

	double getDensity();

	double getSpecificHeat();
	~GlobalData();
};