#include <iostream>
using namespace std;
#include <vector>
#include <math.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <fstream>

double initFunc (double x)
// Sets initial value for y(x)
{
	double val = exp (-1*(double)16*x*x);
	return val;
}

double func (double y)
// Function that defines Burgers' equation
{
	double val = (y*y)/2.0;
	return val;
}

double llf(double yLeft, double yRight, double velocity)
// Gives Local Lax Friedrich value for boundary, taking in values from either side
{
	double funcValLeft = func(yLeft);
	double funcValRight = func(yRight);

	double lffVal = (funcValLeft+funcValRight)/2.0 - (velocity/2.0)*(yRight-yLeft);

	return lffVal;
}

double deltaVal(double valLeft, double valRight, double deltaStep)
// Gives a simple derivative of a linear slope, eg (delta y)/(delta x)
{
	double deriv = (valRight+valLeft)/(deltaStep);

	return deriv;
}

// Runge-Kutta funcitons, using second order and midpoint 
double rk1 (double val, double deltaVal, double deltaTime)
{
	double halfStepVal = val - deltaVal*(deltaTime/2.0);
	return halfStepVal;
}

double rk2 (double val, double deltaHalfStepVal, double deltaTime)
{
	double nextVal = val - deltaHalfStepVal*deltaTime;
	return nextVal;
}

double euler ( double val, double slope, double deltaTime)
{
	double nextVal = val + slope*deltaTime;
}

double minmod (double preVal, double centVal, double nextVal, double cellWidth)
//uses minmod numerical flux limiter to determine slope of cell from neighbor cells
{

	double a = (centVal-preVal)/cellWidth;
	double b = (nextVal-centVal)/cellWidth;

	if (std::abs(a) < std::abs(b) && (a*b)> (double) 0)
	{
		return a;
	}
	else if (std::abs(a) > std::abs(b) && (a*b) > (double) 0)
	{
		return b;
	}
	else
	{
		return (double) 0;
	}
}

double fromm ( double preVal, double postVal, double cellWidth)
// Calculates using center slope, using Fromm's method
{
	double val = (postVal-preVal)/(2.0*cellWidth);
	return val;
}

double beamWarming ( double preVal, double currentVal, double cellWidth)
// Calculates slope using upwind slope, Beam-Warming method
{
	double val = (currentVal-preVal)/(cellWidth);
	return val;
}

double laxWendroff ( double currentVal, double postVal, double cellWidth)
// Calculates slope using downwind slope, Lax-Wendroff method
{
	double val = (postVal-currentVal)/(cellWidth);
}


class FieldArray
// Manages vectors for the Grid class below. Still needs to be refined.
{
	public:
		FieldArray(int cells);
		vector<double> getField();
		double getFieldValue(int i);
		void setFieldValue(int i, double x);
	private:
		std::vector<double> field={};
		int cellNum;
};

FieldArray::FieldArray(int cells)
{
	field.resize(cells);
	for (int i=0; i <cells;i++){
		field[i]=0;
	}
}

vector<double> FieldArray::getField()
{
	return field;
}

double FieldArray::getFieldValue(int i)
{
	return field[i];
}

void FieldArray::setFieldValue(int i, double x)
{
	field[i]=x;
}


class Grid
// Contains information needed to do a Berger's equation shock simulation
// For each cell, coordinate value and field value are for center of cell.
// Slope of cell passes through field value at center of cell.
// Uses 2nd order Runge-Kutta, Courant-Friedrichs-Lewy condition, minmod,
// and Local Lax Friedrich method to get simulation
{
	public:
		Grid( int cells, double start, double end, std::string slopeMethod);
		void updateSlopes(std::string, int n);
		void updateCellEdgeFields(int n);
		void updateMaxSpeed(int n);
		void timeStep(std::string, std::string);
		double getField(int n);

	private:
		int cellNum;
		double startValue;
		double endValue;
		double cellWidth;
		double maxSpeed;
		double maxSpeedHalfStep;
		FieldArray *coordinate;
		FieldArray *deltaFuncVal;

		FieldArray *field;
		FieldArray *slope;
		FieldArray *leftCellEdgeField;
		FieldArray *rightCellEdgeField;
		FieldArray *leftCellEdgeFieldFunc;
		FieldArray *rightCellEdgeFieldFunc;

		// Value for half a time step, used for runge-kutta
		FieldArray *fieldHalfStep;
		FieldArray *slopeHalfStep;
		FieldArray *leftCellEdgeFieldHalfStep;
		FieldArray *rightCellEdgeFieldHalfStep;
		FieldArray *leftCellEdgeFieldFuncHalfStep;
		FieldArray *rightCellEdgeFieldFuncHalfStep;

};

Grid::Grid(int cells, double start, double end, std::string slopeMethod)
{
	cout << "Grid Constructor Called" << endl;

	this->cellNum = cells;
	this->startValue = start;
	this->endValue = end;
	this->cellWidth = (end-start)/((double) cells);
	this->maxSpeed = 0.0;
	this->maxSpeedHalfStep = 0.0;
	
	cout << "Building new FieldArrays in Grid Constructor" << endl;

	field = new FieldArray(this->cellNum);
	coordinate = new FieldArray(this->cellNum);
	slope = new FieldArray(this->cellNum);
	leftCellEdgeField= new FieldArray(this->cellNum);
	rightCellEdgeField = new FieldArray(this->cellNum);
	leftCellEdgeFieldFunc = new FieldArray(this->cellNum);
	rightCellEdgeFieldFunc = new FieldArray(this->cellNum);
	deltaFuncVal = new FieldArray(this->cellNum);


	fieldHalfStep = new FieldArray(this->cellNum);
	slopeHalfStep = new FieldArray(this->cellNum);
	leftCellEdgeFieldHalfStep = new FieldArray(this->cellNum);
	leftCellEdgeFieldFuncHalfStep = new FieldArray(this->cellNum);
	rightCellEdgeFieldHalfStep = new FieldArray(this->cellNum);
	rightCellEdgeFieldFuncHalfStep = new FieldArray(this->cellNum);

	//field (this->cellNum);
	//coordinate (this->cellNum);
	//slope (this->cellNum);
	//leftCellEdgeField(this->cellNum);
	//rightCellEdgeField(this->cellNum);
	//deltaFuncVal(this->cellNum);


	//fieldHalfStep(this->cellNum);
	//slopeHalfStep(this->cellNum);
	//leftCellEdgeFieldHalfStep(this->cellNum);
	//rightCellEdgeFieldHalfStep(this->cellNum);
	
	cout << "Setting initial values for FieldArrays" << endl;
	for (int i = 0; i <cells; i++){
		
		coordinate->setFieldValue(i, start + this->cellWidth/2 + i*(this->cellWidth));
		field->setFieldValue(i,initFunc(coordinate->getFieldValue(i)));
		fieldHalfStep->setFieldValue(i,0);

		if (std::abs(field->getFieldValue(i))>maxSpeed){
			maxSpeed = std::abs(field->getFieldValue(i));
		}
		
	}

	cout << "Updating Slopes" << endl;
	updateSlopes(slopeMethod,1);
	updateSlopes(slopeMethod,2);
	cout << "Updating Cell Edges" << endl;
	updateCellEdgeFields(1);
	updateCellEdgeFields(2);
	cout << "Updating Max Speed" << endl;
	updateMaxSpeed(1);
	updateMaxSpeed(2);

}

double Grid::getField(int n)
{
	return field->getFieldValue(n);
}

void Grid::updateMaxSpeed(int n)
// Update max speed for field and fieldHalfStep. Can choose which ones to update:
// n=0, update both maxSpeed and maxSpeedHalfStep
// n=1, update just maxSpeed
// n=2, update just maxSpeedHalfStep
{
	if (n==0){
		for (int i = 0; i<this->cellNum;i++){
			if (std::abs(field->getFieldValue(i))>this->maxSpeed){
				maxSpeed = std::abs(field->getFieldValue(i));
			}
			if (std::abs(fieldHalfStep->getFieldValue(i))>this->maxSpeedHalfStep){
				maxSpeedHalfStep = std::abs(fieldHalfStep->getFieldValue(i));
			}
		}
	}
	else if (n==1){
		for (int i = 0; i<this->cellNum;i++){
			if (std::abs(field->getFieldValue(i))>this->maxSpeed){
				maxSpeed = std::abs(field->getFieldValue(i));
			}
		}
	}
	else if (n==2){
		for (int i =0; i<this->cellNum;i++){
			if (std::abs(fieldHalfStep->getFieldValue(i))>this->maxSpeedHalfStep){
				maxSpeedHalfStep = std::abs(fieldHalfStep->getFieldValue(i));
			}
		}
	}
}

void Grid::updateSlopes(std::string slopeMethod= "minmod", int n = 1)
// Updates slopes from cell value and neighbors. Uses minmod to calculate. 
// Grid is sinosoidal, so first and last cell are connected.
// Can choose to update slope, slopeHalfStep, or both by index n
// If n=0, update both slope and slopeHalfStep
// If n=1, update just slope
// if n=2, update just slopeHalfStep

{
	
	double preVal, centVal, postVal;
	double preHalfVal, centHalfVal, postHalfVal;

	if (slopeMethod == "minmod"){


		if (n==1){
	
			for (int i = 0; i < cellNum; i++){
				// If cell is first or last, grabs the other as neighbor cell.
				// If cell is in the middle, just updates neighbors and uses those values.
				if (i == 0){
					preVal = field->getFieldValue(cellNum-1);
					centVal = field->getFieldValue(i);
					postVal = field->getFieldValue(i+1);
					slope->setFieldValue(i, minmod(preVal,centVal,postVal,cellWidth));
				}
				else if (i == (cellNum-1)){
					preVal = centVal;
					centVal = postVal;
					postVal = field->getFieldValue(0);
					slope->setFieldValue(i, minmod(preVal, centVal, postVal, cellWidth));
				}
				else{
					preVal = centVal;
					centVal = postVal;
					postVal = field->getFieldValue(i+1);
					slope->setFieldValue(i, minmod(preVal, centVal, postVal, cellWidth));
				}
			}
		}
		else if (n==2){
	
			for (int i = 0; i < cellNum; i++){
				// If cell is first or last, grabs the other as neighbor cell.
				// If cell is in the middle, just updates neighbors and uses those values.
				if (i == 0){
					preHalfVal = fieldHalfStep->getFieldValue(cellNum-1);
					centHalfVal = fieldHalfStep->getFieldValue(i);
					postHalfVal = fieldHalfStep->getFieldValue(i+1);
					slope->setFieldValue(i, minmod(preVal,centVal,postVal,cellWidth));
				}
				else if (i == (cellNum-1)){
					preHalfVal = centHalfVal;
					centHalfVal = postHalfVal;
					postHalfVal = fieldHalfStep->getFieldValue(0);
					slopeHalfStep->setFieldValue(i, minmod(preVal, centVal, postVal, cellWidth));
				}
				else{
					preHalfVal = centHalfVal;
					centHalfVal = postHalfVal;
					postHalfVal = fieldHalfStep->getFieldValue(i+1);
					slopeHalfStep->setFieldValue(i, minmod(preVal, centVal, postVal, cellWidth));
				}
			}
		}

	}

	else if (slopeMethod == "fromm"){
		if (n==1){
			for (int i = 0; i < cellNum; i++){
				if (i == 0)
				{
					preVal = field->getFieldValue(cellNum-1);
					postVal = field->getFieldValue(i+1);
				}
	
				else if (i == (cellNum-1))
				{
					preVal = field->getFieldValue(i-1);
					postVal = field->getFieldValue(0);
				}
	
				else
				{
					preVal = field->getFieldValue(i-1);
					postVal = field->getFieldValue(i+1);
				}
	
				slope->setFieldValue(i, fromm(preVal, postVal, cellWidth));
			}
		}
		if (n==2){
			for (int i = 0; i < cellNum; i++){
				if (i == 0)
				{
					preVal = fieldHalfStep->getFieldValue(cellNum-1);
					postVal = fieldHalfStep->getFieldValue(i+1);
				}
	
				else if (i == (cellNum-1))
				{
					preVal = fieldHalfStep->getFieldValue(i-1);
					postVal = fieldHalfStep->getFieldValue(0);
				}
	
				else
				{
					preVal = fieldHalfStep->getFieldValue(i-1);
					postVal = fieldHalfStep->getFieldValue(i+1);
				}
	
				slopeHalfStep->setFieldValue(i, fromm(preVal, postVal, cellWidth));
			}
		}
	}

	else if (slopeMethod == "beamWarming"){
		if (n==1){
			for (int i =0; i< cellNum; i++)
			{
				if (i==0){
					preVal = field->getFieldValue(cellNum-1);
					centVal = field->getFieldValue(i);
				}
				else {
					preVal = centVal;
					centVal = field->getFieldValue(i);
				}
			}
		}

		if (n==2){
			for (int i =0; i< cellNum; i++)
			{
				if (i==0){
					preVal = fieldHalfStep->getFieldValue(cellNum-1);
					centVal = fieldHalfStep->getFieldValue(i);
				}
				else {
					preVal = centVal;
					centVal = fieldHalfStep->getFieldValue(i);
				}
			slopeHalfStep->setFieldValue(i, beamWarming(preVal, centVal, cellWidth));
			}
		}
	}


	else if (slopeMethod == "laxWendroff"){
		if (n==1){
			for (int i =0; i<cellNum; i++){
				if (i==(cellNum-1)){
					centVal = postVal;
					postVal =  field->getFieldValue(0);
				}
				else if (i==0){
					centVal = field->getFieldValue(i);
					postVal = field->getFieldValue(i+1);
				}
				else{
					centVal = postVal;
					postVal = field->getFieldValue(i+1);
				}
			slope->setFieldValue(i, laxWendroff(centVal, postVal, cellWidth));
			}
		}

		if (n==2){
			for (int i =0; i<cellNum; i++){
				if (i==(cellNum-1)){
					centVal = postVal;
					postVal =  fieldHalfStep->getFieldValue(0);
				}
				else if (i==0){
					centVal = fieldHalfStep->getFieldValue(i);
					postVal = fieldHalfStep->getFieldValue(i+1);
				}
				else{
					centVal = postVal;
					postVal = fieldHalfStep->getFieldValue(i+1);
				}
			slopeHalfStep->setFieldValue(i, laxWendroff(centVal, postVal, cellWidth));
			}
		}
		
	}
}

void Grid::updateCellEdgeFields(int n)
// Extrapolates cell field edge values from center value by slope
// Use index to choose to update cell edges for field or field at half time step, or both
// n=0: update both left/rightCellEdgeField and left/rightCellEdgeFieldHalfStep
// n=1: update just left/rightCellEdgeField
// n=2: update just left/rightCellEdgeFieldHalfStep
{
	if (n==0){
		for (int i =0; i<cellNum; i++){
			leftCellEdgeField->setFieldValue(i,(field->getFieldValue(i)-slope->getFieldValue(i)*(cellWidth/2)));
			rightCellEdgeField->setFieldValue(i,(field->getFieldValue(i) + slope->getFieldValue(i)*(cellWidth/2)));
			leftCellEdgeFieldFunc->setFieldValue(i, func(leftCellEdgeField->getFieldValue(i)));
			rightCellEdgeFieldFunc->setFieldValue(i, func(rightCellEdgeField->getFieldValue(i)));
			leftCellEdgeFieldHalfStep->setFieldValue(i,(fieldHalfStep->getFieldValue(i)-slopeHalfStep->getFieldValue(i)*(cellWidth/2)));
			rightCellEdgeFieldHalfStep->setFieldValue(i,(fieldHalfStep->getFieldValue(i) + slopeHalfStep->getFieldValue(i)*(cellWidth/2)));
			leftCellEdgeFieldFuncHalfStep->setFieldValue(i, func(leftCellEdgeFieldHalfStep->getFieldValue(i)));
			rightCellEdgeFieldFuncHalfStep->setFieldValue(i, func(rightCellEdgeFieldHalfStep->getFieldValue(i)));
		}
	}

	else if (n==1){
		for (int i =0; i<cellNum; i++){

			leftCellEdgeField->setFieldValue(i,(field->getFieldValue(i)-slope->getFieldValue(i)*(cellWidth/2)));
			rightCellEdgeField->setFieldValue(i,(field->getFieldValue(i) + slope->getFieldValue(i)*(cellWidth/2)));
			leftCellEdgeFieldFunc->setFieldValue(i, func(leftCellEdgeField->getFieldValue(i)));
			rightCellEdgeFieldFunc->setFieldValue(i, func(rightCellEdgeField->getFieldValue(i)));
		}
	}

	else if (n==2){
		for (int i =0; i<cellNum; i++){

			leftCellEdgeFieldHalfStep->setFieldValue(i,(fieldHalfStep->getFieldValue(i)-slopeHalfStep->getFieldValue(i)*(cellWidth/2)));
			rightCellEdgeFieldHalfStep->setFieldValue(i,(fieldHalfStep->getFieldValue(i) + slopeHalfStep->getFieldValue(i)*(cellWidth/2)));
			leftCellEdgeFieldFuncHalfStep->setFieldValue(i, func(leftCellEdgeFieldHalfStep->getFieldValue(i)));
			rightCellEdgeFieldFuncHalfStep->setFieldValue(i, func(rightCellEdgeFieldHalfStep->getFieldValue(i)));
		}
	}
}	

void Grid::timeStep(std::string slopeMethod = "minmod", std::string stepMethod = "rk")
// function to push the grid through one time step. Needs some major debugging before ready.
//
// CREATE FUNCTIONS TO CALL VALUES FROM GRID
//
// PARTICULARLY: getField, getDeltaField, setField, setFieldHalfStep, getFieldHalfStep, getDeltaFieldHalfStep
//
// God help you.
{
	double deltaT = 0.9*(this->cellWidth)/(this->maxSpeed);
	double nextLLF = 0.0;
	double prevLLF = 0.0;

	if ((slopeMethod == "minmod") && (stepMethod == "rk"))
	{

		for (int i=0; i<(this->cellNum); i++){
			if (i==0){
				prevLLF = llf(rightCellEdgeField->getFieldValue((this->cellNum)-1),leftCellEdgeField->getFieldValue(i),this->maxSpeed);
				nextLLF = llf(rightCellEdgeField->getFieldValue(i), leftCellEdgeField->getFieldValue(i+1),this->maxSpeed);
				deltaFuncVal->setFieldValue(i,((nextLLF-prevLLF)/(this->cellWidth)));
			}
			else if (i==(this->cellNum-1)){
				prevLLF = nextLLF;
				nextLLF = llf(rightCellEdgeField->getFieldValue(i), leftCellEdgeField->getFieldValue(0), this->maxSpeed);
				deltaFuncVal->setFieldValue(i,((nextLLF-prevLLF)/(this->cellWidth)));
			}
			else {
				prevLLF = nextLLF;
				nextLLF = llf(rightCellEdgeField->getFieldValue(i), leftCellEdgeField->getFieldValue(i+1),this->maxSpeed);
				deltaFuncVal->setFieldValue(i,((nextLLF-prevLLF)/(this->cellWidth)));
			}
		}
		for (int i=0; i<(this->cellNum); i++){
			fieldHalfStep->setFieldValue(i,rk1(field->getFieldValue(i),deltaFuncVal->getFieldValue(i),deltaT));
		}
		if (slopeMethod == "minmod"){
			updateSlopes("minmod",2);
		}
		else if (slopeMethod == "fromm"){
			updateSlopes("fromm",2);
		}
		else if (slopeMethod == "beamWarming"){
			updateSlopes("beamWarming",2);
		}
		else if (slopeMethod == "laxWendroff"){
			updateSlopes("laxWendroff",2);
		}

		// Setup for variable slope calculations


		updateCellEdgeFields(2);
		updateMaxSpeed(2); //Make Function to update max speed.
		deltaT = 0.9*(this->cellWidth)/(this->maxSpeedHalfStep);

		for (int i=0; i<(this->cellNum); i++){
			if (i==0){
				prevLLF = llf(rightCellEdgeFieldHalfStep->getFieldValue(this->cellNum-1),leftCellEdgeFieldHalfStep->getFieldValue(i),this->maxSpeedHalfStep);
				nextLLF = llf(rightCellEdgeFieldHalfStep->getFieldValue(i), leftCellEdgeFieldHalfStep->getFieldValue(i+1),this->maxSpeedHalfStep);
				deltaFuncVal->setFieldValue(i,((nextLLF-prevLLF)/(this->cellWidth)));
			}
			else if (i==(this->cellNum-1)){
				prevLLF = nextLLF;
				nextLLF = llf(rightCellEdgeFieldHalfStep->getFieldValue(i), leftCellEdgeFieldHalfStep->getFieldValue(0), this->maxSpeedHalfStep);
				deltaFuncVal->setFieldValue(i,((nextLLF-prevLLF)/(this->cellWidth)));
			}
			else {
				prevLLF = nextLLF;
				nextLLF = llf(rightCellEdgeFieldHalfStep->getFieldValue(i), leftCellEdgeFieldHalfStep->getFieldValue(i+1),this->maxSpeedHalfStep);
				deltaFuncVal->setFieldValue(i,((nextLLF-prevLLF)/(this->cellWidth)));
			}
		}
		for (int i=0; i<(this->cellNum); i++){
			field->setFieldValue(i,rk2(field->getFieldValue(i),deltaFuncVal->getFieldValue(i),deltaT));
		}
	}

	else if (stepMethod == "euler")
	{
		for (int i = 0; i<(this->cellNum); i++){
			if (i==0){
				prevLLF = llf(rightCellEdgeField->getFieldValue(this->cellNum-1),leftCellEdgeField->getFieldValue(i),this->maxSpeed);
				nextLLF = llf(rightCellEdgeField->getFieldValue(i), leftCellEdgeField->getFieldValue(i+1), this->maxSpeed);
				deltaFuncVal->setFieldValue(i,((nextLLF-prevLLF)/(this->cellWidth)));
			}
			else if (i==(this->cellNum-1)){
				prevLLF = nextLLF;
			nextLLF = llf(rightCellEdgeField->getFieldValue(i),leftCellEdgeField->getFieldValue(0), this->maxSpeed);
			deltaFuncVal->setFieldValue(i,((nextLLF-prevLLF)/(this->cellWidth)));
			}
			else {
				prevLLF = nextLLF;
				nextLLF = llf(rightCellEdgeField->getFieldValue(i),leftCellEdgeField->getFieldValue(i+1),this->maxSpeed);
				deltaFuncVal->setFieldValue(i,((nextLLF-prevLLF)/(this->cellWidth)));
			}
		}
		for (int i=0; i<(this->cellNum); i++){
			field->setFieldValue(i, euler(field->getFieldValue(i), -1.0*(deltaFuncVal->getFieldValue(i)), deltaT));
		}
	}
	if (slopeMethod == "minmod"){
		updateSlopes("minmod",1);
	}

	else if (slopeMethod == "fromm"){
		updateSlopes("fromm");
	}
	else if (slopeMethod == "beamWarming"){
		updateSlopes("beamWarming");
	}
	else if (slopeMethod == "laxWendroff"){
		updateSlopes("laxWendroff");
	}



		updateCellEdgeFields(1);
		updateMaxSpeed(1);
	
}
int main()
{
	cout << "Main initialized" << endl;
	Grid grid (1000,-1.0,1.0, "minmod");
	ofstream testData ("testData.dat");
	std::string slopeStyle = "minmod";
	std::string stepStyle = "rk";
	for (int i = 0; i < 10000; i++){
		grid.timeStep(slopeStyle,stepStyle);
		if (testData.is_open()){
			for (int j=0; j <1000; j++){
				testData << grid.getField(j);
				if (j != 1000-1){
					testData << ',';
				}
			}
			testData << "\n";
		}
	}
	testData.close();
	return 0;
}

// Set up measurement of error and convergence between different methods
