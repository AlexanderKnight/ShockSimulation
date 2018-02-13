#include <iostream>
using namespace std;
#include <vector>
#include <math.h>
#include <algorithm>
#include <cmath>

double func (double x)
// Sets initial value for y(x)
{
	double val = exp (-1*(double)16*x*x);
	return val;
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
		Grid( int cells, double start, double end);
		void updateSlopes();
		void updateCellEdgeFields();

	private:
		int cellNum;
		double startValue;
		double endValue;
		double cellWidth;
		FieldArray *field;
		FieldArray *coordinate;
		FieldArray *slope;
		FieldArray *leftCellEdgeField;
		FieldArray *rightCellEdgeField;
	
};

Grid::Grid(int cells, double start, double end)
{
	cellNum = cells;
	startValue = start;
	endValue = end;
	cellWidth = (end-start)/((double) cells);
		
	FieldArray field (cellNum);
	FieldArray coordinate (cellNum);
	FieldArray slope (cellNum);
	FieldArray leftCellEdgeField(cellNum);
	FieldArray rightCellEdgeField(cellNum);


	for (int i = 0; i <cells; i++){
		coordinate.setFieldValue(i, start + cellWidth/2 + i*cellWidth);
		field.setFieldValue(i,func(coordinate.getFieldValue(i)));
		
	}

	updateSlopes();
	updateCellEdgeFields();

}

void Grid::updateSlopes()
// Updates slopes from cell value and neighbors. Uses minmod to calculate. 
// Grid is sinosoidal, so first and last cell are connected.

{
	
	double preVal, centVal, postVal;
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

void Grid::updateCellEdgeFields()
// Extrapolates cell field edge values from center value by slope
{
	for (int i =0; i<cellNum; i++){

		leftCellEdgeField->setFieldValue(i,(field->getFieldValue(i)-slope->getFieldValue(i)*(cellWidth/2)));
		rightCellEdgeField->setFieldValue(i,(field->getFieldValue(i) + slope->getFieldValue(i)*(cellWidth/2)));
	}
}	

int main()
{
	Grid grid1(10, -1, 1);
	return 0;
}


