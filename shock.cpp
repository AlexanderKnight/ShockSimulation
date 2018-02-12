#include <iostream>
using namespace std;
#include <vector>
#include <math.h>
#include <algorithm>
#include <cmath>

double func (double x)
{
	double val = exp (-1*(double)16*x*x);
	return val;
}

double minmod (double preVal, double centVal, double nextVal, double cellWidth)
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
{
	public:
		FieldArray(int cells);
		getField();
		getFieldValue(int i);
		setFieldValue(int i, double x);
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
{
	public:
		Grid( int cells, double start, double end);
		updateSlopes();
		updateCellEdgeFields();

	private:
		//std::vector<double> gridArray={};
		//std::vector<double> gridCoor = {};
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
	
	//gridArray.resize(cells);
	//gridCoor.resize(cells);
	
	FieldArray field (cellNum);
	FieldArray coordinate (cellNum);
	FieldArray slope (cellNum);
	FieldArray leftCellEdgeField(cellNum);
	FieldArray rightCellEdgeField(cellNum);


	for (int i = 0; i <cells; i++){
		coordinate.setFieldValue(i, start + cellWidth/2 + i*cellWidth);
		field.setFieldValue(i,func(gridCoor[i]));
	}

	updateSlopes();
	updateCellEdgeFields();

}

void Grid::updateSlopes()
{
	
	double preVal, centVal, postVal;
	for (int i = 0; i < cellNum; i++){
		if (i == 0){
			preVal = field.getFieldValue(cellNum-1);
			centVal = field.getFieldValue(i);
			postVal = field.getFieldValue(i+1);
			slope.setFieldValue(i, minmod(preVal,centVal,postVal,cellWidth);
		}
		else if (i == (cellNum-1)){
			preVal = centVal;
			centVal = postVal;
			postVal = field.getFieldValue(0);
			slope.setFieldValue(i, minmod(preVal, centVal, postVal, cellwidth);
		}
		else{
			preVal = centVal;
			centVal = postVal;
			postVal = field.getFieldValue(i+1);
			slope.setFieldValue(i, minmod(preVal, centVal, postVal, cellwidth);
		}
	}
}

void Grid::updateCellEdgeFields()
{
	for (int i =0; i<cellNum; i++){

		leftCellEdgeField.setFieldValue(i,(field.getFieldValue(i)-slope.getFieldValue(i)*(cellWidth/2)));
		rightCellEdgeField.setFieldValue(i,(field.getFieldValue(i) + slope.getFieldValue(i)*(cellwidth/2)));
	}
}	

int main()
{
	Grid grid1(10, -1, 1);
	return 0;
}


