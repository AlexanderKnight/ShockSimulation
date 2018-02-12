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
	else if (std::abs(a) > std::abs(b) && (a*b) > 0)
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
		FieldArray(int cells, double start, double end);
	private:
		std::vector<double> field={};
		int cellNum;
		double startValue;
		double endValue;
		double cellWidth;
};

FieldArray::FieldArray(int cells, double start, double end, double (*func)())
{
	field.resize(cells);
}

class Grid
{
	public:
		Grid( int cells, double start, double end);


	private:
		std::vector<double> gridArray={};
		std::vector<double> gridCoor = {};
		int cellNum;
		double startValue;
		double endValue;
		double cellWidth;
	
};

Grid::Grid(int cells, double start, double end)
{
	gridArray.resize(cells);
	gridCoor.resize(cells);

	startValue = start;
	endValue = end;
	cellWidth = (end-start)/((double) cells);

	for (int i = 0; i <cells; i++){
		gridCoor[i] = start + cellWidth/2 + i*cellWidth;
		gridArray[i] = func(gridCoor[i]);
	}
}

int main()
{
	Grid grid1(10, -1, 1);
	return 0;
}


