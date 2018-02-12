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

double Grid::minmod(int cell)
{
	double a = (gridArray[cell]-gridArray[cell-1])/cellWidth;
	double b = (gridArray[cell+1]-gridArray[cell])/cellWidth;

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
int main()
{
	Grid grid1(10, -1, 1);
	return 0;
}


