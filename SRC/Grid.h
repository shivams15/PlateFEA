#ifndef __GRID_H
#define __GRID_H
#define INTERNAL 0
#define BOUNDARY_LEFT 1
#define BOUNDARY_RIGHT 2
#define BOUNDARY_BOTTOM 3
#define BOUNDARY_TOP 4
#define BOUNDARY_HOLE 5

struct Point{
	double x;
	double y;
};

struct Element{
	int nodes[9];
	int boundary[4] = {INTERNAL, INTERNAL, INTERNAL, INTERNAL};
};

class Grid{
private:
	double l;		//length of plate
	double w ;		//plate width		
	int nr;			//number of elements in radial direction
	int na;			//half of the number of elements in the angular direction
	void SetupGrid();
public:
	int m;			//total number of elements
	int n;			//total number of nodes
	Point *pts;
	Element *eList;
	Grid(int nr, int na, double w);
};

#endif