#include "Grid.h"
#include <cmath>
#include <iostream>

using namespace std;

Grid::Grid(int nr, int na, double w){
	this->nr = nr;
	this->na = na;
	this->w = w;
	l = 2*w;
	SetupGrid();
	printf("Grid generated.\n");
}

//Generates grid points and element list based on user input from the command line
void Grid::SetupGrid(){
	n = (4*na+1)*(2*nr+1) + 2*na*(2*na+1);
	m = 2*na*nr + na*na;
	pts = new Point[n];
	eList = new Element[m];
	double da = M_PI/(4*na);
	double dr = 0;

	for(int i = 0; i<=na; i++){
		dr = (w/cos(i*da) - 1)/nr;
		for(int j = 0; j<=2*nr; j++){
			pts[(2*i)*(2*nr+1)+j].y = (0.5*j*dr+1)*cos(i*da);
			pts[(2*i)*(2*nr+1)+j].x = (0.5*j*dr+1)*sin(i*da);
			pts[(4*na-2*i)*(2*nr+1)+j].x = pts[(2*i)*(2*nr+1)+j].y;
			pts[(4*na-2*i)*(2*nr+1)+j].y = pts[(2*i)*(2*nr+1)+j].x;
			if(i > 0){
				pts[(2*i-1)*(2*nr+1)+j].y = 0.5*(pts[(2*i)*(2*nr+1)+j].y + pts[(2*i-2)*(2*nr+1)+j].y);
				pts[(2*i-1)*(2*nr+1)+j].x = 0.5*(pts[(2*i)*(2*nr+1)+j].x + pts[(2*i-2)*(2*nr+1)+j].x);
				pts[(4*na-(2*i-1))*(2*nr+1)+j].x = pts[(2*i-1)*(2*nr+1)+j].y;
				pts[(4*na-(2*i-1))*(2*nr+1)+j].y = pts[(2*i-1)*(2*nr+1)+j].x;
			}
		}
	}

	for(int i = 0; i<2*na; i++){
		for(int j = 0; j<nr; j++){
			eList[nr*i+j].nodes[0] = (2*i)*(2*nr+1) + 2*j;
			eList[nr*i+j].nodes[1] = (2*i+1)*(2*nr+1) + 2*j;
			eList[nr*i+j].nodes[2] = (2*i+2)*(2*nr+1) + 2*j;
			eList[nr*i+j].nodes[3] = (2*i)*(2*nr+1) + 2*j + 1;
			eList[nr*i+j].nodes[4] = (2*i+1)*(2*nr+1) + 2*j + 1;
			eList[nr*i+j].nodes[5] = (2*i+2)*(2*nr+1) + 2*j + 1;
			eList[nr*i+j].nodes[6] = (2*i)*(2*nr+1) + 2*j + 2;
			eList[nr*i+j].nodes[7] = (2*i+1)*(2*nr+1) + 2*j + 2;
			eList[nr*i+j].nodes[8] = (2*i+2)*(2*nr+1) + 2*j + 2;
			if(i==0)
				eList[nr*i+j].boundary[0] = BOUNDARY_LEFT;
			if(i==2*na-1)
				eList[nr*i+j].boundary[1] = BOUNDARY_BOTTOM;
			if(j==0)
				eList[nr*i+j].boundary[2] = BOUNDARY_HOLE;
			if(i<na && j==nr-1)
				eList[nr*i+j].boundary[3] = BOUNDARY_TOP;
		}
	}

	int N0 = (4*na+1)*(2*nr+1);
	double dl = (l-w)/na;
	for(int i = 0; i<=2*na; i++){
		for(int j = 0; j<2*na; j++){
			pts[N0+i*(2*na)+j].x = w + 0.5*(j+1)*dl;
			pts[N0+i*(2*na)+j].y = pts[(2*na+i+1)*(2*nr+1)-1].y;
		}
	}

	int n0 = 2*na*nr;
	for(int i = 0; i<na; i++){
		for(int j = 0; j<na; j++){
			if(j==0){
				eList[n0+na*i+j].nodes[6] = (2*na+2*i+1)*(2*nr+1)-1;
				eList[n0+na*i+j].nodes[7] = N0+(2*i)*(2*na);
				eList[n0+na*i+j].nodes[8] = N0+(2*i)*(2*na) + 1;
				eList[n0+na*i+j].nodes[3] = (2*na+2*i+2)*(2*nr+1)-1;
				eList[n0+na*i+j].nodes[4] = N0+(2*i+1)*(2*na);
				eList[n0+na*i+j].nodes[5] = N0+(2*i+1)*(2*na) + 1;
				eList[n0+na*i+j].nodes[0] = (2*na+2*i+3)*(2*nr+1)-1;
				eList[n0+na*i+j].nodes[1] = N0+(2*i+2)*(2*na);
				eList[n0+na*i+j].nodes[2] = N0+(2*i+2)*(2*na) + 1;
			}
			else{
				eList[n0+na*i+j].nodes[6] = N0+(2*i)*(2*na)+2*j-1;
				eList[n0+na*i+j].nodes[7] = N0+(2*i)*(2*na)+2*j;
				eList[n0+na*i+j].nodes[8] = N0+(2*i)*(2*na)+2*j+1;
				eList[n0+na*i+j].nodes[3] = N0+(2*i+1)*(2*na)+2*j-1;
				eList[n0+na*i+j].nodes[4] = N0+(2*i+1)*(2*na)+2*j;
				eList[n0+na*i+j].nodes[5] = N0+(2*i+1)*(2*na)+2*j+1;
				eList[n0+na*i+j].nodes[0] = N0+(2*i+2)*(2*na)+2*j-1;
				eList[n0+na*i+j].nodes[1] = N0+(2*i+2)*(2*na)+2*j;
				eList[n0+na*i+j].nodes[2] = N0+(2*i+2)*(2*na)+2*j+1;
			}
			if(i==0)
				eList[n0+na*i+j].boundary[3] = BOUNDARY_TOP;
			if(i==na-1)
				eList[n0+na*i+j].boundary[2] = BOUNDARY_BOTTOM;
			if(j==na-1)
				eList[n0+na*i+j].boundary[1] = BOUNDARY_RIGHT;
		}
	}
}