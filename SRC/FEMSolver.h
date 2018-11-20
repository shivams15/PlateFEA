#ifndef __FEMSOLVER_H
#define __FEMSOLVER_H
#include "petsc.h"
#include "Grid.h"

class FEMSolver{
private:
	Vec u, sxx, syy, sxy;	//Solution vectors for displacements and stresses
	Vec f, fxx, fyy, fxy;	//Nodal force vectors
	Mat K, L;				//Stiffness matrices 
	KSP uSolver;			//KSP solver for displacements
	KSP sSolver;			//KSP solver for stresses
	Grid *grid;
	double nu;				//Poisson's ratio
	double E;				//Young's modulus
	void SolverSetup();
	void ExportData();
	void AssembleStiffnessMatrix();
	void SetupStressSolver();
	void CalculateNodalStresses();
	double* N_x(double xi, double eta, int* nodes, double* N0_xi, double* N0_eta);
	double* N_y(double xi, double eta, int* nodes, double* N0_xi, double* N0_eta);
	double x_xi(double xi, double eta, int* nodes, double* N0_xi);
	double x_eta(double xi, double eta, int* nodes, double* N0_eta);
	double y_xi(double xi, double eta, int* nodes, double* N0_xi);
	double y_eta(double xi, double eta, int* nodes, double* N0_eta);
	double J(double xi, double eta, int* nodes, double* N0_xi, double* N0_eta);
	double J_abs(double xi, double eta, int* nodes, double* N0_xi, double* N0_eta);
	void EvaluateAtGaussPoints(double** x, double* (*f)(double, double));
	void EvaluateAtGaussPoints(double* x, int* nodes, double** N_Xi, double** N_Eta, double (FEMSolver::*f)(double, double, int*, double*, double*));
	void EvaluateAtGaussPoints(double** x, int* nodes, double** N_Xi, double** N_Eta, double* (FEMSolver::*f)(double, double, int*, double*, double*));
	static double* N(double xi, double eta);
	static double* N_xi(double xi, double eta);
	static double* N_eta(double xi, double eta);
	static double GaussIntegralAB(double* J0, double** A, double** B, int j, int k);
	static double ShapeFunctionIntegral_xi(double eta, int i);
	static double ShapeFunctionIntegral_eta(double xi, int i);
public:
	FEMSolver(double nu, double E, Grid* grid);
	void Solve();
};

#endif
