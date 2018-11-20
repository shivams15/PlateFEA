#include "FEMSolver.h"
#include "Grid.h"
#include "petsc.h"
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

FEMSolver::FEMSolver(double nu, double E, Grid *grid){
	this->nu = nu;
	this->E  = E;
	this->grid = grid;

	MatCreate(PETSC_COMM_SELF, &K);
	MatSetType(K, MATAIJ);
	MatSetSizes(K, PETSC_DECIDE, PETSC_DECIDE, 2*grid->n, 2*grid->n);
	MatSetFromOptions(K);
	MatSetUp(K);

	MatCreate(PETSC_COMM_SELF, &L);
	MatSetType(L, MATAIJ);
	MatSetSizes(L, PETSC_DECIDE, PETSC_DECIDE, grid->n, grid->n);
	MatSetFromOptions(L);
	MatSetUp(L);

	VecCreateSeq(PETSC_COMM_SELF, 2*grid->n, &f);
	VecCreateSeq(PETSC_COMM_SELF, 2*grid->n, &u);
	VecCreateSeq(PETSC_COMM_SELF, grid->n, &sxx);
	VecCreateSeq(PETSC_COMM_SELF, grid->n, &fxx);
	VecCreateSeq(PETSC_COMM_SELF, grid->n, &sxy);
	VecCreateSeq(PETSC_COMM_SELF, grid->n, &fxy);
	VecCreateSeq(PETSC_COMM_SELF, grid->n, &syy);
	VecCreateSeq(PETSC_COMM_SELF, grid->n, &fyy);
	VecSet(f,0.0);
	VecSet(u,0.0);
	VecSet(fxx,0.0);
	VecSet(sxx,0.0);
	VecSet(fxy,0.0);
	VecSet(sxy,0.0);
	VecSet(fyy,0.0);
	VecSet(syy,0.0);

	SolverSetup();
	cout<<"Solver Setup Complete.\n";
}

void FEMSolver::Solve(){
	KSPSolve(uSolver,f,u);
	CalculateNodalStresses();
	cout<<"Solve Complete.\n";
	ExportData();
}

void FEMSolver::ExportData(){
	cout<<"Writing Data.\n";
	FILE *f1 = fopen("u.csv","w");
	PetscScalar *u0 = new PetscScalar[2*grid->n];
	VecGetArray(u, &u0);

	for(int i = 0; i<grid->n; i++){
		fprintf(f1, "%lf,%lf,%lf\n", grid->pts[i].x, grid->pts[i].y, u0[2*i]);
	}
	fclose(f1);

	f1 = fopen("v.csv","w");
	for(int i = 0; i<grid->n; i++){
		fprintf(f1, "%lf,%lf,%lf\n", grid->pts[i].x, grid->pts[i].y, u0[2*i+1]);
	}
	fclose(f1);

	u0 = new PetscScalar[grid->n];
	f1 = fopen("sxx.csv","w");
	VecGetArray(sxx, &u0);
	for(int i = 0; i<grid->n; i++){
		fprintf(f1, "%lf,%lf,%lf\n", grid->pts[i].x, grid->pts[i].y, u0[i]);
	}
	fclose(f1);

	f1 = fopen("syy.csv","w");
	VecGetArray(syy, &u0);
	for(int i = 0; i<grid->n; i++){
		fprintf(f1, "%lf,%lf,%lf\n", grid->pts[i].x, grid->pts[i].y, u0[i]);
	}
	fclose(f1);

	f1 = fopen("sxy.csv","w");
	VecGetArray(sxy, &u0);
	for(int i = 0; i<grid->n; i++){
		fprintf(f1, "%lf,%lf,%lf\n", grid->pts[i].x, grid->pts[i].y, u0[i]);
	}
	fclose(f1);
}

void FEMSolver::SolverSetup(){
	AssembleStiffnessMatrix();
	PC pc;
	KSPCreate(PETSC_COMM_SELF, &uSolver);
	KSPSetOperators(uSolver, K, K);
	KSPSetType(uSolver, KSPBCGSL);
	KSPGetPC(uSolver, &pc);
	PCSetType(pc, PCBJACOBI);
	KSPSetTolerances(uSolver, 1.e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
	KSPSetUp(uSolver);
}

//Calculates the stresses at the element nodes
void FEMSolver::CalculateNodalStresses(){
	SetupStressSolver();
	KSPSolve(sSolver,fxx,sxx);
	KSPSolve(sSolver,fyy,syy);
	KSPSolve(sSolver,fxy,sxy);
}

/*
Evaluates element stiffness matrices and assembles them into a global stiffness matrix
Also enforces the appropriate boundary conditions
*/
void FEMSolver::AssembleStiffnessMatrix(){
	printf("Assembling Stiffness Matrix.\n");
	double r0 = 1.0/sqrt(3);
	double r = sqrt(0.6);
	double d[4];
	double J0[9];
	double **N_X = new double*[9];
	double **N_Y = new double*[9];
	double **N_Xi = new double*[9];
	double **N_Eta = new double*[9];
	int l1, l2;
	double w;
	EvaluateAtGaussPoints(N_Xi, N_xi);
	EvaluateAtGaussPoints(N_Eta, N_eta);

	for(int i=0; i<grid->m; i++){
		EvaluateAtGaussPoints(J0, grid->eList[i].nodes, N_Xi, N_Eta, &FEMSolver::J_abs);
		EvaluateAtGaussPoints(N_X, grid->eList[i].nodes, N_Xi, N_Eta, &FEMSolver::N_x);
		EvaluateAtGaussPoints(N_Y, grid->eList[i].nodes, N_Xi, N_Eta, &FEMSolver::N_y);

		for(int j=0; j<9; j++){
			//Evaluating the nodal forces for nodes located on the boundary where the external stress is applied
			if(j%3==2 && grid->eList[i].boundary[1]==BOUNDARY_RIGHT){
				l1 = 2*grid->eList[i].nodes[j];
				w = ShapeFunctionIntegral_eta(1, j);
				w = 0.5*w*((grid->pts[grid->eList[i].nodes[8]]).y - (grid->pts[grid->eList[i].nodes[2]]).y);
				VecSetValues(f, 1, &l1, &w, ADD_VALUES);
			}
			for(int k=0; k<9; k++){
				//Calculating the element stiffness matrices
				d[0] = GaussIntegralAB(J0, N_X, N_X, j, k) + 0.5*(1-nu)*GaussIntegralAB(J0, N_Y, N_Y, j, k);
				d[1] = nu*GaussIntegralAB(J0, N_X, N_Y, j, k) + 0.5*(1-nu)*GaussIntegralAB(J0, N_Y, N_X, j, k);
				d[2] = nu*GaussIntegralAB(J0, N_Y, N_X, j, k) + 0.5*(1-nu)*GaussIntegralAB(J0, N_X, N_Y, j, k);
				d[3] = GaussIntegralAB(J0, N_Y, N_Y, j, k) + 0.5*(1-nu)*GaussIntegralAB(J0, N_X, N_X, j, k);
				
				//Performing matrix assembly
				l1 = 2*grid->eList[i].nodes[j]; l2 = 2*grid->eList[i].nodes[k];
				MatSetValues(K,1,&l1,1,&l2,(PetscScalar *)&d[0],ADD_VALUES);
				l2 = l2+1;
				MatSetValues(K,1,&l1,1,&l2,(PetscScalar *)&d[1],ADD_VALUES);
				l1 = l1+1; l2 = l2-1;
				MatSetValues(K,1,&l1,1,&l2,(PetscScalar *)&d[2],ADD_VALUES);
				l2 = l2+1;
				MatSetValues(K,1,&l1,1,&l2,(PetscScalar *)&d[3],ADD_VALUES);	
			}
		}
	}

	MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);

	//Setting up the appropriate boundary conditions at the left and bottom boundaries
	for(int i=0; i<grid->m; i++){
		for(int j=0; j<9; j++){
			if(j%3==0 && grid->eList[i].boundary[0]==BOUNDARY_LEFT)
				l1 = 2*grid->eList[i].nodes[j]; 
			else if(j%3==2 && grid->eList[i].boundary[1]==BOUNDARY_BOTTOM)
				l1 = 2*grid->eList[i].nodes[j] + 1;	
			else if(j<3 && grid->eList[i].boundary[2]==BOUNDARY_BOTTOM)
				l1 = 2*grid->eList[i].nodes[j] + 1;
			else continue;
			MatZeroRows(K,1,&l1,1.0,NULL,NULL);
		}	
	}

	MatScale(K, E/(1-pow(nu,2)));
	MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
}

//Sets up the system of equations for calculating stress values at the element nodes
void FEMSolver::SetupStressSolver(){
	double r = sqrt(0.6);
	double J0[9];
	double d[2];
	double **N0 = new double*[9];
	double **N_X = new double*[9];
	double **N_Y = new double*[9];
	double **N_Xi = new double*[9];
	double **N_Eta = new double*[9];
	double w0, w1, w2, w3;
	PetscScalar *u0 = new PetscScalar[2*grid->n];
	VecGetArray(u, &u0);

	EvaluateAtGaussPoints(N0, N);
	EvaluateAtGaussPoints(N_Xi, N_xi);
	EvaluateAtGaussPoints(N_Eta, N_eta);

	for(int i=0; i<grid->m; i++){
		EvaluateAtGaussPoints(J0, grid->eList[i].nodes, N_Xi, N_Eta, &FEMSolver::J_abs);
		EvaluateAtGaussPoints(N_X, grid->eList[i].nodes, N_Xi, N_Eta, &FEMSolver::N_x);
		EvaluateAtGaussPoints(N_Y, grid->eList[i].nodes, N_Xi, N_Eta, &FEMSolver::N_y);

		for(int j=0; j<9; j++){
			w1 = 0;
			w2 = 0;
			w3 = 0;
			for(int k=0; k<9; k++){
				w0 = GaussIntegralAB(J0, N0, N0, j, k);
				MatSetValues(L,1,&grid->eList[i].nodes[j],1,&grid->eList[i].nodes[k],(PetscScalar *)&w0,ADD_VALUES);
				d[0] = GaussIntegralAB(J0, N0, N_X, j, k);
				d[1] = GaussIntegralAB(J0, N0, N_Y, j, k);
				w1 = w1 + d[0]*u0[2*grid->eList[i].nodes[k]] + nu*d[1]*u0[2*grid->eList[i].nodes[k]+1];
				w2 = w2 + nu*d[0]*u0[2*grid->eList[i].nodes[k]] + d[1]*u0[2*grid->eList[i].nodes[k]+1];
				w3 = w3 + 0.5*(1-nu)*(d[1]*u0[2*grid->eList[i].nodes[k]] + d[0]*u0[2*grid->eList[i].nodes[k]+1]);
			}
			VecSetValues(fxx, 1, &grid->eList[i].nodes[j], &w1, ADD_VALUES);
			VecSetValues(fyy, 1, &grid->eList[i].nodes[j], &w2, ADD_VALUES);
			VecSetValues(fxy, 1, &grid->eList[i].nodes[j], &w3, ADD_VALUES);
		}
	}
	VecScale(fxx, E/(1-pow(nu,2)));
	VecScale(fxy, E/(1-pow(nu,2)));
	VecScale(fyy, E/(1-pow(nu,2)));
	MatAssemblyBegin(L, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(L, MAT_FINAL_ASSEMBLY);

	PC pc;
	KSPCreate(PETSC_COMM_SELF, &sSolver);
	KSPSetOperators(sSolver, L, L);
	KSPSetType(sSolver, KSPBCGSL);
	KSPGetPC(sSolver, &pc);
	PCSetType(pc, PCBJACOBI);
	KSPSetTolerances(sSolver, 1.e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
	KSPSetUp(sSolver);
}

//Evaluates a function at the Gaussian quadrature points
void FEMSolver::EvaluateAtGaussPoints(double* x, int* nodes, double** N_Xi, double** N_Eta, double (FEMSolver::*f)(double, double, int*, double*, double*)){
	double r = sqrt(0.6);
	for(int i=0; i<9; i++){
		x[i] = (this->*f)((i%3-1)*r, (i/3-1)*r, nodes, N_Xi[i], N_Eta[i]);
	}
}

//Evaluates a function at the Gaussian quadrature points
void FEMSolver::EvaluateAtGaussPoints(double** x, int* nodes, double** N_Xi, double** N_Eta, double* (FEMSolver::*f)(double, double, int*, double*, double*)){
	double r = sqrt(0.6);
	for(int i=0; i<9; i++){
		x[i] = (this->*f)((i%3-1)*r, (i/3-1)*r, nodes, N_Xi[i], N_Eta[i]);
	}
}

//Evaluates a function at the Gaussian quadrature points
void FEMSolver::EvaluateAtGaussPoints(double** x, double* (*f)(double, double)){
	double r = sqrt(0.6);
	for(int i=0; i<9; i++){
		x[i] = (*f)((i%3-1)*r, (i/3-1)*r);
	}
}

//Calculates shape function values at a specified location
double* FEMSolver::N(double xi, double eta){
	double *N0 = new double[9];
	N0[0] = 0.25*xi*eta*(xi-1)*(eta-1);
	N0[1] = 0.5*eta*(1-pow(xi,2))*(eta-1);
	N0[2] = 0.25*xi*eta*(xi+1)*(eta-1);
	N0[3] = 0.5*xi*(xi-1)*(1-pow(eta,2));
	N0[4] = (1-pow(xi,2))*(1-pow(eta,2));
	N0[5] = 0.5*xi*(xi+1)*(1-pow(eta,2));
	N0[6] = 0.25*xi*eta*(xi-1)*(eta+1);
	N0[7] = 0.5*eta*(1-pow(xi,2))*(eta+1);
	N0[8] = 0.25*xi*eta*(xi+1)*(eta+1);
	return N0;
}

//Calculates the derivative of shape functions with repsect to parameter xi
double* FEMSolver::N_xi(double xi, double eta){
	double *N0_xi = new double[9];
	N0_xi[0] = 0.25*eta*(eta-1)*(2*xi-1);
	N0_xi[1] = -xi*eta*(eta-1);
	N0_xi[2] = 0.25*eta*(eta-1)*(2*xi+1);
	N0_xi[3] = 0.5*(1-pow(eta,2))*(2*xi-1);
	N0_xi[4] = -2*xi*(1-pow(eta,2));
	N0_xi[5] = 0.5*(1-pow(eta,2))*(2*xi+1);
	N0_xi[6] = 0.25*eta*(eta+1)*(2*xi-1);
	N0_xi[7] = -xi*eta*(eta+1);
	N0_xi[8] = 0.25*eta*(eta+1)*(2*xi+1);
	return N0_xi;
}

//Calculates the derivative of shape functions with repsect to parameter eta
double* FEMSolver::N_eta(double xi, double eta){
	double *N0_eta = new double[9];
	N0_eta[0] = 0.25*xi*(xi-1)*(2*eta-1);
	N0_eta[1] = 0.5*(1-pow(xi,2))*(2*eta-1);
	N0_eta[2] = 0.25*xi*(xi+1)*(2*eta-1);
	N0_eta[3] = -xi*eta*(xi-1);
	N0_eta[4] = -2*eta*(1-pow(xi,2));
	N0_eta[5] = -xi*eta*(xi+1);
	N0_eta[6] = 0.25*xi*(xi-1)*(2*eta+1);
	N0_eta[7] = 0.5*(1-pow(xi,2))*(2*eta+1);
	N0_eta[8] = 0.25*xi*(xi+1)*(2*eta+1);
	return N0_eta;
}

//Calculates the derivative of global coordinate x with repsect to local parameter xi
double FEMSolver::x_xi(double xi, double eta, int* nodes, double* N0_xi){
	if(!N0_xi) N0_xi = N_xi(xi, eta);
	double x0_xi = 0;
	for(int i=0; i<9; i++){
		x0_xi = x0_xi + N0_xi[i]*(grid->pts[nodes[i]]).x; 
	}
	return x0_xi;
}

//Calculates the derivative of global coordinate x with repsect to local parameter eta
double FEMSolver::x_eta(double xi, double eta, int* nodes, double* N0_eta){
	if(!N0_eta) N0_eta = N_eta(xi, eta);
	double x0_eta = 0;
	for(int i=0; i<9; i++){
		x0_eta = x0_eta + N0_eta[i]*(grid->pts[nodes[i]]).x; 
	}
	return x0_eta;
}

//Calculates the derivative of global coordinate y with repsect to local parameter xi
double FEMSolver::y_xi(double xi, double eta, int* nodes, double* N0_xi){
	if(!N0_xi) N0_xi = N_xi(xi, eta);
	double y0_xi = 0;
	for(int i=0; i<9; i++){
		y0_xi = y0_xi + N0_xi[i]*(grid->pts[nodes[i]]).y; 
	}
	return y0_xi;
}

//Calculates the derivative of global coordinate y with repsect to local parameter eta
double FEMSolver::y_eta(double xi, double eta, int* nodes, double* N0_eta){
	if(!N0_eta) N0_eta = N_eta(xi, eta);
	double y0_eta = 0;
	for(int i=0; i<9; i++){
		y0_eta = y0_eta + N0_eta[i]*(grid->pts[nodes[i]]).y; 
	}
	return y0_eta;
}

//Calculates the derivative of shape functions with respect to global coordinate x
double* FEMSolver::N_x(double xi, double eta, int* nodes, double* N0_xi, double* N0_eta){
	if(!N0_xi) N0_xi = N_xi(xi, eta);
	if(!N0_eta) N0_eta = N_eta(xi, eta);
	double *N0_x = new double[9];
	double y0_xi = y_xi(xi, eta, nodes, N0_xi);
	double y0_eta = y_eta(xi, eta, nodes, N0_eta);
	double j = J(xi, eta, nodes, N0_xi, N0_eta);
	for(int i=0; i<9; i++){
		N0_x[i] = (N0_xi[i]*y0_eta - N0_eta[i]*y0_xi)/j;
	}
	return N0_x;
}

//Calculates the derivative of shape functions with respect to global coordinate y
double* FEMSolver::N_y(double xi, double eta, int* nodes, double* N0_xi, double* N0_eta){
	if(!N0_xi) N0_xi = N_xi(xi, eta);
	if(!N0_eta) N0_eta = N_eta(xi, eta);
	double *N0_y = new double[9];
	double x0_xi = x_xi(xi, eta, nodes, N0_xi);
	double x0_eta = x_eta(xi, eta, nodes, N0_eta);
	double j = J(xi, eta, nodes, N0_xi, N0_eta);
	for(int i=0; i<9; i++){
		N0_y[i] = (N0_eta[i]*x0_xi - N0_xi[i]*x0_eta)/j;
	}
	return N0_y;
}

//Evaluates the determinant of the Jacobian at a specified location
double FEMSolver::J(double xi, double eta, int* nodes, double* N0_xi, double* N0_eta){
	double j = x_xi(xi, eta, nodes, N0_xi)*y_eta(xi, eta, nodes, N0_eta) - x_eta(xi, eta, nodes, N0_eta)*y_xi(xi, eta, nodes, N0_xi);
	return j;
}

//Gives the absolute value of the Jacobian at a specified location
double FEMSolver::J_abs(double xi, double eta, int* nodes, double* N0_xi, double* N0_eta){
	double j = abs(J(xi, eta, nodes, N0_xi, N0_eta));
	return j;
}

//Evaluates the approximate integral of the product of the j-th and k-th shape functions(or their derivatives) using values calculated at the Gaussian quadrature points
double FEMSolver::GaussIntegralAB(double* J0, double** A, double** B, int j, int k){
	double I = (25*(J0[0]*A[0][j]*B[0][k] + J0[2]*A[2][j]*B[2][k] + J0[6]*A[6][j]*B[6][k] + J0[8]*A[8][j]*B[8][k]) + 64*J0[4]*A[4][j]*B[4][k] + 40*(J0[1]*A[1][j]*B[1][k] + J0[3]*A[3][j]*B[3][k] + J0[5]*A[5][j]*B[5][k] + J0[7]*A[7][j]*B[7][k]))/81;
	return I;
}

//Calculates the approximate integral of the i-th shape function along the xi-direction using Gaussian quadrature
double FEMSolver::ShapeFunctionIntegral_xi(double eta, int i){
	double r = sqrt(0.6);
	double *N0 = N(r, eta);
	double I = 5*N0[i]/9;
	N0 = N(0, eta);
	I = I + 8*N0[i]/9;
	N0 = N(-r, eta);
	I = I + 5*N0[i]/9;
	return I;
}

//Calculates the approximate integral of the i-th shape function along the eta-direction using Gaussian quadrature
double FEMSolver::ShapeFunctionIntegral_eta(double xi, int i){
	double r = sqrt(0.6);
	double *N0 = N(xi, r);
	double I = 5*N0[i]/9;
	N0 = N(xi, 0);
	I = I + 8*N0[i]/9;
	N0 = N(xi, -r);
	I = I + 5*N0[i]/9;
	return I;
}