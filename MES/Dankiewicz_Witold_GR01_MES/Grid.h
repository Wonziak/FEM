#pragma once
#ifndef _GRID_H
#include "Element.h"
#include "Node.h"
#include "Matrix.h"
#include "Matricies.h"
using namespace std;
struct Grid {

	Matricies matricies;
	Node * nodes;
	Element* elements;
	int numEl, numNodes;
	Grid();
	void setNodeBC(int nr);
	
	void setBC();
	void getCoordiantes();
	void createXpYp();
	void createJac1();
	void createJacInv();
	void createDetJac();
	void createdNdX();
	void createdNdY();
	void createPTab();
	void createCpTab();
	void createHLocal();
	void createCLocal();
	void calculateBorderCond();
	void createPLocal();
	void calculate();
	void calculateSolution();
	void showGrid();
	void GaussEl();
	int getBC(int nr);
	double getNodeY(int nr);
	double getNodeX(int nr);
	void setNodeX(int nr, double X);
	void setNodeY(int nr, double Y);
	int* getElementId(int nr);
	void createElementMatricies();
	void aggregateHCP();
	void matrixToFile();
	Matrix matrixH;
	Matrix matrixC;
	Matrix PVecG;
	Matrix t0Vec;
	Matrix calcMatrix;
};

	
#define _GRID_H
#endif