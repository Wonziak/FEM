#pragma once
#pragma once
#ifndef _MATRICIES_H
#include "Matrix.h"
using namespace std;
struct Matricies {
	Matricies();
	void showMatricies();
	double Xi = 0, Eta = 0;
	Matrix matrixN; //matrix n1 n2 n3 n4
	Matrix matrixXi; //matrix dN/dXi
	Matrix matrixEta; //matrix dN/dEta
};


#define _MATRICIES_H
#endif