#include "Universal.h"
#include "Matrix.h"
#include <cmath>
#include<iostream>
using namespace std;
Universal::Universal() {
	Xi.createTable(4, 1);
	Eta.createTable(4, 1);
	ip.createTable(2, 1);
	ip.matrixTable[0][0] = -1 / sqrt(3);
	ip.matrixTable[1][0] = 1 / sqrt(3);

	Xi.matrixTable[0][0] = Eta.matrixTable[0][0] = ip.matrixTable[0][0];
	Xi.matrixTable[1][0] = ip.matrixTable[1][0];
	Xi.matrixTable[2][0] = Eta.matrixTable[2][0] = ip.matrixTable[1][0];
	Xi.matrixTable[3][0] = ip.matrixTable[0][0];
	Eta.matrixTable[1][0]= ip.matrixTable[0][0];
	Eta.matrixTable[3][0]= ip.matrixTable[1][0];

}