#pragma once
#ifndef _Element_H
#include "Matrix.h"
	struct Element {
		int elId;
		int nodesId[4];
		Matrix coordiantes;
		Matrix XpYp;
		Matrix dX;
		Matrix dY;
		Matrix Jac1;
		Matrix JacInv;
		Matrix detJac;
		Matrix HL;
		Matrix CL;
		Matrix BC;
		Matrix p1, p2, p3, p4;
		Matrix Cp1, Cp2, Cp3, Cp4;
		Matrix CpTab[4] = { Cp1,Cp2,Cp3,Cp4 };
		Matrix pTab[4] = { p1,p2,p3,p4 };
		Matrix PVecL;
	};
#define _Element_H
#endif