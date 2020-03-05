#include "Grid.h"
#include "GData.h"
#include "Node.h"
#include "Element.h"
#include "Matricies.h"
#include "Matrix.h"
#include "Universal.h"
#include <iostream>
#include <omp.h>

Grid::Grid() {
	omp_set_num_threads(8);
	numEl = gdata.numEl;
	numNodes = gdata.numNodes;
	nodes = new Node[numNodes];
	elements = new Element[numEl];
	t0Vec.createTable(numNodes,1);
	int nH = gdata.nH;
	int nW = gdata.nW;
	double dx = gdata.dx;
	double dy = gdata.dy;
	double x = 0.0, y = 0.0;
	int nodeNr = 0;
	for (int i = 0; i < numNodes; i++) {
		t0Vec.matrixTable[i][0] = gdata.initTemp;
	}
	for (int i = 0; i < numEl; i++) {
		if ((nodeNr+1)%gdata.nH==0) {
			cout << "node " << nodeNr << endl;
			nodeNr+=1; //next column
		}
		elements[i].elId = i;
		elements[i].nodesId[0] = nodeNr;
		elements[i].nodesId[1] = nodeNr + nH;
		elements[i].nodesId[2] = nodeNr + nH + 1;
		elements[i].nodesId[3] = nodeNr + 1;
		nodeNr++;
		}
	
	for (int i = 0; i < numNodes; i++) {
		if ((i + 1) % nH == 0 && i != 0) {
			nodes[i].id = i + 1;
			setNodeX(i, x);
			setNodeY(i, y);
			x += dx;
			y = 0;
		}
		else {
			setNodeX(i,x);
			setNodeY(i, y);
			y += dy;
		}
	}
	setBC();
	
	createElementMatricies();
	aggregateHCP();

}
void Grid::setBC() {
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numNodes; i++) {
			if (getNodeX(i) == 0 ) {
				setNodeBC(i);
			}
			if (gdata.W<= getNodeX(i)+0.000001) {
				setNodeBC(i);
			}
			if (getNodeY(i) == 0) {
				setNodeBC(i);
			}
			if (gdata.H <= getNodeY(i) + 0.000001) {
				setNodeBC(i);
			}
		}
}
void Grid::getCoordiantes() {
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numEl; i++) {
		elements[i].coordiantes.createTable(4, 2);

		for (int j = 0; j < 4; j++) {
			elements[i].coordiantes.matrixTable[j][0] = getNodeX(elements[i].nodesId[j]);
			elements[i].coordiantes.matrixTable[j][1] = getNodeY(elements[i].nodesId[j]);;
		}
	}
}

void Grid::createXpYp() {
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numEl; i++) {
		elements[i].XpYp.createTable(4, 2);
		elements[i].XpYp.matrixTable = matricies.matrixN.matrixMultipication(elements[i].coordiantes);
	}
}

void Grid::createJac1() {
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numEl; i++) {
		elements[i].Jac1.createTable(4, 4);

		for (int k = 0; k < 4; k++) {
			for (int j = 0; j < 4; j++) {
				elements[i].Jac1.matrixTable[0][k] += matricies.matrixXi.matrixTable[0][j] * elements[i].coordiantes.matrixTable[j][0];
				elements[i].Jac1.matrixTable[1][k] += matricies.matrixXi.matrixTable[1][j] * elements[i].coordiantes.matrixTable[j][1];
				elements[i].Jac1.matrixTable[2][k] += matricies.matrixEta.matrixTable[0][j] * elements[i].coordiantes.matrixTable[j][0];
				elements[i].Jac1.matrixTable[3][k] += matricies.matrixEta.matrixTable[1][j] * elements[i].coordiantes.matrixTable[j][1];
			}
		}
	}
}

void Grid::createJacInv() {
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numEl; i++) {
		elements[i].JacInv.createTable(4, 4);
		for (int j = 0; j < 4; j++) {
			elements[i].JacInv.matrixTable[0][j] = elements[i].Jac1.matrixTable[3][j] / elements[i].detJac.matrixTable[j][0];
			elements[i].JacInv.matrixTable[1][j] = elements[i].Jac1.matrixTable[1][j] / elements[i].detJac.matrixTable[j][0];
			elements[i].JacInv.matrixTable[2][j] = elements[i].Jac1.matrixTable[2][j] / elements[i].detJac.matrixTable[j][0];
			elements[i].JacInv.matrixTable[3][j] = elements[i].Jac1.matrixTable[0][j] / elements[i].detJac.matrixTable[j][0];
		}
	}
}

void Grid::createDetJac() {
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numEl; i++) {
		elements[i].detJac.createTable(4, 1);

		for (int j = 0; j < 4; j++) {
			elements[i].detJac.matrixTable[j][0] = elements[i].Jac1.matrixTable[0][j] * elements[i].Jac1.matrixTable[3][j] - elements[i].Jac1.matrixTable[1][j] * elements[i].Jac1.matrixTable[2][j];
		}
	}
}


void Grid::createdNdX() {
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numEl; i++) {
		elements[i].dX.createTable(4, 4);
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				elements[i].dX.matrixTable[j][k] = (elements[i].JacInv.matrixTable[0][j] * matricies.matrixXi.matrixTable[j][k]) + (elements[i].JacInv.matrixTable[1][j] * matricies.matrixEta.matrixTable[j][k]);

			}
		}
	}
}

void Grid::createdNdY() {
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numEl; i++) {
		elements[i].dY.createTable(4, 4);
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				elements[i].dY.matrixTable[j][k] = (elements[i].JacInv.matrixTable[2][j] * matricies.matrixXi.matrixTable[j][k]) + (elements[i].JacInv.matrixTable[3][j] * matricies.matrixEta.matrixTable[j][k]);
			}
		}
	}
}

void Grid::createPTab() {
	Matrix tmpX(1, 4);
	Matrix tmp2X(4, 1);
	Matrix tmpY(1, 4);
	Matrix tmp2Y(4, 1);
	Matrix tmp3(4, 4);
	for (int i = 0; i < numEl; i++) {
		for (int j = 0; j < 4; j++) {
			elements[i].pTab[j].createTable(4, 4);
				
				tmpX.matrixTable = elements[i].dX.getRow(j);
				tmpY.matrixTable = elements[i].dY.getRow(j);
				for (int k = 0; k < 4; k++) {
					tmp2X.matrixTable[k][0] = tmpX.matrixTable[0][k];
					tmp2Y.matrixTable[k][0] = tmpY.matrixTable[0][k];
				}
				elements[i].pTab[j].matrixTable = tmp2X.matrixMultipication(tmpX);

				tmp3.matrixTable= tmp2Y.matrixMultipication(tmpY);

				elements[i].pTab[j].matrixTable = tmp3.matrixSum(elements[i].pTab[j]);
				elements[i].pTab[j].matrixTable = elements[i].pTab[j].scalarMultipication(elements[i].detJac.matrixTable[j][0]);
				elements[i].pTab[j].matrixTable = elements[i].pTab[j].scalarMultipication(gdata.conductivity);
		}
	}

}

void Grid::createCpTab() {
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numEl; i++) {
		for (int j = 0; j < 4; j++) {
			elements[i].CpTab[j].createTable(4, 4);
			Matrix tmpN(1, 4);
			Matrix tmpN2(4, 1);
			tmpN.matrixTable= matricies.matrixN.getRow(j);
			for (int k = 0; k < 4; k++) {
				tmpN2.matrixTable[k][0] = tmpN.matrixTable[0][k];
			}
			double cRoJ = gdata.cRo * elements[i].detJac.matrixTable[j][0];
			elements[i].CpTab[j].matrixTable = tmpN2.matrixMultipication(tmpN);
			elements[i].CpTab[j].matrixTable = elements[i].CpTab[j].scalarMultipication(cRoJ);

		}
	}
}
void Grid::createCLocal() {
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numEl; i++) {
		elements[i].CL.createTable(4, 4);
		for (int j = 0; j < 4; j++) {
			elements[i].CL.matrixTable = elements[i].CpTab[j].matrixSum(elements[i].CL);
		}
	}
}
void Grid::createHLocal() {
#pragma omp parallel for schedule(static) 
	for (int i = 0; i < numEl; i++) {
		elements[i].HL.createTable(4, 4);
		for (int j = 0; j < 4; j++) {
			elements[i].HL.matrixTable = elements[i].pTab[j].matrixSum(elements[i].HL);
		}
	}
}

void Grid::calculateBorderCond() {
	double detJ1D;
	Matrix nTmp(2, 4);
	Matrix tmpB(4, 4);
	Matrix tmp(4, 1);
	Matrix tmp2(1, 4);
	Universal uni;
		for (int i = 0; i < numEl; i++) {
			elements[i].BC.createTable(4, 4);
			if (getBC(elements[i].nodesId[0]) == 1 && getBC(elements[i].nodesId[1]) == 1) {//pow 1
				detJ1D = (nodes[elements[i].nodesId[1]].x - nodes[elements[i].nodesId[0]].x) / 2;
				nTmp.setZero();
				tmpB.setZero();
				tmp.setZero();
				tmp2.setZero();

				nTmp.matrixTable[0][0] = 0.25 * (1 - uni.ip.matrixTable[0][0]) * 2;
				nTmp.matrixTable[0][1] = 0.25 * (1 + uni.ip.matrixTable[0][0]) * 2;
				nTmp.matrixTable[1][0] = 0.25 * (1 - uni.ip.matrixTable[1][0]) * 2;
				nTmp.matrixTable[1][1] = 0.25 * (1 + uni.ip.matrixTable[1][0]) * 2;

				tmp2.matrixTable = nTmp.getRow(0);
				for (int k = 0; k < 4; k++) {
					tmp.matrixTable[k][0] = tmp2.matrixTable[0][k];
				}
				tmpB.matrixTable = tmp.matrixMultipication(tmp2);
				tmpB.matrixTable = tmpB.scalarMultipication(gdata.alfa);
				tmpB.matrixTable = tmpB.scalarMultipication(detJ1D);
				elements[i].BC.matrixTable=elements[i].BC.matrixSum(tmpB);

				tmp2.matrixTable = nTmp.getRow(1);
				for (int k = 0; k < 4; k++) {
					tmp.matrixTable[k][0] = tmp2.matrixTable[0][k];
				}
				tmpB.matrixTable = tmp.matrixMultipication(tmp2);
				tmpB.matrixTable = tmpB.scalarMultipication(gdata.alfa);
				tmpB.matrixTable = tmpB.scalarMultipication(detJ1D);
				elements[i].BC.matrixTable = elements[i].BC.matrixSum(tmpB);
			}

			if (getBC(elements[i].nodesId[1]) == 1 && getBC(elements[i].nodesId[2]) == 1) {//pow 2
				detJ1D = (nodes[elements[i].nodesId[2]].y - nodes[elements[i].nodesId[1]].y) / 2;
				nTmp.setZero();
				tmpB.setZero();
				tmp.setZero();
				tmp2.setZero();

				nTmp.matrixTable[0][1] = 0.25 * 2*(1- uni.ip.matrixTable[0][0]);
				nTmp.matrixTable[0][2] = 0.25 * 2*(1 +uni.ip.matrixTable[0][0]);
				nTmp.matrixTable[1][1] = 0.25 * 2*(1 - uni.ip.matrixTable[1][0]);
				nTmp.matrixTable[1][2] = 0.25 * 2*(1 + uni.ip.matrixTable[1][0]);

				tmp2.matrixTable = nTmp.getRow(0);
				for (int k = 0; k < 4; k++) {
					tmp.matrixTable[k][0] = tmp2.matrixTable[0][k];
				}
				tmpB.matrixTable = tmp.matrixMultipication(tmp2);
				tmpB.matrixTable = tmpB.scalarMultipication(gdata.alfa);
				tmpB.matrixTable = tmpB.scalarMultipication(detJ1D);
				elements[i].BC.matrixTable = elements[i].BC.matrixSum(tmpB);

				tmp2.matrixTable = nTmp.getRow(1);
				for (int k = 0; k < 4; k++) {
					tmp.matrixTable[k][0] = tmp2.matrixTable[0][k];
				}
				tmpB.matrixTable = tmp.matrixMultipication(tmp2);
				tmpB.matrixTable = tmpB.scalarMultipication(gdata.alfa);
				tmpB.matrixTable = tmpB.scalarMultipication(detJ1D);
				elements[i].BC.matrixTable = elements[i].BC.matrixSum(tmpB);
			}

			if (getBC(elements[i].nodesId[2]) == 1 && getBC(elements[i].nodesId[3]) == 1) {//pow 3
				detJ1D = (nodes[elements[i].nodesId[2]].x - nodes[elements[i].nodesId[3]].x) / 2;
				nTmp.setZero();
				tmpB.setZero();
				tmp.setZero();
				tmp2.setZero();

				nTmp.matrixTable[0][2] = 0.25 * (1 + uni.ip.matrixTable[1][0]) * 2;
				nTmp.matrixTable[0][3] = 0.25 * (1 - uni.ip.matrixTable[1][0]) * 2;
				nTmp.matrixTable[1][3] = 0.25 * (1 + uni.ip.matrixTable[0][0]) * 2;
				nTmp.matrixTable[1][4] = 0.25 * (1 - uni.ip.matrixTable[0][0]) * 2;

				tmp2.matrixTable = nTmp.getRow(0);
				for (int k = 0; k < 4; k++) {
					tmp.matrixTable[k][0] = tmp2.matrixTable[0][k];
				}
				tmpB.matrixTable = tmp.matrixMultipication(tmp2);
				tmpB.matrixTable = tmpB.scalarMultipication(gdata.alfa);
				tmpB.matrixTable = tmpB.scalarMultipication(detJ1D);
				elements[i].BC.matrixTable = elements[i].BC.matrixSum(tmpB);

				tmp2.matrixTable = nTmp.getRow(1);
				for (int k = 0; k < 4; k++) {
					tmp.matrixTable[k][0] = tmp2.matrixTable[0][k];
				}
				tmpB.matrixTable = tmp.matrixMultipication(tmp2);
				tmpB.matrixTable = tmpB.scalarMultipication(gdata.alfa);
				tmpB.matrixTable = tmpB.scalarMultipication(detJ1D);
				elements[i].BC.matrixTable = elements[i].BC.matrixSum(tmpB);
			}

			if (getBC(elements[i].nodesId[3]) == true && getBC(elements[i].nodesId[0]) == true) {//pow 4
				detJ1D = (nodes[elements[i].nodesId[3]].y - nodes[elements[i].nodesId[0]].y) / 2;
				nTmp.setZero();
				tmpB.setZero();
				tmp.setZero();
				tmp2.setZero();

				nTmp.matrixTable[0][0] = 0.25 * 2*(1 - uni.ip.matrixTable[1][0]);
				nTmp.matrixTable[0][3] = 0.25 * 2*(1 + uni.ip.matrixTable[1][0]);
				nTmp.matrixTable[1][0] = 0.25 * 2*(1 - uni.ip.matrixTable[0][0]);
				nTmp.matrixTable[1][3] = 0.25 * 2*(1 + uni.ip.matrixTable[0][0]);

				tmp2.matrixTable = nTmp.getRow(0);
				for (int k = 0; k < 4; k++) {
					tmp.matrixTable[k][0] = tmp2.matrixTable[0][k];
				}
				tmpB.matrixTable = tmp.matrixMultipication(tmp2);
				tmpB.matrixTable = tmpB.scalarMultipication(gdata.alfa);
				tmpB.matrixTable = tmpB.scalarMultipication(detJ1D);
				elements[i].BC.matrixTable = elements[i].BC.matrixSum(tmpB);

				tmp2.matrixTable = nTmp.getRow(1);
				for (int k = 0; k < 4; k++) {
					tmp.matrixTable[k][0] = tmp2.matrixTable[0][k];
				}
				tmpB.matrixTable = tmp.matrixMultipication(tmp2);
				tmpB.matrixTable = tmpB.scalarMultipication(gdata.alfa);
				tmpB.matrixTable = tmpB.scalarMultipication(detJ1D);
				elements[i].BC.matrixTable = elements[i].BC.matrixSum(tmpB);
			}
			elements[i].HL.matrixTable = elements[i].HL.matrixSum(elements[i].BC);
		}
}

void Grid::createElementMatricies() {
	getCoordiantes();
	createXpYp();
	createJac1();
	createDetJac();
	createJacInv();
	createdNdX();
	createdNdY();
	createPTab();
	createHLocal();
	createCpTab();
	createCLocal();
	createPLocal();
	calculateBorderCond();
}

void Grid::aggregateHCP() {
	matrixH.createTable(numNodes, numNodes);
	matrixC.createTable(numNodes, numNodes);
	PVecG.createTable(numNodes, 1);
	int tabID[4] = { 0,0,0,0 };
	for (int i = 0; i < numEl; i++) {
		for (int j = 0; j < 4; j++) {
			tabID[j] = elements[i].nodesId[j];
		}
		for (int j = 0; j < 4; j++) {
			PVecG.matrixTable[tabID[j]][0] += elements[i].PVecL.matrixTable[j][0];
			for (int k = 0; k < 4; k++) {
				matrixH.matrixTable[tabID[j]][tabID[k]] += elements[i].HL.matrixTable[j][k];
				matrixC.matrixTable[tabID[j]][tabID[k]] += elements[i].CL.matrixTable[j][k];
			}
		}
	}
}

void Grid::createPLocal(){
	double detJ1D=0;
	double alfaTemp = gdata.alfa * gdata.ambTemp;
	Universal uni;

	for (int i = 0; i < numEl; i++) {
		Matrix p(4, 1);
		elements[i].PVecL.createTable(4, 1);
		if (getBC(elements[i].nodesId[0]) == 1 && getBC(elements[i].nodesId[1]) == 1) {//pow 1
			detJ1D = (gdata.dx) / 2;

			 p.matrixTable[0][0] += 0.5 * (1 - uni.ip.matrixTable[0][0]) * detJ1D * alfaTemp;//-
			 p.matrixTable[1][0] += 0.5 * (1 - uni.ip.matrixTable[1][0]) * detJ1D * alfaTemp;

			 p.matrixTable[0][0] += 0.5 * (1 + uni.ip.matrixTable[0][0]) * detJ1D * alfaTemp;
			 p.matrixTable[1][0] += 0.5 * (1 + uni.ip.matrixTable[1][0]) * detJ1D * alfaTemp;

		}

		if (getBC(elements[i].nodesId[1]) == 1 && getBC(elements[i].nodesId[2]) == 1) {//pow 2
			detJ1D = (gdata.dy) / 2;
			p.matrixTable[1][0] += 0.5 * (1 - uni.ip.matrixTable[0][0]) * detJ1D * alfaTemp;
			p.matrixTable[2][0] += 0.5 * (1 - uni.ip.matrixTable[1][0]) * detJ1D * alfaTemp;

			p.matrixTable[1][0] += 0.5 * (1 + uni.ip.matrixTable[0][0]) * detJ1D * alfaTemp;
			p.matrixTable[2][0] += 0.5 * (1 + uni.ip.matrixTable[1][0]) * detJ1D * alfaTemp;


		}

		if (getBC(elements[i].nodesId[2]) == 1 && getBC(elements[i].nodesId[3]) == 1) {//pow 3
			detJ1D = (gdata.dx) / 2;

			p.matrixTable[2][0] += 0.5 * (1 + uni.ip.matrixTable[1][0]) * detJ1D * alfaTemp;
			p.matrixTable[3][0] += 0.5 * (1 + uni.ip.matrixTable[0][0]) * detJ1D * alfaTemp;
			
			p.matrixTable[2][0] += 0.5 * (1 - uni.ip.matrixTable[1][0]) * detJ1D * alfaTemp;
			p.matrixTable[3][0] += 0.5 * (1 - uni.ip.matrixTable[0][0]) * detJ1D * alfaTemp;
		


			
		}

		if (getBC(elements[i].nodesId[3]) == true && getBC(elements[i].nodesId[0]) == true) {//pow 4
			detJ1D = (gdata.dy) / 2;
			p.matrixTable[3][0] += 0.5 * (1 + uni.ip.matrixTable[1][0])*detJ1D* alfaTemp;
			p.matrixTable[0][0] += 0.5 * (1 + uni.ip.matrixTable[0][0]) * detJ1D * alfaTemp;

			p.matrixTable[3][0] += 0.5 * (1 - uni.ip.matrixTable[1][0]) * detJ1D * alfaTemp;
			p.matrixTable[0][0] += 0.5 * (1 - uni.ip.matrixTable[0][0]) * detJ1D * alfaTemp;

		}
		for (int j = 0; j < 4; j++)
		{
			elements[i].PVecL.matrixTable[j][0] = p.matrixTable[j][0];
		}
	}
}


void Grid::calculate() {
	for (int i = 0; i < numNodes; i++) {
		for (int j = 0; j < numNodes; j++) {
			matrixH.matrixTable[i][j] += matrixC.matrixTable[i][j] / gdata.dT;
		}

	}

	Matrix tmp2(numNodes, numNodes);

	for (int i = 0; i < numNodes; i++) {
		for (int j = 0; j < numNodes; j++) {
			tmp2.matrixTable[i][j] = matrixC.matrixTable[i][j] / gdata.dT;
		}
	}
	Matrix tmp(numNodes, 1);
	tmp.matrixTable = tmp2.matrixMultipication(t0Vec);

	for (int i = 0; i < numNodes; i++) {
			PVecG.matrixTable[i][0] += tmp.matrixTable[i][0];
	}
	cout << "wektor P\n";
	PVecG.showMatrix();
	
}

void Grid::calculateSolution() {
	calcMatrix.createTable(numNodes, numNodes + 1);
	for (int k = 0; k <= gdata.duration; k += gdata.dT) {
		calculate();

		for (int i = 0; i < numNodes; i++) {
			for (int j = 0; j < numNodes; j++) {
				calcMatrix.matrixTable = matrixH.matrixTable;	
			}
			calcMatrix.matrixTable[i][numNodes] = PVecG.matrixTable[i][0];
		}
	
		double max = 0;
		double min = 2000;
		GaussEl();
		for (int i = 0; i < numNodes; i++) {
			if (t0Vec.matrixTable[i][0] > max) {
				max = t0Vec.matrixTable[i][0];
			}
			if (t0Vec.matrixTable[i][0] < min) {
				min = t0Vec.matrixTable[i][0];
			}
		}
		if (k%int(gdata.dT)==0&&k!=0) {
			cout << "Po uplywie "<<k<<" sekund, temperatura minimalna: " << min << "\tTemperatura maksymalna : " <<max<<endl;
		}
	}
}
void Grid::GaussEl() {
	
	for (int i = 0; i < numNodes; i++) {
		// Search for maximum in this column
		double maxEl = abs(calcMatrix.matrixTable[i][i]);
		int maxRow = i;
		for (int k = i + 1; k < numNodes; k++) {
			if (abs(calcMatrix.matrixTable[k][i]) > maxEl) {
				maxEl = abs(calcMatrix.matrixTable[k][i]);
				maxRow = k;
			}
		}
		// Swap maximum row with current row (column by column)
		for (int k = i; k < numNodes + 1; k++) {
			double tmp = calcMatrix.matrixTable[maxRow][k];
			calcMatrix.matrixTable[maxRow][k] = calcMatrix.matrixTable[i][k];
			calcMatrix.matrixTable[i][k] = tmp;
		}

		// Make all rows below this one 0 in current column
		for (int k = i + 1; k < numNodes; k++) {
			double c = -calcMatrix.matrixTable[k][i] / calcMatrix.matrixTable[i][i];
			for (int j = i; j < numNodes + 1; j++) {
				if (i == j) {
					calcMatrix.matrixTable[k][j] = 0;
				}
				else {
					calcMatrix.matrixTable[k][j] += c * calcMatrix.matrixTable[i][j];
				}
			}
		}
	}

	// Solve equation Ax=b for an upper triangular matrix A
	for (int i = numNodes - 1; i >= 0; i--) {
		t0Vec.matrixTable[i][0] = calcMatrix.matrixTable[i][numNodes] / calcMatrix.matrixTable[i][i];
		for (int k = i - 1; k >= 0; k--) {
			calcMatrix.matrixTable[k][numNodes] -= calcMatrix.matrixTable[k][i] * t0Vec.matrixTable[i][0];
		}
	}


}
void Grid::showGrid() {
	for (int i = 0; i < numEl; i++) {
		cout << "Element " << i  << ":  nodes [" << elements[i].nodesId[0] << ", " << elements[i].nodesId[1] << ", "
			<< elements[i].nodesId[2] << ", " << elements[i].nodesId[3] << "] " << endl;
	}
}

int Grid::getBC(int nr) {
	return nodes[nr].BC;
}
void Grid::setNodeBC(int nr) {

	nodes[nr].BC = 1;
}
double Grid::getNodeY( int nr) {
	return nodes[nr].y; // zwraca wspó³rzêdn¹ danego wêz³a
}

void Grid::setNodeX(int nr,double X) {
	nodes[nr].x = X;
}
void Grid::setNodeY(int nr,double Y) {
	nodes[nr].y = Y;
}
double Grid::getNodeX(int nr) {

	return nodes[nr].x; // zwraca wspó³rzêdn¹ danego wêz³a
}
int* Grid::getElementId(int nr) {
	return elements[nr].nodesId; //zwraca ID - wêz³y dla elementu o danym numerze
}


void Grid::matrixToFile() {
	ofstream file("matrixH_MatrixC.txt");
	if (file.is_open()) {
		file << "matrix H\n";
		for (int i = 0; i < numNodes; i++) {
			for (int j = 0; j < numNodes; j++)
			{
				file << matrixH.matrixTable[i][j] << " ";
			}
			file << "\n";
		}
		file << "matrix C\n";
		for (int i = 0; i < numNodes; i++) {
			for (int j = 0; j < numNodes; j++)
			{
				file << matrixC.matrixTable[i][j] << " ";
			}
			file << "\n";
		}
		file << "vector P\n";
		for (int i = 0; i < numNodes; i++) {
				file << PVecG.matrixTable[i][0] << " ";
		}
		file << "\n";
		file << "matrix calculated\n";
		for (int i = 0; i < numNodes; i++) {
			for (int j = 0; j < numNodes+1; j++)
			{
				file << calcMatrix.matrixTable[i][j] << " ";
			}
			file << "\n";
		}
	}
}
