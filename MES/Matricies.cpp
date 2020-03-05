#include "Matrix.h"
#include "Matricies.h"
#include <iostream>
#include "Universal.h"
using namespace std;

Matricies::Matricies() {

	matrixN.createTable(4, 4);
	matrixEta.createTable(4, 4);
	matrixXi.createTable(4, 4);
	Universal uni;
	for (int i = 0; i < 4; i++) {

		matrixN.matrixTable[i][0] = (0.25 * (1 - uni.Xi.matrixTable[i][0]) * (1 - uni.Eta.matrixTable[i][0])); 
		matrixN.matrixTable[i][1] = (0.25 * (1 + uni.Xi.matrixTable[i][0]) * (1 - uni.Eta.matrixTable[i][0]));
		matrixN.matrixTable[i][2] = (0.25 * (1 + uni.Xi.matrixTable[i][0]) * (1 + uni.Eta.matrixTable[i][0]));
		matrixN.matrixTable[i][3] = (0.25 * (1 - uni.Xi.matrixTable[i][0]) * (1 + uni.Eta.matrixTable[i][0]));

		matrixXi.matrixTable[i][0] = (-0.25 * (1 - uni.Eta.matrixTable[i][0]));
		matrixXi.matrixTable[i][1] = (0.25 * (1 - uni.Eta.matrixTable[i][0]));
		matrixXi.matrixTable[i][2] = (0.25 * (1 + uni.Eta.matrixTable[i][0]));
		matrixXi.matrixTable[i][3] = (-0.25 * (1 + uni.Eta.matrixTable[i][0]));

		matrixEta.matrixTable[i][0] = (-0.25 * (1 - uni.Xi.matrixTable[i][0]));
		matrixEta.matrixTable[i][1] = (-0.25 * (1 + uni.Xi.matrixTable[i][0]));
		matrixEta.matrixTable[i][2] = (0.25 * (1 + uni.Xi.matrixTable[i][0]));
		matrixEta.matrixTable[i][3] = (0.25 * (1 - uni.Xi.matrixTable[i][0]));
	}

}

void Matricies::showMatricies() {
	cout << endl;
	cout << "matrix N\n      N1            N2              N3              N4\n";

	matrixN.showMatrix();
	cout << endl;

	cout << "matrix dN/dXi\n        N1              N2             N3               N4\n";

	matrixXi.showMatrix();
	cout << endl;

	cout << "matrix dN/dEta\n       N1             N2            N3               N4\n";

	matrixEta.showMatrix();
}