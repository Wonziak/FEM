#include <iostream>
#include <string>
#include <math.h>
#include "Element.h"
#include "Grid.h"
#include "Node.h"
#include "Universal.h"
#include "Matrix.h"
#include "Matricies.h"
#include <chrono>
using namespace std;



int main() {
	Grid siatka;
	cout << endl;
	siatka.showGrid();
	cout << endl;
	Matricies matricies;
	matricies.showMatricies();
	siatka.calculateSolution();
	siatka.matrixToFile();
	system("pause");
	return 0;
}