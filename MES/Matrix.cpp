#include "Matrix.h"
#include <fstream>
#include <iostream>
using namespace std;
Matrix::Matrix(int x, int y) {
	rows = x;
	cols = y;
	matrixTable = new double*[x];
	for (int i = 0; i < x; i++) {
		matrixTable[i] = new double[y];
	}
	setZero();
}

Matrix::Matrix() {
}

void Matrix::createTable(int x, int y) {
	rows = x;
	cols = y;
	matrixTable = new double* [x];
	for (int i = 0; i < x; i++) {
		matrixTable[i] = new double[y];
	}
	setZero();
}
double Matrix::det() {
	if (rows != cols) {
		string ex = "cannot calculate determinant";
		throw ex;
	}
	if(rows==1)
	{
		return matrixTable[0][0];
	}
	if (rows == 2)
	{
		return (matrixTable[0][0]* matrixTable[1][1])-(matrixTable[0][1]* matrixTable[1][0]);
	}
}

double ** Matrix::matrixMultipication( Matrix m2) {
	if (this->cols != m2.rows) {
		string ex = "cannot multipicate matrices";
		throw ex;
	}
	Matrix result(this->rows, m2.cols);
	for (int i = 0; i < this->rows; i++)
	{
		for (int j = 0; j < m2.cols; j++)
		{
			for (int k = 0; k < this->cols; k++)
			{
				result.matrixTable[i][j] = result.matrixTable[i][j] + this->matrixTable[i][k] * m2.matrixTable[k][j];
			}
		}
	}
	return result.matrixTable;
}

double** Matrix::scalarMultipication(double scalar) {
	Matrix result;
	result.matrixTable = this->matrixTable;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			result.matrixTable[i][j]=result.matrixTable[i][j]*scalar;
		}
	}
	return result.matrixTable;
}

double** Matrix::scalarDivision(Matrix m, double scalar) {
	if (scalar == 0) { string ex = "division by 0"; throw ex; }
	Matrix result = m;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			result.matrixTable[i][j] = result.matrixTable[i][j] / scalar;
		}
	}
	return result.matrixTable;
}
double** Matrix::matrixSum(Matrix m) {
	if ((this->rows != m.rows) || (this->cols != m.cols)) {
		string ex = "cannot sum matrices";
		throw ex;
	}
	Matrix result(rows, cols);
	for (int i = 0; i < this->rows; i++) {
		for (int j = 0; j < this->cols; j++) {
			result.matrixTable[i][j] = this->matrixTable[i][j] + m.matrixTable[i][j];
		}
	}
	return result.matrixTable;
}
double** Matrix::getCol(int ColN) {
	Matrix result(rows,1);
	for (int i = 0; i < rows; i++) {
		result.matrixTable[i][0] = this->matrixTable[i][ColN];
	}
	return result.matrixTable;
}

double** Matrix::getRow(int RowN) {
	Matrix result(1,cols);
	for (int i = 0; i < cols; i++) {
		result.matrixTable[0][i] = matrixTable[RowN][i];
	}
	return result.matrixTable;
}

void Matrix::setZero() {
	for (int i = 0; i < rows;i++) {
		for (int j = 0; j < cols; j++) {
			this->matrixTable[i][j] = 0;
		}
	}
}
void Matrix::showMatrix() {
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			cout << matrixTable[i][j] << "         ";
		}
		cout << endl;
	}
	cout << endl << endl;
}


