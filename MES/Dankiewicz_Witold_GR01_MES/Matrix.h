#pragma once
#ifndef _Mateix_H
struct Matrix {
	int x, y;
	double rows, cols;
	double** matrixTable;
	Matrix(int x,int y);
	Matrix();
	double det();
	void createTable(int x, int y);
	double ** matrixMultipication(Matrix m2);
	double** matrixSum(Matrix m2);
	void setZero();
	void showMatrix();
	double** getRow(int RowN);
	double** getCol(int ColN);
	double** scalarMultipication(double scalar);
	double** scalarDivision(Matrix m, double scalar);
	
};
#define _Matrix_H
#endif