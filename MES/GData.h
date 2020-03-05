#pragma once
#ifndef _GData_H
#include <fstream>
#include<iostream>
using namespace std;

struct GlobalData {
	double H, W; //wysoko�� i szeroko�� w metrach
	int nH, nW; //ilo�� w�z��w w pionie/poziomie

	int numNodes, numEl; //ilo�� wszystkich w�z��w/element�w
	double dx, dy; //odleg�o�ci w�z��w
	double initTemp = 100;
	double dT = 50;
	double duration = 500;
	double alfa = 300;
	double ambTemp = 1200;
	double conductivity=25;
	double ro = 7800;
	double specHeat = 700;
	double cRo = ro * specHeat;
	GlobalData() {
		fstream file;
		file.open("mes.txt.txt");
		if (file.good()) {

			while (!file.eof()) {
				file >> H >> W >> nH >> nW;
			}
			numNodes = nH * nW;
			numEl = (nH - 1) * (nW - 1);
			dx = (W / (nW - 1.0));
			dy = (H / (nH - 1.0));
			cout << H <<" "<< W << " " << nH << " " << nW << " " <<dx<<" "<<dy<<" "<< numNodes << " " << numEl << endl;
		}
		else
			file.close();
	};
}gdata;
#define _GData_H
#endif



