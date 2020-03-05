#pragma once
#ifndef _UNIVERSAL_H
#include "Matrix.h"
struct Universal {
	Matrix ip; //integration point array
	Matrix Xi;
	Matrix Eta;

	Universal();
};
#define _UNIVERSAL_H
#endif
