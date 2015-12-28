


#ifndef FUNCTION_H
#define FUNCTION_H

#include <iostream>
#include <vector>
#include <fstream>

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <freeglut.h>
#include <math.h>

using namespace std;
using namespace Eigen;



struct location{
	float x;
	float y;
	float z;
	location(float x_, float y_, float z_) :x(x_), y(y_), z(z_) {}
};

// test input
// Point's (x,y,z)  Control lattice's = n,m,d 
inline float basisf(int i, float t)
{

	float b = 0;

	if (i == 1)
	{
		b = (((1 - t)*(1 - t)*(1 - t)) / 6);
	}

	if (i == 2)
	{
		b = (3 * t*t*t - 6 * t*t + 4) / 6;
	}

	if (i == 3)
	{
		b = (-3 * t*t*t + 3 * t*t + 3 * t + 1) / 6;
	}

	if (i == 4)
	{
		b = (t*t*t) / 6;
	}

	return b;
}

MatrixXf controlsetlo(float x, float y, float z, int n, int m, int d);
MatrixXf controlsetwkl(float x, float y, float z, int n, int m, int d);


inline float surface(MatrixXf xc, int pp, int oo, float u, float v)
{

	return 	basisf(1, u)* xc(pp, oo) + basisf(2, u)* xc(pp, oo + 1) + basisf(3, u)* xc(pp, oo + 2) + basisf(4, u)* xc(pp, oo + 3);

}


inline MatrixXf cal1(MatrixXf t, float u)
{
	t(0, 0) = u*u*u;
	t(0, 1) = u*u;
	t(0, 2) = u;
	t(0, 3) = 1;
	return t;
}

inline MatrixXf cal2(MatrixXf t, float u)
{
	t(0, 0) = u*u*u;
	t(1, 0) = u*u;
	t(2, 0) = u;
	t(3, 0) = 1;
	return t;
}



void refine(MatrixXf cps, MatrixXf & rcps);
void controlvalue(MatrixXf wkl, MatrixXf lot, MatrixXf& pz, int n, int m);
void valuez(int xsize, int ysize, MatrixXf xpos, MatrixXf ypos, MatrixXf & zpos, MatrixXf ct, int d1, float minz);
float diff(float x, float y, MatrixXf pz, int d1);

class controlset
{
public:

	float xpos;
	float ypos;
	float zpos;

	MatrixXf wkl, lot, pz;

	controlset();

};


struct Orientation
{
	float xpos;
	float ypos;
	float zpos;

};

float Fileload(char *filename, vector<Orientation>& points_1);

float Filesave(char *filename, vector<Orientation>& points_1);
float Filesave2(char *filename, MatrixXf contpz, MatrixXf contpx, MatrixXf contpy, float value);
float Filesave3(char *filename, vector<location> points_1, MatrixXf contpz, MatrixXf number, float value);



#endif
