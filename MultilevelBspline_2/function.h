


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
float basisf(int i, float t);
MatrixXf controlsetlo(float x, float y, float z, int n, int m, int d);
MatrixXf controlsetwkl(float x, float y, float z, int n, int m, int d);

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
