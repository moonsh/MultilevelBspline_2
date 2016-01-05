

#include <iostream>
#include <istream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <freeglut.h>
#include "function.h"
#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>



using namespace std;

typedef float vec3_t[3];
#define X_AXIS                  0
#define Y_AXIS                  1
#define Z_AXIS                  2

enum {
	UP = 1,
	DOWN,
};
int oldX = -4;
int oldY = -4;
int mState = UP;
vec3_t gRot = { 0, 0, 0 };


// Domain ohm -> { x, y | 0 < x < m / d, 0 < y < n / d }

int n;  // y range
int m;  // x range 

int d1 = 10;  // x,y resolution   control set #1
float v = 0.005;  // parameter interval 
float RMSE;
int latticen = 2;

int sn;  // Scatter points number;
int mAction = 0;
bool gouraud = false;

float n1;  // y range
float m1;  // x range 

int d2, i, j, k, l, oo, pp, im, in, check;
float o, p, s, t, d, c, rmse, maxx, maxy, miny, minx, minz;
float xcp, ycp;
float s1, t1, u1;

vector <float> xs, ys;

MatrixXf x5, y5, z5, pz1, pz2, z5_loca, pz1_loca;
MatrixXf xc, yc, bsm, tpa,contx,conty,contz,uc,vc,xcp2,ycp2,zcp2;

float zp[32000];

float XUP[3] = { 1, 0, 0 }, XUN[3] = { -1, 0, 0 },
YUP[3] = { 0, 1, 0 }, YUN[3] = { 0, -1, 0 },
ZUP[3] = { 0, 0, 1 }, ZUN[3] = { 0, 0, -1 },
ORG[3] = { 0, 0, 0 };

GLsizei lastWidth, lastHeight;


void display();
void SetupRC(); void Triad(void);
void clamp(vec3_t v);
void glutMotion(int x, int y);
void glutMouse(int button, int state, int x, int y);
void DoMenu(int value);

vector<Orientation> points_1;
vector<location> l1;

int main(int argc, char* argv[])
{
	bsm.resize(4, 4);
	contx.resize(4, 4);
	conty.resize(4, 4);
	contz.resize(4, 4);
	xcp2.resize(1, 1);
	ycp2.resize(1, 1);
	zcp2.resize(1, 1);
	uc.resize(1, 4);
	vc.resize(4, 1);

	uc.setZero();
	vc.setZero();
	bsm.setZero();
	bsm <<  -1, 3, -3, 1, 
			3, -6, 3, 0,
			-3, 0, 3, 0,
			1, 4, 1, 0 ;
	bsm = bsm / 6.0;

	

	cout << "Max RMSE  :";
	cin >> RMSE;

	clock_t start, finish;
	double duration;
	start = clock();


	cout << " data loading " << endl;


	Fileload("Result1.obj", points_1);
	cout << " data loaded " << endl;

	sn = points_1.size();

	minx = 100;
	miny = 100;
	minz = 0;

	for (i = 0; i < sn; i = i + 1)
	{
		if (minx > points_1.at(i).xpos)
		{
			minx = points_1.at(i).xpos;
		}
		if (miny > points_1.at(i).ypos)
		{
			miny = points_1.at(i).ypos;
		}
	}

	for (i = 0; i < sn; i = i + 1)
	{
		minz = points_1.at(i).zpos + minz;
	}

	minz = minz / sn;
	minz = 0;
	im = minx * d1 - 4;
	in = miny * d1 - 4;

	maxx = 0.002;
	maxy = 0.002;

	for (i = 0; i < sn; i = i + 1)
	{
		if (maxx < points_1.at(i).xpos)
		{
			maxx = points_1.at(i).xpos;
		}
		if (maxy < points_1.at(i).ypos)
		{
			maxy = points_1.at(i).ypos;
		}
	}


	m = maxx * d1 + 1;
	n = maxy * d1 + 1;

	n1 = (float)n / (float)d1;  // y range
	m1 = (float)m / (float)d1;  // x range 

	// Make contol points set to get contol values

	MatrixXf lo1(n + 3, m + 3);
	MatrixXf lot(n + 3, m + 3);
	lot.setZero();

	MatrixXf wkl1(n + 3, m + 3);
	MatrixXf wklt(n + 3, m + 3);
	wklt.setZero();


	for (i = 0; i < sn; i = i + 1)
	{
		points_1.at(i).zpos = points_1.at(i).zpos - minz;
	}


	for (i = 0; i < sn; i = i + 1)
	{
		wkl1 = controlsetwkl(points_1.at(i).xpos, points_1.at(i).ypos, points_1.at(i).zpos, n, m, d1);
		wklt = wkl1 + wklt;

		lo1 = controlsetlo(points_1.at(i).xpos, points_1.at(i).ypos, points_1.at(i).zpos, n, m, d1);
		lot = lo1 + lot;
	}


	// Assign control points values
	pz1.resize(n + 3, m + 3);
	controlvalue(wklt, lot, pz1, n, m);

	// ----------------------------------------------------------------------------------------------- 
	//  control lattice  반복 구문

	while (1)
	{

		MatrixXf pz1_1(2 * n + 3, 2 * m + 3);
		refine(pz1, pz1_1);

		for (i = 0; i < sn; i = i + 1)
		{
			zp[i] = points_1.at(i).zpos - diff(points_1.at(i).xpos, points_1.at(i).ypos, pz1, d1);
		}

		lo1.resize(2 * n + 3, 2 * m + 3);
		lot.resize(2 * n + 3, 2 * m + 3);
		lot.setZero();

		wkl1.resize(2 * n + 3, 2 * m + 3);
		wklt.resize(2 * n + 3, 2 * m + 3);
		wklt.setZero();

		for (i = 0; i < sn; i = i + 1)
		{
			wkl1 = controlsetwkl(points_1.at(i).xpos, points_1.at(i).ypos, zp[i], 2 * n, 2 * m, 2 * d1);
			wklt = wkl1 + wklt;

			lo1 = controlsetlo(points_1.at(i).xpos, points_1.at(i).ypos, zp[i], 2 * n, 2 * m, 2 * d1);
			lot = lo1 + lot;
		}

		pz1.resize(2 * n + 3, 2 * m + 3);
		pz1_loca.resize(2 * n + 3, 2 * m + 3);
		pz2.resize(2 * n + 3, 2 * m + 3);
		xc.resize(2 * n + 3, 2 * m + 3);
		yc.resize(2 * n + 3, 2 * m + 3);

		pz1.setZero();
		pz1_loca.setZero();
		xc.setZero();
		yc.setZero();
		controlvalue(wklt, lot, pz1, 2 * n, 2 * m);

		pz1 = pz1 + pz1_1;
		pz2 = pz1;

		rmse = 0;
		for (i = 0; i < sn; i = i + 1)
		{
			zp[i] = zp[i] * zp[i];
			rmse = zp[i] + rmse;
		}

		rmse = rmse / sn;
		rmse = sqrt(rmse);
		cout << "RMSE  = " << rmse * 70 << "mm " << endl;
		cout << "Resolution  = " << d1 * 2 << endl;

		if ((rmse <= RMSE) && (rmse >= (-RMSE)))
		{
			break;
		}

		d1 = d1 * 2;
		n = n * 2;
		m = m * 2;

		d2 = d1;
		latticen = latticen + 1;

	}

	cout << "Lattice number  = " << latticen << endl;


	// ---------------------------------------------------------Drawing graph-------------------------------------- 


	for (oo = 0; oo < xc.cols(); oo = oo + 1)
	{
		for (pp = 0; pp < xc.rows(); pp = pp + 1)
		{
			xc(pp, oo) = (v*10 * (oo)-v*10);
			yc(pp, oo) = (v*10 * (pp)-v*10);
		}
	}




	for (o = 0; o < (m1 - 0.0001); o = o + v)
	{
		xs.push_back(o);
	}
	for (p = 0; p < (n1 - 0.0001); p = p + v)
	{
		ys.push_back(p);
	}

	x5.resize(ys.size(), xs.size());
	y5.resize(ys.size(), xs.size());

	for (oo = 0; oo < xs.size(); oo = oo + 1)
	{
		for (pp = 0; pp < ys.size(); pp = pp + 1)
		{
			x5(pp, oo) = (v * (oo));
			y5(pp, oo) = (v * (pp));
		}
	}

	z5.resize(ys.size(), xs.size());
	z5.setZero();
	z5_loca.resize(ys.size(), xs.size());
	z5_loca.setZero();
	valuez(xs.size(), ys.size(), x5, y5, z5, pz1, d1 * 2, minz);     // z5 값이 바로 계산된 함수값.


	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;

	cout << " Construction time : " << duration << endl;

	cout << " OBJ File saving " << endl;

	int numbering;
	numbering = 1;
	for (pp = 0; pp < ys.size(); pp = pp + 1)
	{
		for (oo = 0; oo < xs.size(); oo = oo + 1)
		{
			if (z5(pp, oo) > 0.1)
			{
				l1.push_back({ x5(pp, oo), y5(pp, oo), z5(pp, oo) });
				z5_loca(pp, oo) = numbering;
				numbering = numbering + 1;
			}
		}
	}


	for (pp = 0; pp < xc.rows(); pp = pp + 1)
	{
		for (oo = 0; oo < xc.cols(); oo = oo + 1)
		{
			if (pz2(pp, oo) > 0.1)
			{
				pz1_loca(pp, oo) = 1;
			}


		}
	}




	check = 1;
	for (pp = 0; pp < xc.rows(); pp = pp + 1)
	{
		for (oo = 0; oo < xc.cols(); oo = oo + 1)
		{
			check = 1;
			if (pz1_loca(pp, oo) == 0)
			{
				
				
					for (i = oo; i < xc.cols(); i = i + 1)
					{
						if ((pz1_loca(pp, i) == 1) && (check == 1))
						{
							xcp = xc(pp,i);
							ycp = pz2(pp,i);
							check = check + 1;
						}
					}

					if (check==2)
					{ 
				xc(pp, oo) = xcp;
				pz2(pp, oo) = ycp;
					}
				
			}

		}
	}



	for (pp = xc.rows()-1; pp >= 0; pp = pp - 1)
	{
		for (oo = xc.cols()-1; oo >= 0; oo = oo - 1)
		{
			check = 1;
			if (pz1_loca(pp, oo) == 0)
			{


				for (i = oo; i >= 0; i = i - 1)
				{
					if ((pz1_loca(pp, i) == 1) && (check == 1))
					{
						xcp = xc(pp, i);
						ycp = pz2(pp, i);
						check = check + 1;
					}
				}

				if (check == 2)
				{
					xc(pp, oo) = xcp;
					pz2(pp, oo) = ycp;
				}

			}

		}
	}




	for (pp = xc.rows() - 1; pp >= 0; pp = pp - 1)
	{
		check = 1;

		for (oo = 0; oo < xc.cols(); oo = oo + 1)
		{
			if (pz1_loca(pp, oo) == 0)
			{

				for (i = oo; i < xc.cols(); i = i + 1)
				{
					if ((pz1_loca(pp, i) == 1) && (check == 1))
					{

						check = check + 1;
					}
				}
			}


		}

		if (check == 1)
		{
			for (oo = 0; oo < xc.cols(); oo = oo + 1)
			{
				xc(pp, oo) = xc(pp + 1, oo);
				yc(pp, oo) = yc(pp + 1, oo);
				pz2(pp, oo) = pz2(pp + 1, oo);
			}
		}


	}

	cout << cal1(uc, 0.1) << endl;
	cout << cal2(vc, 0.1) << endl;
	Filesave2("Data3.obj", pz2, xc, yc, 1);
	Filesave3("Data2.obj", l1, z5, z5_loca, 0.1);

	cout << " m1 " << pz2.cols() << endl;
	cout << " n1 " << pz2.rows() << endl;



	glutInit(&argc, argv);
	glutInitWindowSize(1200, 900);
	glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH);
	glutCreateWindow("Multilevel B spline");
	glutDisplayFunc(display);
	glutMotionFunc(glutMotion);
	glutCreateMenu(DoMenu);
	glutAddMenuEntry("Gouraud Shading on/off", 1);
	glutAttachMenu(GLUT_RIGHT_BUTTON);

	glutMouseFunc(glutMouse);
	glutMainLoop();

	return 0;


}



void display()
{

	SetupRC();
	GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat mat_shininess[] = { 20.0 };
	GLfloat light_position[] = { -0.6, 1.0, 1.0, 0.0 };
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glShadeModel(GL_SMOOTH);

	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);

	glRotatef(gRot[0], 1.0, 0.0, 0.0);
	glRotatef(gRot[1], 0.0, 1.0, 0.0);



	Triad();

	//glBegin(GL_POINTS);



	glColor3f(1.0, 1.0, 1.0);
	glPointSize(1.0);
	glEnable(GL_NORMALIZE);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);
	//cout << bsm.cols() << endl;
	for (oo = 0; oo < xc.cols()-3 ; ++oo)
	{
		for (pp = 0; pp < xc.rows()-3 ; ++pp) 
		{ 

			for (i = 0; i < 4; i = i + 1)

			{
				for (j = 0; j < 4; j = j + 1)
				{

					contx(i, j) = xc(pp + i, oo + j);
					conty(i, j) = yc(pp + i, oo + j);
					contz(i, j) = pz2(pp + i, oo + j);
				}

			}

			for (t = 0; t < 1; t = t + 0.25)
			{

				for (u1 = 0; u1 < 1; u1 = u1 + 0.25)
				{ 
					glBegin(GL_POLYGON);

					glNormal3f(
						((cal1(uc, u1 + 0.5)*bsm*conty *bsm.transpose() * cal2(vc, t + 0.5))(0, 0) - (cal1(uc, u1 + 0.5)*bsm*conty *bsm.transpose() * cal2(vc, t + 0.5))(0, 0)) * ((cal1(uc, u1)*bsm*contz *bsm.transpose() * cal2(vc, t))(0, 0) - (cal1(uc, u1)*bsm*contz *bsm.transpose() * cal2(vc, t + 0.5))(0, 0))
						- ((cal1(uc, u1 + 0.5)*bsm*contz *bsm.transpose() * cal2(vc, t + 0.5))(0, 0) - (cal1(uc, u1)*bsm*contz *bsm.transpose() * cal2(vc, t + 0.5))(0, 0)) * ((cal1(uc, u1)*bsm*conty *bsm.transpose() * cal2(vc, t))(0, 0) - (cal1(uc, u1)*bsm*conty *bsm.transpose() * cal2(vc, t + 0.5))(0, 0)),
						
						((cal1(uc, u1 + 0.5)*bsm*contz *bsm.transpose() * cal2(vc, t + 0.5))(0, 0) - (cal1(uc, u1)*bsm*contz *bsm.transpose() * cal2(vc, t + 0.5))(0, 0)) * ((cal1(uc, u1)*bsm*contx *bsm.transpose() * cal2(vc, t))(0, 0) - (cal1(uc, u1)*bsm*contx *bsm.transpose() * cal2(vc, t + 0.5))(0, 0))
						- ((cal1(uc, u1 + 0.5)*bsm*contx *bsm.transpose() * cal2(vc, t + 0.5))(0, 0) - (cal1(uc, u1)*bsm*contx *bsm.transpose() * cal2(vc, t + 0.5))(0, 0))* ((cal1(uc, u1)*bsm*contz *bsm.transpose() * cal2(vc, t))(0, 0) - (cal1(uc, u1)*bsm*contz *bsm.transpose() * cal2(vc, t + 0.5))(0, 0)),

						((cal1(uc, u1 + 0.5)*bsm*contx *bsm.transpose() * cal2(vc, t + 0.5))(0, 0) - (cal1(uc, u1)*bsm*contx *bsm.transpose() * cal2(vc, t + 0.5))(0, 0))* ((cal1(uc, u1)*bsm*conty *bsm.transpose() * cal2(vc, t))(0, 0) - (cal1(uc, u1)*bsm*conty *bsm.transpose() * cal2(vc, t + 0.5))(0,0))
						- ((cal1(uc, u1 + 0.5)*bsm*conty *bsm.transpose() * cal2(vc, t + 0.5))(0, 0) - (cal1(uc, u1)*bsm*conty *bsm.transpose() * cal2(vc, t + 0.5))(0, 0)) * ((cal1(uc, u1)*bsm*contx *bsm.transpose() * cal2(vc, t))(0, 0) - (cal1(uc, u1)*bsm*contx *bsm.transpose() * cal2(vc, t + 0.5))(0, 0))


						);


					glVertex3f( (cal1(uc, u1)*bsm*contx *bsm.transpose() * cal2(vc, t))(0, 0),
								(cal1(uc, u1)*bsm*conty *bsm.transpose() * cal2(vc, t))(0, 0),
								(cal1(uc, u1)*bsm*contz *bsm.transpose() * cal2(vc, t))(0, 0)
							  );

					glVertex3f((cal1(uc, u1)*bsm*contx *bsm.transpose() * cal2(vc, t + 0.5))(0, 0),
						(cal1(uc, u1)*bsm*conty *bsm.transpose() * cal2(vc, t + 0.5))(0, 0),
						(cal1(uc, u1)*bsm*contz *bsm.transpose() * cal2(vc, t + 0.5))(0, 0)
						);


					glVertex3f((cal1(uc, u1 + 0.5)*bsm*contx *bsm.transpose() * cal2(vc, t + 0.5))(0, 0),
						(cal1(uc, u1 + 0.5)*bsm*conty *bsm.transpose() * cal2(vc, t + 0.5))(0, 0),
						(cal1(uc, u1 + 0.5)*bsm*contz *bsm.transpose() * cal2(vc, t + 0.5))(0, 0)
						);

					glVertex3f((cal1(uc, u1+0.5)*bsm*contx *bsm.transpose() * cal2(vc, t ))(0, 0),
						(cal1(uc, u1+0.5)*bsm*conty *bsm.transpose() * cal2(vc, t ))(0, 0),
						(cal1(uc, u1+0.5)*bsm*contz *bsm.transpose() * cal2(vc, t ))(0, 0)
						);
					


					glEnd();
			
				}


			}
		}
	}

	

	if (gouraud)
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_DEPTH_TEST);

		glEnable(GL_NORMALIZE);


		for (oo = 0; oo < xs.size() - 1; oo = oo + 1)
		{
			for (pp = 0; pp < ys.size() - 1; pp = pp + 1)
			{


				if (z5(pp, oo)> 0.1)
				{
					glBegin(GL_POLYGON);
					glNormal3f((y5(pp + 1, oo + 1) - y5(pp, oo)) * (z5(pp + 1, oo) - z5(pp, oo)) - (z5(pp + 1, oo + 1) - z5(pp, oo)) * (y5(pp + 1, oo) - y5(pp, oo)),
						(z5(pp + 1, oo + 1) - z5(pp, oo)) * (x5(pp + 1, oo) - x5(pp, oo)) - (x5(pp + 1, oo + 1) - x5(pp, oo)) * (z5(pp + 1, oo) - z5(pp, oo)),
						(x5(pp + 1, oo + 1) - x5(pp, oo)) * (y5(pp + 1, oo) - y5(pp, oo)) - (y5(pp + 1, oo + 1) - y5(pp, oo)) * (x5(pp + 1, oo) - x5(pp, oo))
						);
					glVertex3f(x5(pp, oo), y5(pp, oo), z5(pp, oo));
					glVertex3f(x5(pp + 1, oo + 1), y5(pp + 1, oo + 1), z5(pp + 1, oo + 1));
					glVertex3f(x5(pp + 1, oo), y5(pp + 1, oo), z5(pp + 1, oo));
					glEnd();
				}


				if (z5(pp, oo)> 0.1)
				{

					glBegin(GL_POLYGON);
					glNormal3f((y5(pp + 1, oo + 1) - y5(pp, oo + 1)) * (z5(pp, oo) - z5(pp, oo + 1)) - (z5(pp + 1, oo + 1) - z5(pp, oo + 1)) * (y5(pp, oo) - y5(pp, oo + 1)),
						(z5(pp + 1, oo + 1) - z5(pp, oo + 1)) * (x5(pp, oo) - x5(pp, oo + 1)) - (x5(pp + 1, oo + 1) - x5(pp, oo + 1)) * (z5(pp, oo) - z5(pp, oo + 1)),
						(x5(pp + 1, oo + 1) - x5(pp, oo + 1)) * (y5(pp, oo) - y5(pp, oo + 1)) - (y5(pp + 1, oo + 1) - y5(pp, oo + 1)) * (x5(pp, oo) - x5(pp, oo + 1))
						);
					glVertex3f(x5(pp, oo), y5(pp, oo), z5(pp, oo));
					glVertex3f(x5(pp, oo + 1), y5(pp, oo + 1), z5(pp, oo + 1));
					glVertex3f(x5(pp + 1, oo + 1), y5(pp + 1, oo + 1), z5(pp + 1, oo + 1));

					glEnd();



					glEnd();
				}

			}
		}


	}


	glFlush();

}


void SetupRC()
{
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1, 3, -1, 3, -4.0, 12.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

}


void Triad(void)
{
	glColor3f(1.0, 1.0, 1.0);

	glBegin(GL_LINES);
	glVertex3fv(ORG); glVertex3fv(XUP);
	glVertex3fv(ORG); glVertex3fv(YUP);
	glVertex3fv(ORG); glVertex3fv(ZUP);
	glEnd();

	glRasterPos3f(1.1, 0.0, 0.0);
	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'X');

	glRasterPos3f(0.0, 1.1, 0.0);
	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'Y');

	glRasterPos3f(0.0, 0.0, 1.1);
	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'Z');
}




void clamp(vec3_t v)
{
	int i;

	for (i = 0; i < 3; i++)
		if (v[i] > 360 || v[i] < -360)
			v[i] = 0.0f;
}



void glutMotion(int x, int y)
{
	if (mState == DOWN)
	{
		gRot[0] -= ((oldY - y) * 180.0f) / 100.0f;
		gRot[1] -= ((oldX - x) * 180.0f) / 100.0f;
		clamp(gRot);
		glutPostRedisplay();
	}
	oldX = x;
	oldY = y;
}


void glutMouse(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN)
	{
		switch (button)
		{
		case GLUT_LEFT_BUTTON:
		case GLUT_RIGHT_BUTTON:
			mState = DOWN;
			oldX = x;
			oldY = y;
			break;
		}
	}
	else if (state == GLUT_UP)
		mState = UP;
}




///////////////////////////////////////////////////////////
// Window has changed size, recalculate projection
void Changeviewport(int w, int h)
{
	lastWidth = w;
	lastHeight = h;

	switch (mAction) {

	case 1:
		if (gouraud == true)
		{

			gouraud = false;
			break;
		}
		gouraud = true;

		break;

	case 2:
		glMatrixMode(GL_MODELVIEW);

		glViewport(150, 150, 1200, 800);
		break;

	case 3:
		glMatrixMode(GL_MODELVIEW);
		glViewport(-150, -150, 1200, 800);
		break;

	}
}



void DoMenu(int value)
{
	mAction = value;
	Changeviewport(lastWidth, lastHeight);
	glutPostRedisplay();
}

