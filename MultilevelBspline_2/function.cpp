

#include "function.h"

#include <vector>


#define MAX_SIZE 2560

char InputString[MAX_SIZE];


float diff(float x, float y, MatrixXf pz, int d1)
{
	float s, t, g, a, x1, y1;
	int i, j, k, l;

	x1 = x*d1;
	y1 = y*d1;
	s = x1 - floor(x1);
	t = y1 - floor(y1);

	j = floor(y1) - 1;
	i = floor(x1) - 1;

	g = 0;
	for (k = 1; k <= 4; k = k + 1)
	{
		for (l = 1; l <= 4; l = l + 1)
		{

			a = basisf(k, s)*basisf(l, t) * pz(j + l, i + k);
			g = g + a;
		}
	}

	return g;
}





void valuez(int xsize, int ysize, MatrixXf xpos, MatrixXf ypos, MatrixXf & zpos, MatrixXf ct, int d1, float minz)
{

	int oo, pp, k, l;
	float j, i, c, d, s, t, t1, s1;
	c = minz;
	for (oo = 0; oo < xsize; oo = oo + 1)
	{
		for (pp = 0; pp < ysize; pp = pp + 1)
		{

			for (k = 1; k <= 4; k = k + 1)
			{
				for (l = 1; l <= 4; l = l + 1)
				{

					s1 = xpos(pp, oo) * d1;
					s = s1 - floor(s1);
					t1 = ypos(pp, oo) * d1;
					t = t1 - floor(t1);

					i = floor(s1) - 1;
					j = floor(t1) - 1;

					d = basisf(k, s)*basisf(l, t) * ct(j + l, i + k);
					c = d + c;
				}
			}
			zpos(pp, oo) = c;
			c = minz;
		}
	}




}


void refine(MatrixXf cps, MatrixXf & rcps)
{
	for (int i = -1; i < rcps.rows(); i = i + 1)
	{
		for (int j = -1; j < rcps.cols(); j = j + 1)
		{
			if ((i * 2 + 1 >= 0) && (j * 2 + 1 >= 0) && (i * 2 + 1 < rcps.rows()) && (j * 2 + 1 < rcps.cols()))
			{
				rcps(i * 2 + 1, j * 2 + 1) = 1 / 64.0 * (cps(i, j) + cps(i, j + 2) + cps(i + 2, j) + cps(i + 2, j + 2) + 6 * (cps(i, j + 1) + cps(i + 1, j) + cps(i + 1, j + 2) + cps(i + 2, j + 1)) + 36 * cps(i + 1, j + 1));
			}
			if ((i * 2 + 1 >= 0) && (j * 2 + 2 >= 0) && (i * 2 + 1 < rcps.rows()) && (j * 2 + 2 < rcps.cols()))
			{
				rcps(i * 2 + 1, j * 2 + 2) = 1 / 16.0 * (cps(i, j + 1) + cps(i, j + 2) + cps(i + 2, j + 1) + cps(i + 2, j + 2) + 6 * (cps(i + 1, j + 1) + cps(i + 1, j + 2)));
			}
			if ((i * 2 + 2 >= 0) && (j * 2 + 1 >= 0) && (i * 2 + 2 < rcps.rows()) && (j * 2 + 1 < rcps.cols()))
			{
				rcps(i * 2 + 2, j * 2 + 1) = 1 / 16.0 * (cps(i + 1, j) + cps(i + 1, j + 2) + cps(i + 2, j) + cps(i + 2, j + 2) + 6 * (cps(i + 1, j + 1) + cps(i + 2, j + 1)));
			}
			if ((i * 2 + 2 >= 0) && (j * 2 + 2 >= 0) && (i * 2 + 2 < rcps.rows()) && (j * 2 + 2 < rcps.cols()))
			{
				rcps(i * 2 + 2, j * 2 + 2) = 1 / 4.0 * (cps(i + 1, j + 1) + cps(i + 1, j + 2) + cps(i + 2, j + 1) + cps(i + 2, j + 2));
			}

		}

	}



}



void controlvalue(MatrixXf wkl, MatrixXf lot, MatrixXf& pz, int n, int m)
{
	int i, j;

	for (i = 0; i < n + 3; i = i + 1)
	{
		for (j = 0; j < m + 3; j = j + 1)
		{
			if (wkl(i, j) == 0)
			{
				pz(i, j) = 0;
			}

			else
			{
				pz(i, j) = (lot(i, j) / wkl(i, j));
			}

		}
	}


}


MatrixXf controlsetlo(float x, float y, float z, int n, int m, int d)
{
	int i, j, k, l;
	MatrixXf lo(n + 3, m + 3);
	lo.setZero();

	float x1 = x*d;
	float y1 = y*d;

	float s = x1 - floor(x1);
	float t = y1 - floor(y1);

	float c = 0;
	float a = 0;

	for (k = 1; k <= 4; k = k + 1)
	{
		for (l = 1; l <= 4; l = l + 1)
		{
			a = basisf(k, s)*basisf(l, t) * basisf(k, s)*basisf(l, t);
			c = a + c;
		}
	}

	i = floor(x1) - 1;
	j = floor(y1) - 1;


	for (k = 1; k <= 4; k = k + 1)
	{
		for (l = 1; l <= 4; l = l + 1)
		{
			lo(j + l, i + k) = (basisf(k, s)*basisf(l, t) * basisf(k, s)*basisf(l, t) *  (basisf(k, s)*basisf(l, t) * z) / c);
		}
	}

	return lo;
}



MatrixXf controlsetwkl(float x, float y, float z, int n, int m, int d)
{
	int i, j, k, l;
	MatrixXf wkl(n + 3, m + 3);
	wkl.setZero();

	float x1 = x*d;
	float y1 = y*d;

	float s = x1 - floor(x1);
	float t = y1 - floor(y1);

	float c = 0;
	float a = 0;

	for (k = 1; k <= 4; k = k + 1)
	{
		for (l = 1; l <= 4; l = l + 1)
		{
			a = basisf(k, s)*basisf(l, t) * basisf(k, s)*basisf(l, t);
			c = a + c;
		}
	}


	i = floor(x1) - 1;
	j = floor(y1) - 1;


	for (k = 1; k <= 4; k = k + 1)
	{
		for (l = 1; l <= 4; l = l + 1)
		{
			wkl(j + l, i + k) = (basisf(k, s)*basisf(l, t) * basisf(k, s)*basisf(l, t));
		}
	}


	return wkl;


}




float Fileload(char *filename, vector<Orientation>& points_1)
{


	ifstream inFile_1(filename);

	int count_m = 0;
	int count_m2 = 0;
	//vector<Orientation1> points_1;

	while (!inFile_1.eof())
	{

		float xpp1 = 0;
		float ypp1 = 0;
		float zpp1 = 0;

		inFile_1.getline(InputString, 256);
		if (InputString[0] == 'v' && InputString[1] == ' ')
		{
			sscanf_s(InputString, "v %f %f %f", &xpp1, &ypp1, &zpp1);
			count_m++;
			struct Orientation ori;
			ori.xpos = xpp1;
			ori.ypos = ypp1;
			ori.zpos = zpp1;
			points_1.push_back(ori);
		}
	}
	inFile_1.close();

	return 0;
}



float Filesave(char *filename, vector<Orientation>& points_1)
{
	MatrixXd New(3, points_1.size());
	New.setZero();
	MatrixXd Orientation(3, points_1.size());
	for (int i = 0; i<points_1.size(); i++)
	{
		Orientation(0, i) = points_1.at(i).xpos;
		Orientation(1, i) = points_1.at(i).ypos;
		Orientation(2, i) = points_1.at(i).zpos;
	}
	New = Orientation;
	//cout<<"Registrated point="<<New<<endl;
	ofstream outfile(filename);

	for (int j = 0; j<points_1.size(); j++)
	{
		outfile << "v " << New(0, j) << " " << New(1, j) << " " << New(2, j) << endl;
	}

	outfile.close();
	return 0;
}



float Filesave3(char *filename, vector<location> points_1, MatrixXf contpz, MatrixXf number, float value)
{
	MatrixXd New(3, points_1.size());
	New.setZero();
	MatrixXd point(3, points_1.size());
	for (int i = 0; i<points_1.size(); i++)
	{
		point(0, i) = points_1.at(i).x;
		point(1, i) = points_1.at(i).y;
		point(2, i) = points_1.at(i).z;
	}
	New = point;
	//cout<<"Registrated point="<<New<<endl;

	ofstream outfile(filename);

	for (int j = 0; j<points_1.size(); j++)
	{
		outfile << "v " << New(0, j)*100-84.7441-0.1 << " " << -(New(1, j)*100-50.9349-0.1) << " " << -(New(2, j)*100-1713.4) << endl;
	}

	for (int j = 0; j< contpz.rows() - 1; j++)
	{
		for (int i = 0; i < contpz.cols() - 1; i++)
		{
			if ((contpz(j, i) > value) && (contpz(j + 1, i + 1) > value) && (contpz(j, i + 1) > value) && (contpz(j + 1, i) > value))
			{
				outfile << "f " << number(j, i) << " " << number(j + 1, i + 1) << " " << number(j + 1, i) << endl;
			}
		}
	}

	for (int j = 0; j< contpz.rows() - 1; j++)
	{
		for (int i = 0; i < contpz.cols() - 1; i++)
		{
			if ((contpz(j, i) > value) && (contpz(j + 1, i + 1) > value) && (contpz(j, i + 1) > value) && (contpz(j + 1, i) > value))
			{
				outfile << "f " << number(j, i) << " " << number(j, i + 1) << " " << number(j + 1, i + 1) << endl;
			}
		}
	}

	outfile.close();
	return 0;
}





float Filesave2(char *filename, MatrixXf contpz, MatrixXf contpx, MatrixXf contpy, float value)
{
	ofstream outfile(filename);

	for (int j = 0; j<contpz.rows() ; j++)
	{
		for (int i = 0; i < contpz.cols() ; i++)
		{
				outfile << "v " << contpx(j, i) << " " << contpy(j, i) << " " << contpz(j, i) << endl;
		}
	}

	outfile.close();
	return 0;
}