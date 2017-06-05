//
// Created by dilin on 6/3/17.
//

#include "DRLSE_Edge.h"

DRLSE_Edge::DRLSE_Edge(double mu, double lamda, double alpha, double sigma, double timeStep) :
	mu(mu), lamda(lamda), alpha(alpha), sigma(sigma), timeStep(timeStep)
{
	return;
}

DRLSE_Edge::DRLSE_Edge() :
	DRLSE_Edge(0.2, 5, -3, 1.5, 1)
{
	return;
}

Mat DRLSE_Edge::neumannBoundFunc(Mat f)
{
	int nrows = f.rows, ncols = f.cols;
	Mat g = f;

	g.at<double>(0, 0) = g.at <double>(2, 2);
	g.at<double>(0, ncols - 1) = g.at <double>(2, ncols - 3);
	g.at<double>(nrows - 1, 0) = g.at <double>(nrows - 3, 2);
	g.at<double>(nrows - 1, nrows - 1) = g.at <double>(nrows - 3, ncols - 3);

	for (int i = 1; i < ncols; i++)
	{

		g.at<double>(0, i) = g.at <double>(2, i);
		g.at<double>(nrows - 1, i) = g.at <double>(nrows - 3, i);

	}
	for (int i = 1; i < nrows; i++)
	{

		g.at<double>(i, 0) = g.at <double>(i, 2);
		g.at<double>(i, ncols - 1) = g.at <double>(i, ncols - 3);
	}

	return g;
}



Mat DRLSE_Edge::dirac(Mat x, double sigma)
{
	Mat f = x.clone();
	for (int i = 0; i<x.rows; i++)
	{
		for (int j = 0; j<x.cols; j++)
		{
			if (abs(x.at<double>(i, j)) <= sigma)
				f.at<double>(i, j) = (1.0 / (2.0*sigma))*(1 + cos(M_PI*x.at<double>(i, j) / sigma));
			else
				f.at<double>(i, j) = 0;
		}
	}
	return f;
}



Mat DRLSE_Edge::div(Mat nx, Mat ny)
{
	Mat nxx, nxy, nyx, nyy, f;
	gradient(nx, nxx, nxy);
	gradient(ny, nyx, nyy);
	f = nxx + nyy;
	nxx.release();
	nxy.release();
	nyx.release();
	nyy.release();
	return f;
}



Mat DRLSE_Edge::dp(Mat s)
{
	Mat p = s.clone();
	for (int i = 0; i<s.rows; i++)
	{
		for (int j = 0; j<s.cols; j++)
		{
			if (s.at<double>(i, j) == 0)
				p.at<double>(i, j) = 1;
			else if (s.at<double>(i, j) <= 1)
				p.at<double>(i, j) = ((1.0 / (2 * M_PI))*sin(2*M_PI*s.at<double>(i,j)))/ s.at<double>(i, j);
			else
				p.at<double>(i, j) = (s.at<double>(i, j) - 1)/ s.at<double>(i, j);
		}
	}
	return p;
}


Mat DRLSE_Edge::distReg_p2(Mat phi)
{
	Mat dx, dy, dm;
	gradient(phi, dx, dy);
	magnitude(dx, dy, dm);
	Mat dps = dp(dm);
	Mat nx, ny;
	multiply(dps, dx, nx);
	multiply(dps, dy, ny);
	dx.release();
	dy.release();
	dm.release();
	dps.release();
	return div(nx, ny);
}



Mat DRLSE_Edge::edgeT(Mat phi, Mat g, double sigma)
{
	Mat dphix, dphiy, dgx, dgy,dphim;
	gradient(phi, dphix, dphiy);
	gradient(g, dgx, dgy);
	magnitude(dphix, dphiy, dphim);
	dphim = dphim + 1e-10; // to avoid division by zero
	divide(dphix, dphim, dphix);
	divide(dphiy, dphim, dphiy);
	Mat diracPhi = dirac(phi, sigma);
	Mat y1, y2, y3;
	multiply(dphix, dgx, y1);
	multiply(dphiy, dgy, y2);
	multiply(g, div(dphix, dphiy), y3);
	y1 = y1 + y2 + y3;
	multiply(y1, diracPhi, y1);
	dphix.release();
	dphiy.release();
	dgx.release();
	dgy.release();
	dphim.release();
	diracPhi.release();
	y2.release();
	y3.release();

	return y1;

	/*Mat dx, dy, dx2, dy2, dm;
	gradient(phi, dx, dy);
	pow(dx, 2, dx2);
	pow(dy, 2, dy2);
	sqrt(dx2 + dy2, dm);
	divide(dx, dm + 1e-10, dx);
	divide(dy, dm + 1e-10, dy);
	Mat nx, ny;
	multiply(g, dx, nx);
	multiply(g, dy, ny);
	Mat e;
	multiply(dirac(phi, sigma), div(nx, ny), e);
	dx.release();
	dy.release();
	dx2.release();
	dy2.release();
	dm.release();
	nx.release();
	ny.release();
	return e;*/
}



Mat DRLSE_Edge::areaT(Mat phi, Mat g, double sigma)
{
	Mat a;
	multiply(g, dirac(phi, sigma), a);
	return a;
}



void DRLSE_Edge::run(Mat &phi, Mat &g, int iter)
{
	for (int i = 0; i<iter; i++)
	{
		phi = neumannBoundFunc(phi);
		phi = phi + timeStep * (mu*distReg_p2(phi) + lamda*edgeT(phi, g, sigma) + alpha*areaT(phi, g, sigma));
	}
}


void DRLSE_Edge::gradient(Mat &src, Mat &dx, Mat &dy)
{
	Sobel(src, dx, CV_64FC1, 1, 0, 3);
	Sobel(src, dy, CV_64FC1, 0, 1, 3);
}