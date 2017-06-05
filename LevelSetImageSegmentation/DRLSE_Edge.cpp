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
	/*Mat f = x.clone();
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
	return f;*/
	Mat f, b;
	Mat q = M_PI * x / sigma;
	f = (1.0 / (2.0*sigma))*(1.0 + cos2(q));
	b = (x <= sigma) & (x >= -sigma);
	b.convertTo(b,CV_64F);
	b = b / 255;
	f = f.mul(b);
	b.release();
	q.release();
	return f;
}



Mat DRLSE_Edge::div(Mat nx, Mat ny)
{
	Mat nxx, junk, nyy, f;
	gradient(nx, nxx, junk);
	gradient(ny, junk, nyy);
	f = nxx + nyy;
	nxx.release();
	junk.release();
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
			if (s.at<double>(i, j) < 1e-10)
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
	/*Mat dx, dy, dm;
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
	return div(nx, ny);*/

	Mat phi_x, phi_y,s;
	gradient(phi, phi_x, phi_y);
	magnitude(phi_x, phi_y, s);
	Mat a, b,ps,dps,f,delta2;
	a = (s >= 0) & (s <= 1);
	b = (s > 1);
	a.convertTo(a, CV_64F);
	a = a / 255.0;
	b.convertTo(b, CV_64F);
	b = b / 255.0;
	Mat q = 2.0*M_PI*s;
	ps = a.mul(sin2(q) / (2.0 * M_PI)) + b.mul(s-1.0);
	Mat ps_eq_0, ps_neq_0, s_eq_0, s_neq_0;
	ps_eq_0 = (ps == 0);
	ps_neq_0 = (ps != 0);
	s_eq_0 = (s == 0);
	s_neq_0 = (s != 0);
	ps_eq_0.convertTo(ps_eq_0, CV_64F);
	ps_eq_0 = ps_eq_0 / 255;
	ps_neq_0.convertTo(ps_neq_0, CV_64F);
	ps_neq_0 = ps_neq_0 / 255;
	s_eq_0.convertTo(s_eq_0, CV_64F);
	s_eq_0 = s_eq_0 / 255;
	s_neq_0.convertTo(s_neq_0, CV_64F);
	s_neq_0 = s_neq_0 / 255;
	dps = (ps_neq_0.mul(ps) + ps_eq_0) / (s_neq_0.mul(s) + s_eq_0);
	Laplacian(phi, delta2, CV_64F);
	f = div(dps.mul(phi_x) - phi_x, dps.mul(phi_y) - phi_y) + delta2;
	phi_x.release();
	phi_y.release();
	s.release();
	a.release();
	b.release();
	ps.release();
	dps.release();
	delta2.release();
	q.release();
	ps_eq_0.release();
	ps_neq_0.release();
	s_eq_0.release();
	s_neq_0.release();
	return f;
}



Mat DRLSE_Edge::edgeT(Mat phi, Mat g, double sigma)
{
	//Mat dphix, dphiy, dgx, dgy,dphim;
	//gradient(phi, dphix, dphiy);
	//gradient(g, dgx, dgy);
	//magnitude(dphix, dphiy, dphim);
	//dphim = dphim + 1e-10; // to avoid division by zero
	//divide2(dphix, dphim, dphix);
	//divide2(dphiy, dphim, dphiy);
	//Mat diracPhi = dirac(phi, sigma);
	//Mat y1, y2, y3;
	//multiply(dphix, dgx, y1);
	//multiply(dphiy, dgy, y2);
	//multiply(g, div(dphix, dphiy), y3);
	//y1 = y1 + y2 + y3;
	//multiply(y1, diracPhi, y1);
	//dphix.release();
	//dphiy.release();
	//dgx.release();
	//dgy.release();
	//dphim.release();
	//diracPhi.release();
	//y2.release();
	//y3.release();

	//return y1;

	/*Mat dx, dy, dx2, dy2, dm;
	gradient(phi, dx, dy);
	pow(dx, 2, dx2);
	pow(dy, 2, dy2);
	sqrt(dx2 + dy2, dm);
	divide2(dx, dm + 1e-10, dx);
	divide2(dy, dm + 1e-10, dy);
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

	Mat phi_x, phi_y, vx, vy, s;
	gradient(phi, phi_x, phi_y);
	gradient(g, vx, vy);
	magnitude(phi_x, phi_y, s);
	Mat Nx = phi_x / (s + 1e-10);
	Mat Ny = phi_y / (s + 1e-10);
	Mat curvature = div(Nx, Ny);
	Mat diracPhi = dirac(phi, sigma);
	Mat edgeTerm = diracPhi.mul(vx.mul(Nx) + vy.mul(Ny)) +
		diracPhi.mul(g.mul(curvature));
	return edgeTerm;
}



Mat DRLSE_Edge::areaT(Mat phi, Mat g, double sigma)
{
	Mat a,diracPhi;
	diracPhi = dirac(phi, sigma);
	a = g.mul(diracPhi);
	return a;
}



void DRLSE_Edge::run(Mat &phi, Mat &g, int iter)
{
	Mat distTerm, edgeTerm, areaTerm;
	for (int i = 0; i<iter; i++)
	{
		phi = neumannBoundFunc(phi);
		distTerm = distReg_p2(phi);
		edgeTerm = edgeT(phi, g, sigma);
		areaTerm = areaT(phi, g, sigma);
		imwrite("dist.bmp", distTerm);
		phi = phi + timeStep * (mu*distTerm + lamda*edgeTerm + alpha*areaTerm);
	}
}


void DRLSE_Edge::gradient(Mat &src, Mat &dx, Mat &dy)
{
	/*Sobel(src, dx, CV_64FC1, 1, 0, 3);
	Sobel(src, dy, CV_64FC1, 0, 1, 3);*/
	static cv::Mat kernelx = (cv::Mat_<double>(1, 3) << -0.5, 0, 0.5);
	static cv::Mat kernely = (cv::Mat_<double>(3, 1) << -0.5, 0, 0.5);

	cv::filter2D(src, dx, -1, kernelx, cv::Point(-1, -1), 0, cv::BORDER_REPLICATE);
	cv::filter2D(src, dy, -1, kernely, cv::Point(-1, -1), 0, cv::BORDER_REPLICATE);

	dx.col(dx.cols - 1) *= 2;
	dx.col(0) *= 2;
	dy.row(dy.rows - 1) *= 2;
	dy.row(0) *= 2;
}

void DRLSE_Edge::divide2(Mat & src1, Mat & src2, Mat & ans)
{
	ans = src1.clone();
	for (int i = 0; i < src1.rows; i++)
	{
		for (int j = 0; j < src1.cols; j++)
		{
			ans.at<double>(i, j) = src1.at<double>(i, j) / src2.at<double>(i, j);
		}
	}
}

void DRLSE_Edge::divide2(double val1, Mat & src2, Mat & ans)
{
	ans = src2.clone();
	for (int i = 0; i < src2.rows; i++)
	{
		for (int j = 0; j < src2.cols; j++)
		{
			ans.at<double>(i, j) = val1 / src2.at<double>(i, j);
		}
	}
}

Mat DRLSE_Edge::sin2(Mat & x)
{
	Mat ans = x.clone();
	for (int i = 0;i < x.rows; i++)
	{
		for (int j = 0; j < x.cols; j++)
		{
			ans.at<double>(i, j) = sin(x.at<double>(i, j));
		}
	}
	return ans;
}

Mat DRLSE_Edge::cos2(Mat & x)
{
	Mat ans = x.clone();
	for (int i = 0; i < x.rows; i++)
	{
		for (int j = 0; j < x.cols; j++)
		{
			ans.at<double>(i, j) = cos(x.at<double>(i, j));
		}
	}
	return ans;
}
