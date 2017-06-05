//
// Created by dilin on 6/3/17.
//

#include "DRLSE_Edge.h"

DRLSE_Edge::DRLSE_Edge(double mu, double lamda, double alpha, double sigma, double timeStep) :
	mu(mu), lamda(lamda), alpha(alpha), sigma(sigma), timeStep(timeStep)
{
	return;
}

DRLSE_Edge::DRLSE_Edge() 

{
	DRLSE_Edge(0.2, 5, -3, 1.5, 1);
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
			//if (abs(x.at<double>(i, j)) <= sigma)
				f.at<double>(i, j) = (1.0 / (2.0*sigma))*(1 + cos(M_PI*x.at<double>(i, j) / sigma));
			//else f.at<double>(i, j) = 0;*/


		}
	}

	Mat b = ( x <= sigma ) & (x >= -sigma) ;
	Mat ret ;
	multiply( f, b , ret ,1, CV_64FC1);
	f.release();
	return ret;
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
			/*if (s.at<double>(i, j) == 0)
			p.at<double>(i, j) = 1;
			else if (s.at<double>(i, j) <= 1)
			p.at<double>(i, j) = ((1.0 / (2 * M_PI))*sin(2*M_PI*s.at<double>(i,j)))/ s.at<double>(i, j);
			else
			p.at<double>(i, j) = (s.at<double>(i, j) - 1)/ s.at<double>(i, j);
			*/
			p.at<double>(i, j) = sin ( 2 * M_PI * s.at<double>(i,j))/ 2* M_PI;
		}
	}
	return p;
}


Mat DRLSE_Edge::distReg_p2(Mat phi)
{
	Mat dx, dy, dm, dx2 , dy2 ,s;
	gradient(phi, dx, dy);
	pow(dx,2,dx2);
	pow(dy,2,dy2);
	sqrt(dx2+dy2, s);

	Mat a = (s >= 0) & (s<=1);
	Mat b = s > 1;
	Mat sinMat = dp(s);
	Mat ps ,temp, t1, t2;
	a.convertTo(a,CV_64FC1);
	b.convertTo(b,CV_64FC1);

	multiply ( a/255.0 , sinMat , ps, 1,CV_64FC1);
	multiply ( b/255.0, s-1 , temp, 1, CV_64FC1);
	ps = (ps + temp );
	ps = ps/255;
	
	multiply((ps!=0),ps,t1 , 1 ,CV_64FC1);
	add( t1 , (ps == 0) , t1 ,noArray(), CV_64FC1 );
	multiply((s!=0),s,t2,1 , CV_64FC1);
	
	add( t2 , (ps == 0) , t2,noArray(), CV_64FC1 );
	Mat dps = t1 / t2;

	Mat nx, ny, lap;
	Laplacian( phi , lap , CV_64FC1);
	multiply( dps, dx ,nx );
	multiply( dps, dy ,ny);
	nx = nx - dx;
	ny = ny - dy;
	
	dx.release();
	dy.release();
	dx2.release();
	dy2.release();
	dx.release();
	s.release();
	t1.release();
	t2.release();
	dps.release();
	ps.release();
	temp.release();
	return div(nx, ny) + lap;

}



Mat DRLSE_Edge::edgeT(Mat phi, Mat g, double sigma)
{
	Mat dx, dy, dx2, dy2, dm, vx, vy;
	gradient(g, vx, vy);
	gradient(phi, dx, dy);


	pow(dx, 2, dx2);
	pow(dy, 2, dy2);
	sqrt(dx2 + dy2, dm);
	//divide(dx, dm + 1e-10, dx);
	//divide(dy, dm + 1e-10, dy);

	Mat nx, ny;

	nx = dx / (dm + 1e-10);
	ny = dy / (dm + 1e-10);

	Mat curvature = div(nx,ny);
	Mat diracPhi = dirac(phi, sigma);
	Mat t1,t2,t3,t4 ,t5;

	multiply(vx,nx,t1);
	multiply(vy,ny,t2);
	multiply(diracPhi,t1+t2, t3);
	multiply(diracPhi,g, t4);
	multiply(t4,curvature, t5);

	dx.release();
	dy.release();
	dx2.release();
	dy2.release();
	dm.release();
	nx.release();
	ny.release();
	t1.release();
	t2.release();
	t4.release();
	curvature.release();
	diracPhi.release();

	return t3+t5;

}



Mat DRLSE_Edge::areaT(Mat phi, Mat g, double sigma)
{
	Mat a;
	multiply(g, dirac(phi, sigma), a);
	return a;
}



Mat DRLSE_Edge::run(Mat phi, Mat &g, int iter)
{
	Mat distTerm, edgeTerm, areaTerm;
	for (int i = 0; i<iter; i++)
	{
		phi = neumannBoundFunc(phi);
		distTerm = distReg_p2(phi);
		edgeTerm = edgeT(phi, g, sigma);
		areaTerm = areaT(phi, g, sigma);
		phi = phi + timeStep * (mu*distTerm + lamda*edgeTerm + alpha*areaTerm);

		imwrite("dist.bmp", distTerm);
	}
	return phi;
}


void DRLSE_Edge::gradient(Mat &src, Mat &dx, Mat &dy)
{
	Sobel(src, dx, CV_64FC1, 1, 0, 3);
	Sobel(src, dy, CV_64FC1, 0, 1, 3);
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
