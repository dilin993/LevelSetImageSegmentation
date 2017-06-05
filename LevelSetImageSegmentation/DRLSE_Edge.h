#ifndef DRLSE_DRLSE_EDGE_H
#define DRLSE_DRLSE_EDGE_H
#include<opencv2/opencv.hpp>
#include<iostream>
#include<cmath>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795



using namespace std;
using namespace cv;

class DRLSE_Edge
{
public:
	DRLSE_Edge(double mu, double lamda, double alpha, double sigma, double timeStep);
	DRLSE_Edge();
	Mat neumannBoundFunc(Mat f);
	Mat dirac(Mat x, double sigma);
	Mat div(Mat nx, Mat ny);
	Mat distReg_p2(Mat phi);
	Mat dp(Mat s);
	Mat edgeT(Mat phi, Mat g, double sigma);
	Mat areaT(Mat phi, Mat g, double sigma);
	void run(Mat &phi, Mat &g, int iter);
	void gradient(Mat &src, Mat &dx, Mat &dy);

private:
	double mu;
	double lamda;
	double alpha;
	double sigma;
	double timeStep;
};

#endif // !M_PI
#endif //DRLSE_DRLSE_EDGE_H
