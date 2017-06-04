//
// Created by dilin on 6/4/17.
//

#include<opencv2/opencv.hpp>
#include<iostream>
#include "DRLSE_Edge.h"

using namespace std;
using namespace cv;

const string DISPLAY_WINDOW = "Display Window";
const int INNER_ITER = 5;
const int OUTER_ITER = 20;
const double sigma = 0.8;
const double c0 = 2;
const double epsilon = 1.5;
int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		cout << "Pleae provide an input image path as an argument." << endl;
		return -1;
	}

	Mat img, imgSmooth, dx, dy, f, g;

	img = imread(argv[1], cv::ImreadModes::IMREAD_GRAYSCALE);
	//img.convertTo(img,CV_64F);
	GaussianBlur(img, imgSmooth, Size(7, 7), sigma);
	spatialGradient(imgSmooth, dx, dy);
	dx.convertTo(dx, CV_64F);
	dx.convertTo(dy, CV_64F);
	pow(dx, 2.0, dx);
	pow(dy, 2.0, dy);
	f = dx + dy;
	divide(1.0, f, g);

	Mat phi = Mat::ones(img.rows, img.cols, CV_64F);
	phi = c0 * phi;



	for (int i = 24; i < 35; i++)
	{

		for (int j = 19; j < 25; j++)
		{

			phi.at<double>(i, j) = -c0;

		}

		for (int j = 39; j < 50; j++)
		{

			phi.at<double>(i, j) = -c0;

		}
	}

	namedWindow("dx", WINDOW_NORMAL);
	namedWindow("dy", WINDOW_NORMAL);
	namedWindow("dirac", WINDOW_NORMAL);

	DRLSE_Edge drlse;

	Mat ddx, ddy, dirac;
	drlse.gradient(img, ddx, ddy);

	dirac = drlse.edgeT(phi, g, epsilon);
	cout << sum(dirac) << endl;
	ddx.convertTo(ddx, CV_8U);
	ddy.convertTo(ddy, CV_8U);
	dirac.convertTo(dirac, CV_8U);
	imshow("dx", ddx);
	imshow("dy", ddy);
	imshow("dirac", dirac);
	cvWaitKey(0);
	return 0;
}