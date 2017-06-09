#include<opencv2/opencv.hpp>
#include<iostream>
#include "DRLSE_Edge.h"

using namespace std;
using namespace cv;

const string DISPLAY_WINDOW = "Display Window";
const string TRACKBAR_GAMMA = "Gamma Correction";
int gamma;
const int GAMMA_MAX = 100;
const int INNER_ITER = 5;
const int OUTER_ITER = 20;
const double sigma = 0.8;
const double c0 = 2;
Mat img, imgC;

void GammaCorrection(Mat& src, Mat& dst, float fGamma)
{
	CV_Assert(src.data);

	// accept only char type matrices
	CV_Assert(src.depth() != sizeof(uchar));

	// build look up table
	unsigned char lut[256];
	for (int i = 0; i < 256; i++)
	{
		lut[i] = saturate_cast<uchar>(pow((float)(i / 255.0), fGamma) * 255.0f);
	}

	dst = src.clone();
	const int channels = dst.channels();
	switch (channels)
	{
	case 1:
	{

		MatIterator_<uchar> it, end;
		for (it = dst.begin<uchar>(), end = dst.end<uchar>(); it != end; it++)
			//*it = pow((float)(((*it))/255.0), fGamma) * 255.0;
			*it = lut[(*it)];

		break;
	}
	case 3:
	{

		MatIterator_<Vec3b> it, end;
		for (it = dst.begin<Vec3b>(), end = dst.end<Vec3b>(); it != end; it++)
		{

			(*it)[0] = lut[((*it)[0])];
			(*it)[1] = lut[((*it)[1])];
			(*it)[2] = lut[((*it)[2])];
		}

		break;

	}
	}
}

void on_trackbar(int, void*)
{
	float gammaf = 3.0*(float)gamma / (float)GAMMA_MAX;
	GammaCorrection(img, imgC, gammaf);
	cout << "Gamma: " << gammaf << endl;
	imshow(DISPLAY_WINDOW, imgC);
}

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		cout << "Pleae provide an input image path as an argument." << endl;
		return -1;
	}

	img = imread(argv[1], cv::ImreadModes::IMREAD_GRAYSCALE);
	namedWindow(DISPLAY_WINDOW, WINDOW_AUTOSIZE);
	imshow(DISPLAY_WINDOW, img);
	createTrackbar(TRACKBAR_GAMMA, DISPLAY_WINDOW, &gamma, GAMMA_MAX, on_trackbar);
	on_trackbar(gamma, 0);
	waitKey(0);
	return 0;

	/*Mat img,imgSmooth,dx,dy,f,g,dx2,dy2;
	DRLSE_Edge drlse;

	img = imread(argv[1], cv::ImreadModes::IMREAD_GRAYSCALE);
	img.convertTo(img,CV_64F);
	GaussianBlur(img,imgSmooth,Size(7,7),sigma);
	drlse.gradient(img, dx, dy);
	pow(dx,2.0,dx2);
	pow(dy,2.0,dy2);
	f = dx2 + dy2;
	g = 1.0/(f+1);

	Mat phi = Mat::ones( img.rows, img.cols ,  CV_64F);
	phi = c0 * phi ;

	for ( int i = 24 ; i < 35 ; i ++ )
	{
		for ( int j= 19 ; j < 25 ; j ++ )
		{
			phi.at<double>(i,j) = -c0;
		}

		for ( int j= 39 ; j < 50 ; j ++ )
		{
			phi.at<double>(i,j) = -c0;
		}
	}

	namedWindow(DISPLAY_WINDOW, cv::WINDOW_NORMAL);

	for(int i=0;i<OUTER_ITER;i++)
	{
		drlse.run(phi,g,INNER_ITER);
		imshow(DISPLAY_WINDOW, phi);
		cout << "Iter: " << i << " is completed." << endl;
		cvWaitKey(1000);
	}
	cvWaitKey(0);
	return 0;*/
}