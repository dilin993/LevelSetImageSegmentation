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

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cout << "Pleae provide an input image path as an argument." << endl;
        return -1;
    }

    Mat img,imgSmooth,dx,dy,f,g,dx2,dy2;
	DRLSE_Edge drlse;

    img = imread(argv[1], cv::ImreadModes::IMREAD_GRAYSCALE);
    img.convertTo(img,CV_64FC1);
    
	GaussianBlur(img,imgSmooth,Size(7,7),sigma);
	drlse.gradient(img, dx, dy);
    pow(dx,2.0,dx2);
    pow(dy,2.0,dy2);
    f = dx2 + dy2 + 1;
    drlse.divide2(1.0,f,g);

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
        phi = drlse.run(phi,g,INNER_ITER);
        //drawContours(img,phi,0,Scalar(255,0,0));
        imshow(DISPLAY_WINDOW, phi);
        cout << "Iter: " << i << " is completed." << endl;
        cvWaitKey(100);
    }


    cvWaitKey(0);
    return 0;
}