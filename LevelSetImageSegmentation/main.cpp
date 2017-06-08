#include<opencv2/opencv.hpp>
#include<iostream>
#include "DRLSE_Edge.h"

using namespace std;
using namespace cv;

const string DISPLAY_WINDOW = "Display Window";
const int INNER_ITER = 5;
const int OUTER_ITER = 300;
const double sigma = 0.8;
const double c0 = 2;

bool roi_capture = false;
bool got_roi = false;

bool initial_click = false;
Point pt1, pt2;
Mat cap_img;
//Callback for mousclick event, the x-y coordinate of mouse button-up and button-down 
//are stored in two points pt1, pt2.

void performSegmentation(Mat img){

	Mat imgSmooth,dx,dy,f,g,dx2,dy2;
	DRLSE_Edge drlse;

	img.convertTo(img,CV_64F);
	GaussianBlur(img,imgSmooth,Size(7,7),sigma);
	drlse.gradient(img, dx, dy);
	pow(dx,2.0,dx2);
	pow(dy,2.0,dy2);
	f = dx2 + dy2;
	g = 1.0/(f+1);

	Mat phi = Mat::ones( img.rows, img.cols ,  CV_64F);
	phi = c0 * phi ;

	for ( int i = pt1.y ; i < pt2.y; i ++ )
	{
		for ( int j= pt1.x ; j < pt2.x ; j ++ )
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
		cvWaitKey(50);
	}
	//cvWaitKey(0);
}
void mouse_click(int event, int x, int y, int flags, void *param)
{

	switch(event)
	{
	case CV_EVENT_LBUTTONDOWN:
		{
			std::cout<<"Mouse Pressed"<<std::endl;

			if(!got_roi && !initial_click)
			{
				pt1.x = x;
				pt1.y = y;
				initial_click = true;

			}
			else
			{
				std::cout<<"ROI Already Acquired"<<std::endl;
			}
			break;
		}

		
						
	case CV_EVENT_LBUTTONUP:
		{
			if(!got_roi)
			{
				Mat cl;
				std::cout<<"Mouse LBUTTON Released"<<std::endl;

				pt2.x = x;
				pt2.y = y;

				cl = cap_img.clone();
				Mat roi(cl, Rect(pt1, pt2));
				Mat prev_imgT = roi.clone();
				std::cout<<"PT1"<<pt1.x<<", "<<pt1.y<<std::endl;
				std::cout<<"PT2"<<pt2.x<<","<<pt2.y<<std::endl;

				//imshow("Clone",cl);

				got_roi = true;
				performSegmentation(cap_img);
			}
			else
			{
				std::cout<<"ROI Already Acquired"<<std::endl;
			}
			break;  
		}

	}

}




int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		cout << "Pleae provide an input image path as an argument." << endl;
		return -1;
	}


	Mat img;
	img = imread(argv[1], cv::ImreadModes::IMREAD_GRAYSCALE);
	cap_img = img.clone();

	namedWindow("Original", cv::WINDOW_GUI_NORMAL);
	setMouseCallback("Original",mouse_click);
	
	imshow("Original", cap_img);
	cvWaitKey(0);
	return 0;
}