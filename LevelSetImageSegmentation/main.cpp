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
int numberOfAreas=0;

bool roi_capture = false;
bool got_roi = false;

bool finished_click = false;
Point pt1, pt2;
Mat cap_img;

vector<pair<Point,Point>>  initAreas;

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

	Point p1,p2;
	for (int k = 0; k < initAreas.size() ; k++){

		p1 = initAreas[k].first;
		p2= initAreas[k].second;
		for ( int i = p1.y ; i < p2.y; i ++ )
		{
			for ( int j= p1.x ; j < p2.x ; j ++ )
			{
				phi.at<double>(i,j) = -c0;
			}
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

}
void mouse_click(int event, int x, int y, int flags, void *param)
{
	switch(event)
	{
	case CV_EVENT_LBUTTONDOWN:
		{
			std::cout<<"Mouse Pressed"<<std::endl;
			if(!got_roi )
			{
				pt1.x = x;
				pt1.y = y;
				numberOfAreas--;

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

				pair < Point, Point> temp (pt1,pt2);
				initAreas.push_back(temp);

				cl = cap_img.clone();
				Mat roi(cl, Rect(pt1, pt2));
				Mat prev_imgT = roi.clone();

				std::cout<<"PT1"<<pt1.x<<", "<<pt1.y<<std::endl;
				std::cout<<"PT2"<<pt2.x<<", "<<pt2.y<< "/n" <<std::endl;

				if(numberOfAreas == 0){
					got_roi = true;
					performSegmentation(cap_img);
				}

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
		cout << "Please provide an input image path as an argument." << endl;
		return -1;
	}

	cout << "Enter number of initial areas" << endl;
	cin >> numberOfAreas;

	Mat img;
	img = imread(argv[1], cv::ImreadModes::IMREAD_GRAYSCALE);
	cap_img = img.clone();

	namedWindow("Original", cv::WINDOW_GUI_NORMAL);
	setMouseCallback("Original",mouse_click);

	imshow("Original", cap_img);
	cvWaitKey(0);
	return 0;
}