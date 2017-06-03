#include<opencv2\opencv.hpp>
#include<iostream>

using namespace std;

const string DISPLAY_WINDOW = "Display Window";

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		cout << "Pleae provide an input image path as an argument." << endl;
		return -1;
	}

	cv::Mat img;
	img = cv::imread(argv[1], cv::ImreadModes::IMREAD_GRAYSCALE);

	cv::namedWindow(DISPLAY_WINDOW, cv::WINDOW_NORMAL);

	cv::imshow(DISPLAY_WINDOW, img);

	cvWaitKey(0);
	return 0;
}