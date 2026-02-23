#include <iostream>
#include <opencv2/opencv.hpp>
//#include<typeinfo>

#include "CT_main.h"

//#include <opencv2/highgui/highgui_c.h>

using namespace std;
using namespace cv;

int main() {
  //int radon_flag = radonFlags::analytical;  // for shepplogan
  int radon_flag = radonFlags::integral;  // integral:for existing image

  int radon_type = radonTypes::parallel;
  //int radon_type = radonTypes::sector;
  string image_name = "shepplogan_original";
  Mat radon_Img = Mat::zeros(256, 256, COLOR_BGR2GRAY);
  Mat Img;

  ///*
  if (radon_type == radonTypes::parallel) {

    //
    if (radon_flag == radonFlags::analytical) {
      image_name = "shepplogan_analytical";
      int rows_of_shepplogan = 256;
      //
      Img = Create_shepplogan(rows_of_shepplogan);
      write_and_show_Img(Img, image_name, "original");

      //
      radon_Img = Calc_shepplogan_radon(rows_of_shepplogan);
      write_and_show_Img(radon_Img, image_name, "radon");
    }

    //
    if (radon_flag == radonFlags::integral) {
      //Mat inputImg = Mat::zeros(256,256, COLOR_BGR2GRAY);

      int rows_of_shepplogan = 256;
      Img = Create_shepplogan(rows_of_shepplogan);
      write_and_show_Img(Img, image_name, "original");

      Img = cv::imread("./pics/" + image_name + ".png");
      if (Img.data == nullptr) {
        cout << "ͼƬ�ļ�������" << endl;
        return 0;
      }
      cvtColor(Img, Img, COLOR_BGR2GRAY);
      cout << "image type = " << Img.type() << " ,channel = " << Img.channels()
           << " ,size = " << Img.size() << endl;
      //showImg("original image", image);
      //Mat radon_Img = Mat::zeros(256, 256, COLOR_BGR2GRAY);
      radon_Img = Radon_parallel_line(Img);  //�Լ�д��ƽ����ͶӰ����
      write_and_show_Img(radon_Img, image_name, "radon");
    }

    /*
		Mat radon_Img;
		if (LoadData("./txt/" + image_name + "_radon.txt", radon_Img,365,180,1) == 0)
		{
			cout << "load txt to mat SUCCESS!" << endl;
		}
		radon_Img.convertTo(radon_Img, COLOR_BGR2GRAY);

	*/

    //Mat direct_iradon_Img = iRadon_parallel_line(radon_Img);
    //write_and_show_Img(direct_iradon_Img, image_name, "direct_iradon");

    Mat filtered_radon_Img = Filter_radon_parallel_line(radon_Img, SL_filter);
    write_and_show_Img(filtered_radon_Img, image_name, "filtered_radon_SL");
    Mat filtered_iradon_Img = iRadon_parallel_line(filtered_radon_Img);
    write_and_show_Img(filtered_iradon_Img, image_name, "filtered_iradon_SL");

    double error_d = Calculate_error_sqrt(Img, filtered_iradon_Img);
    cout << "sqrt error d = " << error_d << endl;
    double error_r = Calculate_error_abs(Img, filtered_iradon_Img);
    cout << "sqrt error r = " << error_r << endl;

  } else  //radon_type == radonTypes::sector;
  {
    //
    if (radon_flag == radonFlags::analytical) {
      //
    }

    //
    if (radon_flag == radonFlags::integral) {
      //����ͼƬ
      cv::Mat image;
      image = cv::imread("./pics/" + image_name + ".png");
      if (image.data == nullptr) {
        cout << "ͼƬ�ļ�������" << endl;
        return 0;
      }
      cvtColor(image, image, COLOR_BGR2GRAY);
      cout << "image type = " << image.type()
           << " ,channel = " << image.channels() << " ,size = " << image.size()
           << endl;
      //showImg("original image", image);

      radon_Img = Radon_fan_beam(image, 360, 90);
      write_and_show_Img(radon_Img, image_name, "radon_fanbeam");
    }
  }

  //
  //waitKey(0);

  return 0;
}
