#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include "imageprocess.h"
#include "blas.h"

using namespace cv;
using namespace std;

void cvMat2IMG_RGB(const Mat& m, IMG_RGB& i);
void IMG_RGB2cvMat(const IMG_RGB& i, Mat& m);
string type2str(int type);

int main( int argc, char** argv )
{
    // double A[1][4]={{1,2,3,4}}, B[4][1]={{5},{6},{7},{8}}, AB[1][1];
    // mul_MM<double,1,4,1>(A,B,AB);
    // printf("%lf\r\n",AB[0][0]);
    // double A[4][4]={{1,35,23,4},{15,6,7,8},{9,1,11,12},{13,0,20,16}}, Ainverse[4][4];
    // inverseSquareM<double, 4>(A, Ainverse);
    // for(int i=0;i<4;i++){
    //   for(int j=0;j<4;j++){
    //     printf("%lf,",Ainverse[i][j]);
    //   }
    //   printf("\r\n");
    // }
    if( argc != 2)
    {
     cout <<" Usage: display_image ImageToProcessAndDisplay" << endl;
     return -1;
    }

    Mat rawImage;
    rawImage = imread(argv[1], CV_LOAD_IMAGE_COLOR);   // Read the file
    size_t rawSize[2] = {rawImage.rows, rawImage.cols};
    Mat processedImage(rawSize[0], rawSize[1], rawImage.type());

    IMG_RGB processed_IMG_RGB(rawSize), raw_IMG_RGB(rawSize);
    cvMat2IMG_RGB(rawImage,raw_IMG_RGB);
    // twist(raw_IMG_RGB, processed_IMG_RGB, (double)90/180*PI, 250, biCubicRGB);
    double k[3]={0.5, 0.5, 0.5};
    distort(raw_IMG_RGB, processed_IMG_RGB, k, biCubicRGB);
    IMG_RGB2cvMat(processed_IMG_RGB, processedImage);
    //IMG_RGB2cvMat(raw_IMG_RGB, processedImage);

    if(! rawImage.data )                              // Check for invalid input
    {
        cout <<  "Could not open the image: " << argv[1] << std::endl ;
        return -1;
    }

    // string ty =  type2str( rawImage.type() );
    // printf("Matrix: %s %dx%d \n", ty.c_str(), rawImage.cols, rawImage.rows );

    namedWindow( "raw image Display window", WINDOW_AUTOSIZE );// Create a window for display.
    namedWindow( "processed image Display window", WINDOW_AUTOSIZE );// Create a window for display.
    imshow( "raw image Display window", rawImage );                   // Show our image inside it.
    imshow( "processed image Display window", processedImage );

    waitKey(0);                                          // Wait for a keystroke in the window
    return 0;
}

string type2str(int type) {
  string r;

  uchar depth = type & CV_MAT_DEPTH_MASK;
  uchar chans = 1 + (type >> CV_CN_SHIFT);

  switch ( depth ) {
    case CV_8U:  r = "8U"; break;
    case CV_8S:  r = "8S"; break;
    case CV_16U: r = "16U"; break;
    case CV_16S: r = "16S"; break;
    case CV_32S: r = "32S"; break;
    case CV_32F: r = "32F"; break;
    case CV_64F: r = "64F"; break;
    default:     r = "User"; break;
  }

  r += "C";
  r += (chans+'0');

  return r;
}

void IMG_RGB2cvMat(const IMG_RGB& i, Mat& m){
  assert(m.rows == i.size[0] && m.cols == i.size[1]);

  for(int p=0;p<i.size[0];p++){
    for(int q=0;q<i.size[1];q++){
      m.at<Vec3b>(p,q) = Vec3b(i.data[p][q].data[0], i.data[p][q].data[1], i.data[p][q].data[2]);
    }
  }
}

void cvMat2IMG_RGB(const Mat& m, IMG_RGB& i)
{
  assert(m.rows == i.size[0] && m.cols == i.size[1]);

  for(int p=0;p<i.size[0];p++){
    for(int q=0;q<i.size[1];q++){
      for(int c=0;c<3;c++){
        i.data[p][q].data[c] = m.at<Vec3b>(p,q).val[c];
      }
    }
  }
}
