#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include "imageprocess.h"
#include "blas.h"
#include "cvutilities.h"

using namespace cv;
using namespace std;

// Mouse events
void onMouse(int event, int x, int y, int flags, void* param);

// TPS vars
static int TPSstate=0;
ControlPoints controlPs;

int main( int argc, char** argv )
{
    if( argc != 2)
    {
     cout <<" Usage: display_image ImageToProcessAndDisplay" << endl;
     return -1;
    }

    Mat rawImage;
    rawImage = imread(argv[1], CV_LOAD_IMAGE_COLOR);   // Read the raw image

    if(!rawImage.data )                              // Check for invalid input
    {
        cout <<  "Could not open the image: " << argv[1] << std::endl ;
        return -1;
    }

    size_t rawSize[2] = {rawImage.rows, rawImage.cols};

    Mat processedImage(rawSize[0], rawSize[1], rawImage.type());

    // 转换为自己定义的图像数据类型
    IMG_RGB processed_IMG_RGB(rawSize), raw_IMG_RGB(rawSize);
    cvMat2IMG_RGB(rawImage,raw_IMG_RGB);

    // 默认插值方法是最近邻
    InterpolationMethod *interpolation = NearestNeighborRGB;

    char key;
    cout<<"select interpolation method:\r\n0. Nearest Neighbor\r\n1. biLinear\r\n2. biCubic\r\n";
    cin>>key;
    switch(key)
    {
      case '0':
      {
        interpolation = NearestNeighborRGB;
        break;
      }
      case '1':
      {
        interpolation = biLinearRGB;
        break;
      }
      case '2':
      {
        interpolation = biCubicRGB;
        break;
      }
      default:
      {
        interpolation = NearestNeighborRGB;
        break;
      }
    }

    cout<<"select mode:\r\n0. twist image\r\n1. distort an image\r\n2. apply TPS on an image\r\n";
    cin>>key;
    switch(key)
    {
      case '0':
      {
        double theta=0.0;
        cout<<"input angle:\r\n";
        cin>>theta;
        cout<<"\r\ntheta="<<theta<<"\r\n";
        int row = (int)((min(rawSize[0], rawSize[1])-1)/2);
        cout<<"row="<<row<<"\r\n";
        twist(raw_IMG_RGB, processed_IMG_RGB, theta/180.0*PI, row, *interpolation);
        break;
      }
      case '1':
      {
        double k[3]={0.5, 0.5, 0.5};
        cout<<"input 3 k-parameters:(like 0.5 0.5 0.5)\r\n";
        cin>>k[0]>>k[1]>>k[2];
        distort(raw_IMG_RGB, processed_IMG_RGB, k, *interpolation);
        break;
      }
      case '2':
      {
        bool startRender=false;
        namedWindow( "raw image Display window", WINDOW_AUTOSIZE );// Create a window for display.
        setMouseCallback("raw image Display window",onMouse,reinterpret_cast<void*> (&rawImage));

        while(1){
          Mat TPSSelectingImage = rawImage.clone();

          for(int i=0;i<min(controlPs.target.size(), controlPs.source.size());i++){
            circle( TPSSelectingImage, Point( controlPs.source[i].x, controlPs.source[i].y ), 3.0, Scalar( 0, 0, 255 ), 1, 8 );
            circle( TPSSelectingImage, Point( controlPs.target[i].x, controlPs.target[i].y ), 3.0, Scalar( 0, 0, 255 ), 1, 8 );
            line(TPSSelectingImage, Point( controlPs.source[i].x, controlPs.source[i].y ), Point( controlPs.target[i].x, controlPs.target[i].y), Scalar( 110, 220, 0 ),  2, 8 );
          }
          if(controlPs.source.size()>controlPs.target.size())
          {
            int i = controlPs.source.size()-1;
            circle( TPSSelectingImage, Point( controlPs.source[i].x, controlPs.source[i].y ), 3.0, Scalar( 0, 0, 255 ), 1, 8 );
          }

          imshow( "raw image Display window", TPSSelectingImage );                   // Show our image inside it.

          char TPSKey = waitKey(10);                                          // Wait for a keystroke in the window
          switch(TPSKey)
          {
            case 'd': //按下'd'后删除一个点
              switch(TPSstate){
                case 0://去除一个target
                  if(controlPs.target.size()>0){
                    controlPs.target.pop_back();
                  }
                  TPSstate = 1;
                  break;
                case 1://去除一个source
                  if(controlPs.source.size()>0){
                    controlPs.source.pop_back();
                  }
                  TPSstate = 0;
                  break;
              }
              break;

            case '\n':
              if(TPSstate == 0){
                startRender = true;
              }
              break;

            default:
              break;
          }
          if(startRender)
          {
            break;
          }
        }
        TPSdist(raw_IMG_RGB, processed_IMG_RGB, controlPs, *interpolation);
      }
    }

    IMG_RGB2cvMat(processed_IMG_RGB, processedImage);
    namedWindow( "processed image Display window", WINDOW_AUTOSIZE );// Create a window for display.


    imshow( "processed image Display window", processedImage);
    waitKey(0);

  return 0;
}

void onMouse(int event, int x, int y, int flags, void* param)
{
    Mat *im = reinterpret_cast<Mat*>(param);
    switch (event)
    {
        case CV_EVENT_LBUTTONDOWN:     //鼠标左键按下响应
            switch(TPSstate){
              case 0://添加source
                controlPs.source.push_back(mPoint({x,y}));
                TPSstate = 1;
                break;
              case 1://添加target
                controlPs.target.push_back(mPoint({x,y}));
                TPSstate = 0;
                break;
            }
            break;
        default:
            break;
    }
}
