#include "cvutilities.h"

using namespace cv;
using namespace std;

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
      m.at<Vec3b>(q,p) = Vec3b(i.data[p][q].data[0], i.data[p][q].data[1], i.data[p][q].data[2]);
    }
  }
}

void cvMat2IMG_RGB(const Mat& m, IMG_RGB& i)
{
  assert(m.rows == i.size[0] && m.cols == i.size[1]);

  for(int p=0;p<i.size[0];p++){
    for(int q=0;q<i.size[1];q++){
      for(int c=0;c<3;c++){
        i.data[p][q].data[c] = m.at<Vec3b>(q,p).val[c];
      }
    }
  }
}
