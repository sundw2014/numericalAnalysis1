#ifndef _IMAGEPROCESS_H
#define _IMAGEPROCESS_H

#include "stddef.h"
#include "stdint.h"
#include "stdio.h"
#include <vector>

#define MAX_UCHAR 255

// cast 函数， 模仿OpenCV
template<class T>
inline int saturateCastUchar(T x) {
  return (x > MAX_UCHAR ? MAX_UCHAR : (x < 0 ? 0 : x));
}

// 下面是自己定义的一些数据类型
typedef uint8_t PixelInt;
typedef size_t PointInt;

typedef union
{
  struct{PixelInt r,g,b;};
  PixelInt data[3];
}PixelRGB;

typedef union
{
  struct{PointInt x,y;};
  PointInt data[2];
}mPoint;

typedef std::vector<mPoint> PointSet;
typedef struct
{
  PointSet source, target;
}ControlPoints;

class IMG_RGB
{
public:
  PixelRGB **data;
  size_t size[2];
  IMG_RGB(size_t _size[2]){
    size[0] = _size[0];
    size[1] = _size[1];
    data = new PixelRGB* [size[0]];
    for(int i=0;i<size[0];i++){
      data[i] = new PixelRGB[size[1]];
    }
  }
  ~IMG_RGB(){
    for(int i=0;i<size[0];i++){
      delete [] data[i];
    }
    delete [] data;
  }
};

template<class T>
class Array2D
{
public:
  T **data;
  size_t size[2];

  Array2D(size_t m, size_t n){
    size[0] = m;
    size[1] = n;

    data = new T* [m];
    for(size_t i=0;i<m;i++){
      data[i] = new T[n];
    }
  }
  ~Array2D(){
    for(size_t i=0;i<size[0];i++){
      delete [] data[i];
    }
    delete [] data;
  }
  T* operator[] (size_t index) const
  {
    return data[index];
  }
  void print()
  {
    for(int i=0;i<size[0];i++){
      for(int j=0;j<size[1];j++){
        printf("%lf,",data[i][j]);
      }
      printf("\r\n");
    }
  }
};

// 插值函数指针类型
typedef PixelRGB InterpolationMethod(const double axis[2], const IMG_RGB& rawI);

void twist(const IMG_RGB& rawI, const IMG_RGB& result, const double theta, const double row, InterpolationMethod interpolation);
void distort(const IMG_RGB& rawI, const IMG_RGB& result, const double k[3], InterpolationMethod interpolation);
void TPSdist(const IMG_RGB& rawI, const IMG_RGB& result, const ControlPoints control, InterpolationMethod interpolation);
PixelRGB NearestNeighborRGB(const double axis[2], const IMG_RGB& rawI);
PixelRGB biLinearRGB(const double axis[2], const IMG_RGB& rawI);
PixelRGB biCubicRGB(const double axis[2], const IMG_RGB& rawI);
inline double biCubicKernel(const double x, const double a = -1.0);

#endif
