#ifndef _IMAGEPROCESS_H
#define _IMAGEPROCESS_H

#include "stddef.h"
#include "stdint.h"
#include "stdio.h"
#include <vector>

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
}Point;

typedef std::vector<Point> PointSet;
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
    // printf("fuck\r\n");
    delete [] data;
  }
};

template<class T>
class 2DArray
{
public:
  T **data;
  size_t size[2];

  2DArray(size_t m, size_t n){
    size[0] = m;
    size[1] = n;

    data = new T*[m];
    for(size_t i=0;i<m;i++){
      data[i] = new T[n];
    }
  }
  ~2DArray(){
    for(size_t i=0;i<m;i++){
      delete [] data[i];
    }
    delete [] data;
  }
  T* operator[] (size_t index){
    return data[index];
  }
};

typedef PixelRGB InterpolationMethod(const double axis[2], const IMG_RGB& rawI);

void twist(const IMG_RGB& rawI, const IMG_RGB& result, const double theta, const double row, InterpolationMethod interpolation);
void distort(const IMG_RGB& rawI, const IMG_RGB& result, const double k[3], InterpolationMethod interpolation);
PixelRGB NearestNeighborRGB(const double axis[2], const IMG_RGB& rawI);
PixelRGB biLinearRGB(const double axis[2], const IMG_RGB& rawI);
PixelRGB biCubicRGB(const double axis[2], const IMG_RGB& rawI);
inline double biCubicKernel(const double x, const double a = -0.6);

#endif
