#include <cmath>
#include <assert.h>
#include "imageprocess.h"
#include "blas.h"

/*
* twist an image
*/
void twist(const IMG_RGB& rawI, const IMG_RGB& result, const double theta, const double row, InterpolationMethod interpolation)
{
  assert(rawI.size[0] == result.size[0] && rawI.size[1] == result.size[1] && rawI.size[0]>1 && rawI.size[1]>1);
  double offset[2] = {((double)rawI.size[0]-1)/2, ((double)rawI.size[1]-1)/2};

  for(int p=0;p<rawI.size[0];p++){
    for(int q=0;q<rawI.size[1];q++){
      double x = p-offset[0];
      double y = q-offset[1];
      double r = sqrt(x*x + y*y);
      if(r<=row)
      {
        double alpha,garma;
        garma = acos(x/r);
        garma = y<0 ? 2*PI-garma : garma;
        alpha = garma - theta * (row-r) / row;
        double axis[2] = {offset[0]+(double)r*cos(alpha),offset[1]+(double)r*sin(alpha)};
        result.data[p][q] = interpolation(axis, rawI);
      }
      else
      {
        result.data[p][q] = rawI.data[p][q];
      }
    }
  }
}

/*
* distort an image
*/
void distort(const IMG_RGB& rawI, const IMG_RGB& result, const double k[3], InterpolationMethod interpolation)
{
  assert(rawI.size[0] == result.size[0] && rawI.size[1] == result.size[1] && rawI.size[0]>1 && rawI.size[1]>1);

  double offset[2] = {((double)rawI.size[0]-1)/2, ((double)rawI.size[1]-1)/2};

  double row = sqrt(offset[0]*offset[0] + offset[1]*offset[1]);

  for(int p=0;p<rawI.size[0];p++){
    for(int q=0;q<rawI.size[1];q++){
      double x = p-offset[0];
      double y = q-offset[1];
      double r = sqrt(x*x + y*y) / row;

      double r_distorted = r*(1 + r*r*(k[0] + r*r*(k[1] + k[2]*r*r)));
      r_distorted *= row;

      double alpha;
      alpha = acos(x/(r*row));
      alpha = y<0 ? 2*PI-alpha:alpha;
      double axis[2] = {offset[0]+(double)r_distorted*cos(alpha),offset[1]+(double)r_distorted*sin(alpha)};

      if(axis[0] >= 0 && axis[1] >= 0 && axis[0] < rawI.size[0] && axis[1] < rawI.size[1]){
        result.data[p][q] = interpolation(axis, rawI);
      }
      else
      {
        result.data[p][q].data[0] = 0;
        result.data[p][q].data[1] = 0;
        result.data[p][q].data[2] = 0;
      }
    }
  }
}

/*
* apply TPS on an image
*/
void TPSdist(const IMG_RGB& rawI, const IMG_RGB& result, const ControlPoints control, InterpolationMethod interpolation)
{
  assert(rawI.size[0] == result.size[0] && rawI.size[1] == result.size[1] && rawI.size[0]>1 && rawI.size[1]>1);
  assert(control.target.size() == control.source.size());

  const size_t N = control.target.size();

  2DArray<double> A(N+3,N+3), Ainverse(N+3,N+3), b(N+3,1);

  // fill content of A
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      double dx = control.target[i].x - control.target[j].x, dy = control.target[i].y - control.target[j].y;
      double d = sqrt(dx * dx + dy * dy);
      A[i][j] = d*d*log(d);
    }
  }
  for(int j=0;j<N;j++){
    A[N][j] = 1;
    A[N+1][j] = control.target[j].x;
    A[N+2][j] = control.target[j].y;
  }
  for(int i=0;i<N;i++){
    A[i][N] = 1;
    A[i][N+1] = control.target[i].x;
    A[i][N+2] = control.target[i].y;
  }
  for(int i=N;i<N+3;i++){
    for(int j=N;j<N+3;j++){
      A[i][j] = 0.0;
    }
  }

  2DArray<double> wx(N+3,1), wy(N+3,1);

  //calculate wx
  //fill content of b
  for(int i=0;i<N;i++){
    b[i][1] = control.source[i].x - control.target[i].x;
  }
  b[N][1] = b[N+1][1] = b[N+2][1] = 0.0;

  inverseSquareM(A,Ainverse);


  double offset[2] = {((double)rawI.size[0]-1)/2, ((double)rawI.size[1]-1)/2};

  double row = sqrt(offset[0]*offset[0] + offset[1]*offset[1]);

  // for(int p=0;p<rawI.size[0];p++){
  //   for(int q=0;q<rawI.size[1];q++){
  //     double x = p-offset[0];
  //     double y = q-offset[1];
  //     double r = sqrt(x*x + y*y) / row;
  //
  //     double r_distorted = r*(1 + r*r*(k[0] + r*r*(k[1] + k[2]*r*r)));
  //     r_distorted *= row;
  //
  //     double alpha;
  //     alpha = acos(x/(r*row));
  //     alpha = y<0 ? 2*PI-alpha:alpha;
  //     double axis[2] = {offset[0]+(double)r_distorted*cos(alpha),offset[1]+(double)r_distorted*sin(alpha)};
  //
  //     if(axis[0] >= 0 && axis[1] >= 0 && axis[0] < rawI.size[0] && axis[1] < rawI.size[1]){
  //       result.data[p][q] = interpolation(axis, rawI);
  //     }
  //     else
  //     {
  //       result.data[p][q].data[0] = 0;
  //       result.data[p][q].data[1] = 0;
  //       result.data[p][q].data[2] = 0;
  //     }
  //   }
  // }
}

/*
* start of some interpolation methods
*/
PixelRGB NearestNeighborRGB(const double axis[2], const IMG_RGB& rawI)
{
  PixelRGB result;
  assert(axis[0] >= 0 && axis[1] >= 0 && axis[0] < rawI.size[0] && axis[1] < rawI.size[1]);
  result = rawI.data[(size_t)round(axis[0])][(size_t)round(axis[1])];
  return result;
}

PixelRGB biLinearRGB(const double axis[2], const IMG_RGB& rawI)
{
  PixelRGB result;
  assert(axis[0] >= 0 && axis[1] >= 0 && axis[0] < rawI.size[0] && axis[1] < rawI.size[1]);

  size_t i = floor(axis[0]), j = floor(axis[1]);
  double u = axis[0] - i;
  double v = axis[1] - j;
  for(int c=0;c<3;c++){
    result.data[c] = ((1-u) * rawI.data[i][j].data[c] + u * rawI.data[i+1][j].data[c])*(1-v)+((1-u) * rawI.data[i][j+1].data[c] + u * rawI.data[i+1][j+1].data[c])*v;
  }
  return result;
}

PixelRGB biCubicRGB(const double axis[2], const IMG_RGB& rawI)
{
  PixelRGB result;
  assert(axis[0] >= 0 && axis[1] >= 0 && axis[0] < rawI.size[0] && axis[1] < rawI.size[1]);

  PixelRGB** data = rawI.data;

  size_t i = floor(axis[0]), j = floor(axis[1]);
  double u = axis[0] - i;
  double v = axis[1] - j;

  double A[1][4] = {{biCubicKernel(u+1), biCubicKernel(u), biCubicKernel(u-1), biCubicKernel(u-2)}};
  double CT[4][1] = {{biCubicKernel(v+1)},{biCubicKernel(v)},{biCubicKernel(v-1)},{biCubicKernel(v-2)}};
  double B[4][4];
  double AB[1][4];
  for(int c=0;c<3;c++){
    for(size_t p=0;p<4;p++){
      for(size_t q=0;q<4;q++){
        int x = (int)(i-1+p), y = (int)(j-1+q);
        size_t index[2] = {(x < 0)?0:((x>=rawI.size[0])?rawI.size[0]-1:x) ,(y < 0)?0:((y>=rawI.size[1])?rawI.size[1]-1:y)};
        B[p][q] = data[index[0]][index[1]].data[c];
      }
    }
    double tmpResult[1][1]={0};
    mul_MM<double,1,4,4>(A,B,AB);
    mul_MM<double,1,4,1>(AB,CT,tmpResult);
    result.data[c] = (int)tmpResult[0][0];
    // result.data[c] = data[i][j].data[c];
  }
  return result;
}

inline double biCubicKernel(const double x, const double a)
{
  double xabs = std::abs(x);
  if(xabs<=1){
    return (a+2)*xabs*xabs*xabs - (a+3)*xabs*xabs + 1.0;
  }
  else if(1< xabs && xabs < 2){
    return a*xabs*xabs*xabs - 5*a*xabs*xabs + 8*a*xabs - 4.0*a;
  }
  else{
    return 0.0;
  }
}

/*
* end of some interpolation methosd
*/
