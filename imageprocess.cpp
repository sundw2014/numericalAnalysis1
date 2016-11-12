#include <cmath>
#include <assert.h>
#include "imageprocess.h"
#include "blas.h"

// 实现了所有的变换和插值函数

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
        if(axis[0] >= 0 && axis[1] >= 0 && axis[0] <= rawI.size[0]-1 && axis[1] <= rawI.size[1]-1){
          result.data[p][q] = interpolation(axis, rawI);
        }
        else
        {
          result.data[p][q].data[0] = 0;
          result.data[p][q].data[1] = 0;
          result.data[p][q].data[2] = 0;
        }
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

      if(axis[0] >= 0 && axis[1] >= 0 && axis[0] <= rawI.size[0]-1 && axis[1] <= rawI.size[1]-1){
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
  if(control.target.size() < 3){
      printf("please choose 3 couples of control points at least.\r\n");
      exit(1);
  }

  const size_t N = control.target.size();

  // 解方程Aw=b
  Array2D<double> A(N+3,N+3), Ainverse(N+3,N+3), b(N+3,1);

  // 填充矩阵A
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      double dx = (double)control.target[i].x - (double)control.target[j].x, dy = (double)control.target[i].y - (double)control.target[j].y;
      double d = sqrt(dx * dx + dy * dy);
      // 防止NaN
      if(std::abs(d)<1e-40){
        A[i][j] = 0.0;
      }
      else{
        A[i][j] = d*d*log(d);
      }
    }
  }
  // 填充A的剩余部分
  for(int j=0;j<N;j++){
    A[N][j] = 1.0;
    A[N+1][j] = (double)control.target[j].x;
    A[N+2][j] = (double)control.target[j].y;
  }
  for(int i=0;i<N;i++){
    A[i][N] = 1.0;
    A[i][N+1] = (double)control.target[i].x;
    A[i][N+2] = (double)control.target[i].y;
  }
  for(int i=N;i<N+3;i++){
    for(int j=N;j<N+3;j++){
      A[i][j] = 0.0;
    }
  }
  // 求逆
  inverseSquareM<double>(A,Ainverse);

  Array2D<double> wx(N+3,1), wy(N+3,1);

  // calculate wx, wx是X坐标做样条插值的权重
  // 填充b，这个b矩阵对应的是X坐标
  for(int i=0;i<N;i++){
    b[i][0] = (double)control.source[i].x - (double)control.target[i].x;
  }
  b[N][0] = b[N+1][0] = b[N+2][0] = 0.0;

  // wx=A^(-1)*b
  mul_MM<double>(Ainverse, b, wx);

  // 计算y坐标
  // 其他同上
  for(int i=0;i<N;i++){
    b[i][0] = (double)control.source[i].y - (double)control.target[i].y;
  }
  b[N][0] = b[N+1][0] = b[N+2][0] = 0.0;

  mul_MM<double>(Ainverse, b, wy);

  // iteration over all pixels
  for(int p=0;p<rawI.size[0];p++){
    for(int q=0;q<rawI.size[1];q++){
      Array2D<double> X(1,N+3),tmpResult(1,1);
      double axis[2];

      // 填充X，X是当前坐标算出来的与control points的差
      for(int j=0;j<N;j++){
        double dx = (double)p - (double)control.target[j].x, dy = (double)q - (double)control.target[j].y;
        double d = sqrt(dx * dx + dy * dy);
        if(std::abs(d)<1e-40){
          X[0][j] = 0.0;
        }
        else{
          X[0][j] = d*d*log(d);
        }
      }
      X[0][N] = 1.0;
      X[0][N+1] = (double)p;
      X[0][N+2] = (double)q;

      //calculate dx， 然后x坐标就等于当前整数坐标加上dx
      mul_MM<double>(X,wx,tmpResult);
      axis[0] = p + tmpResult[0][0];

      //calculate dy
      mul_MM<double>(X,wy,tmpResult);
      axis[1] = q + tmpResult[0][0];

      if(axis[0] >= 0 && axis[1] >= 0 && axis[0] <= rawI.size[0]-1 && axis[1] <= rawI.size[1]-1){
        result.data[p][q] = interpolation(axis, rawI);
      }
      else // 超过边缘限制的就不插值了，直接黑掉
      {
        result.data[p][q].data[0] = 0;
        result.data[p][q].data[1] = 0;
        result.data[p][q].data[2] = 0;
      }
    }
  }
}

/*
* start of some interpolation methods
*/
PixelRGB NearestNeighborRGB(const double axis[2], const IMG_RGB& rawI)
{
  PixelRGB result;
  assert(axis[0] >= 0 && axis[1] >= 0 && axis[0] <= rawI.size[0]-1 && axis[1] <= rawI.size[1]-1);
  result = rawI.data[(size_t)round(axis[0])][(size_t)round(axis[1])];
  return result;
}

PixelRGB biLinearRGB(const double axis[2], const IMG_RGB& rawI)
{
  PixelRGB result;
  assert(axis[0] >= 0 && axis[1] >= 0 && axis[0] <= rawI.size[0]-1 && axis[1] <= rawI.size[1]-1);

  size_t i = floor(axis[0]), j = floor(axis[1]);
  double u = axis[0] - i;
  double v = axis[1] - j;
  for(int c=0;c<3;c++){
    result.data[c] = saturateCastUchar(((1-u) * rawI.data[i][j].data[c]\
     + u * rawI.data[i+1][j].data[c])*(1-v)\
     + ((1-u) * rawI.data[i][j+1].data[c]\
     + u * rawI.data[i+1][j+1].data[c])*v);
  }
  return result;
}

PixelRGB biCubicRGB(const double axis[2], const IMG_RGB& rawI)
{
  PixelRGB result;
  assert(axis[0] >= 0 && axis[1] >= 0 && axis[0] <= rawI.size[0]-1 && axis[1] <= rawI.size[1]-1);

  PixelRGB** data = rawI.data;

  size_t i = floor(axis[0]), j = floor(axis[1]);
  double u = axis[0] - i;
  double v = axis[1] - j;

  Array2D<double> A(1,4), CT(4,1), B(4,4), AB(1,4);
  for(int i=0;i<4;i++){
    A[0][i] = biCubicKernel(u+1-i);
    CT[i][0] = biCubicKernel(v+1-i);
  }

  for(int c=0;c<3;c++){
    for(size_t p=0;p<4;p++){
      for(size_t q=0;q<4;q++){
        int x = (int)(i-1+p), y = (int)(j-1+q);
        size_t index[2] = {(x < 0)?0:((x>=rawI.size[0])?rawI.size[0]-1:x) ,(y < 0)?0:((y>=rawI.size[1])?rawI.size[1]-1:y)};
        B[p][q] = data[index[0]][index[1]].data[c];
      }
    }

    Array2D<double> tmpResult(1,1);
    mul_MM<double>(A,B,AB);
    mul_MM<double>(AB,CT,tmpResult);
    // cast一下是为了防止下溢出和上溢出
    result.data[c] = saturateCastUchar(tmpResult[0][0]);
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
* end of some interpolation methods
*/
