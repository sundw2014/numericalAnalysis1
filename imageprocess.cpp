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

  const size_t N = control.target.size();

  Array2D<double> A(N+3,N+3), Ainverse(N+3,N+3), b(N+3,1);

  // fill content of A
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      double dx = (double)control.target[i].x - (double)control.target[j].x, dy = (double)control.target[i].y - (double)control.target[j].y;
      //printf("control.target[%d].x=%d\r\n", i, control.target[i].x);
      double d = sqrt(dx * dx + dy * dy);
      if(std::abs(d)<1e-40){
        A[i][j] = 0.0;
      }
      else{
        A[i][j] = d*d*log(d);
      }
      // printf("%lf,%lf\r\n",log(d),A[i][j]);
    }
  }
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
  // printf("A\r\n");
  // A.print();
  // printf("Ainverse\r\n");
  // Ainverse.print();
  inverseSquareM<double>(A,Ainverse);
  // A.print();
  // printf("Ainverse\r\n");
  // Ainverse.print();
  // Array2D<double> tpppl(N+3, N+3);
  // mul_MM(A, Ainverse, tpppl);
  // A.print();
  // Ainverse.print();
  // printf("mul\r\n");
  // tpppl.print();

  Array2D<double> wx(N+3,1), wy(N+3,1);

  //calculate wx
  //fill content of b
  for(int i=0;i<N;i++){
    b[i][0] = (double)control.source[i].x - (double)control.target[i].x;
  }
  b[N][0] = b[N+1][0] = b[N+2][0] = 0.0;
  // printf("\r\nbx\r\n");
  // b.print();

  mul_MM<double>(Ainverse, b, wx);

  // Array2D<double> check(N+3,1);
  // mul_MM<double>(A,wx,check);
  // printf("bx=\r\n");
  // b.print();
  // printf("\r\ncheck=\r\n");
  // check.print();

  //calculate wy
  //fill content of b
  for(int i=0;i<N;i++){
    b[i][0] = (double)control.source[i].y - (double)control.target[i].y;
  }
  b[N][0] = b[N+1][0] = b[N+2][0] = 0.0;
  // printf("\r\nbr\n");
  // b.print();

  mul_MM<double>(Ainverse, b, wy);


  // printf("\r\nwx\r\n");
  // wx.print();
  // printf("\r\nwy\r\n");
  // wy.print();

  // for(int p=0;p<control.target.size();p++){
  //
  //   Array2D<double> X(1,N+3),tmpResult(1,1);
  //   double axis[2];
  //
  //   // fill content of X
  //   for(int j=0;j<N;j++){
  //     double dx = (double)control.target[p].x - (double)control.target[j].x, dy = (double)control.target[p].y - (double)control.target[j].y;
  //     double d = sqrt(dx * dx + dy * dy);
  //     if(std::abs(d)<1e-40){
  //       X[0][j] = 0.0;
  //     }
  //     else{
  //       X[0][j] = d*d*log(d);
  //     }
  //     // printf("%lf,",X[0][j]);
  //   }
  //   // printf("\r\n");
  //   X[0][N] = 1.0;
  //   X[0][N+1] = (double)control.target[p].x;
  //   X[0][N+2] = (double)control.target[p].y;
  //
  //   // X.print();
  //
  //   //calculate dx
  //   mul_MM<double>(X,wx,tmpResult);
  //   axis[0] = control.target[p].x + tmpResult[0][0];
  //   // printf("dx=%lf,",tmpResult[0][0]);
  //   //calculate dy
  //   mul_MM<double>(X,wy,tmpResult);
  //   axis[1] = control.target[p].y + tmpResult[0][0];
  //   // printf("dy=%lf,",tmpResult[0][0]);
  //   printf("(%lf,%lf),(%d,%d)\r\n",axis[0],axis[1],control.source[p].x,control.source[p].y);
  //   // printf("rawDx=%lf,rawDy=%lf\r\n",(double)control.source[p].x - (double)control.target[p].x, (double)control.source[p].y - (double)control.target[p].y);
  // }

  // iteration over all pixels
  for(int p=0;p<rawI.size[0];p++){
    for(int q=0;q<rawI.size[1];q++){
      Array2D<double> X(1,N+3),tmpResult(1,1);
      double axis[2];

      // fill content of X
      for(int j=0;j<N;j++){
        double dx = (double)p - (double)control.target[j].x, dy = (double)q - (double)control.target[j].y;
        double d = sqrt(dx * dx + dy * dy);
        if(std::abs(d)<1e-40){
          X[0][j] = 0.0;
        }
        else{
          X[0][j] = d*d*log(d);
        }
        // printf("%lf,",X[0][j]);
      }
      // printf("\r\n");
      X[0][N] = 1.0;
      X[0][N+1] = (double)p;
      X[0][N+2] = (double)q;

      // X.print();

      //calculate dx
      mul_MM<double>(X,wx,tmpResult);
      axis[0] = p + tmpResult[0][0];
      // printf("dx=%lf,",tmpResult[0][0]);
      //calculate dy
      mul_MM<double>(X,wy,tmpResult);
      axis[1] = q + tmpResult[0][0];
      // printf("dy=%lf\r\n",tmpResult[0][0]);

      // if(p==control.target[0].x && q==control.target[0].y){
      //   printf("(%d,%d),(%lf,%lf)\r\n",control.source[0].x,control.source[0].y,axis[0],axis[1]);
      //
      //   continue;
      // }
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
    result.data[c] = ((1-u) * rawI.data[i][j].data[c] + u * rawI.data[i+1][j].data[c])*(1-v)+((1-u) * rawI.data[i][j+1].data[c] + u * rawI.data[i+1][j+1].data[c])*v;
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
