#ifndef _BLAS_H
#define _BLAS_H

// 实现了矩阵乘法和求逆
#include "stdlib.h"
#include "stddef.h"
#include "stdint.h"
#include "assert.h"
#include "imageprocess.h"

template <class T>
inline void mul_MM(const Array2D<T> &A, const Array2D<T> &B, Array2D<T> &result)
{
  assert(A.size[1] == B.size[0] && A.size[0] == result.size[0] && B.size[1] == result.size[1]);
  size_t M = A.size[0], N = A.size[1], O = B.size[1];

  for(size_t i=0;i<M;i++){
    for(size_t j=0;j<O;j++){
      T tmpResult = 0;
      for(size_t k=0;k<N;k++){
        tmpResult += A[i][k]*B[k][j];
      }
      result[i][j] = tmpResult;
    }
  }
}

// 这个求逆函数来自网络，精度比较差
template <class T>
inline void inverseSquareM(const Array2D<T> &RAW, Array2D<T> &A)
{
  assert(RAW.size[0] == RAW.size[1] && A.size[0] == A.size[1] && RAW.size[0] == A.size[0]);
  size_t n = RAW.size[0];
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      A[i][j] = RAW[i][j];
    }
  }

  int *is,*js,i,j,k,l;
  double temp,fmax;
  is=(int *)malloc(n*sizeof(int));
  js=(int *)malloc(n*sizeof(int));

  for(k=0;k<n;k++)
   {
      fmax=0.0;
      for(i=k;i<n;i++)
       for(j=k;j<n;j++)
        { temp=fabs(A[i][j]);//找最大值
          if(temp>fmax)
           { fmax=temp;
             is[k]=i;js[k]=j;
            }
         }
      if((fmax+1.0)==1.0)
      {
        free(is);free(js);
        printf("no inv");
        return;
      }
    if((i=is[k])!=k)
      for(j=0;j<n;j++)
      {
        double tmp = A[i][j];
        A[i][j] = A[k][j];
        A[k][j] = tmp;
      }
   if((j=js[k])!=k)
     for(i=0;i<n;i++){
       double tmp = A[i][j];
       A[i][j] = A[i][k];
       A[i][k] = tmp;
     }
   A[k][k]=1.0/A[k][k];

   for(j=0;j<n;j++)
     if(j!=k)
      A[k][j]*=A[k][k];
   for(i=0;i<n;i++)
      if(i!=k)
        for(j=0;j<n;j++)
          if(j!=k)
            A[i][j]=A[i][j]-A[i][k]*A[k][j];
   for(i=0;i<n;i++)
     if(i!=k)
      A[i][k] *= -A[k][k];
  }
  for(k=n-1;k>=0;k--)
   {
     if((j=js[k])!=k)
        for(i=0;i<n;i++){
          double tmp = A[j][i];
          A[j][i] = A[k][i];
          A[k][i] = tmp;
        }
     if((i=is[k])!=k)
        for(j=0;j<n;j++)
        {
          double tmp = A[j][i];
          A[j][i] = A[j][k];
          A[j][k] = tmp;
        }
  }
  free(is);
  free(js);
}

#define PI 3.141592653589793238463

#endif
