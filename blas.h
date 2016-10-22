#ifndef _BLAS_H
#define _BLAS_H

#include "stddef.h"
#include "stdint.h"
#include "assert.h"
#include "imageprocess.h"

template <class T, int M, int N, int O>
inline void mul_MM(T A[M][N], T B[N][O], T result[M][O])
{
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

template <class T>
inline void inverseSquareM(const 2DArray<T> &A, 2DArray<T> &result)
{
  assert(A.size[0] == A.size[1] && result.size[0] == result.size[1] && A.szie[0] == result.size[0]);
  size_t N = A.size[0];
  2DArray<T> tmpMatrix(N,2*N);
  // memcpy(tmpMatrix, A, N*N*sizeof(T));
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      tmpMatrix[i][j] = A[i][j];
    }
  }

  for(int i=0;i<N;i++)
   {
      for(int j=N;j<2*N;j++)
      {
          if(i==j-N)
             tmpMatrix[i][j]=1;
         else
             tmpMatrix[i][j]=0;
       }
   }
   for(int i=0;i<N;i++)
   {
      T t = tmpMatrix[i][i];
      for(int j=i;j<2*N;j++){
        tmpMatrix[i][j] = tmpMatrix[i][j]/t;
      }
      for(int j=0;j<N;j++)
      {
         if(i!=j)
         {
            t=tmpMatrix[j][i];
            for(int k=0;k<2*N;k++){
                tmpMatrix[j][k]=tmpMatrix[j][k]-t*tmpMatrix[i][k];
            }
          }
      }
   }
   for(int i=0;i<N;i++){
     for(int j=0;j<N;j++){
       result[i][j] = tmpMatrix[i][j+N];
     }
   }
}

#define PI 3.141592653589793238463

#endif
