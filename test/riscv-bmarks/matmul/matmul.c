// See LICENSE for license details.

#include "dataset.h"
#include "util.h"
#include <stddef.h>

#define STEP 2

void matmul(const size_t coreid, const size_t ncores, const size_t lda,  const data_t A[], const data_t B[], data_t C[])
{
  size_t i, j, k;

  for (i = 0; i < lda; i++) {
    for (j = coreid; j < lda; j += ncores) {
      data_t sum = 0;
      for (k = 0; k < lda; k++)
        sum += A[j*lda + k] * B[k*lda + i];
      C[i + j*lda] = sum;
    }
  }
}

void matmul_opt1(const size_t coreid, const size_t ncores, const size_t lda,  const data_t A[], const data_t B[], data_t C[])
{
  size_t i, j, k, u = coreid*lda/ncores, d = (coreid+1)*lda/ncores;

  for (i = 0; i < lda; i++) {
    for (j = u; j < d; j++) {
      data_t sum = 0;
      for (k = 0; k < lda; k++)
        sum += A[j*lda + k] * B[k*lda + i];
      C[i + j*lda] = sum;
    }
  }
}

void matmul_opt2(const size_t coreid, const size_t ncores, const size_t lda,  const data_t A[], const data_t B[], data_t C[])
{
  size_t i, j, k;

  if(coreid == 0)
  {
    for(i = 0; i < DIM_SIZE/2; i++)
      for(j = 0; j < DIM_SIZE/2; j++){
        data_t sum = 0;
        for(k = 0; k < DIM_SIZE; k++)
          sum += A[j*lda + k] * B[k*lda + i];
        C[i + j*lda] = sum;
      }
    barrier(ncores);
    for(i = DIM_SIZE/2; i < DIM_SIZE; i++)
      for(j = 0; j < DIM_SIZE/2; j++){
        data_t sum = 0;
        for(k = 0; k < DIM_SIZE; k++)
          sum += A[j*lda + k] * B[k*lda + i];
        C[i + j*lda] = sum;
      }
  }
  if(coreid == 1)
  {
    for(i = DIM_SIZE/2; i < DIM_SIZE; i++)
      for(j = DIM_SIZE/2; j < DIM_SIZE; j++){
        data_t sum = 0;
        for(k = 0; k < DIM_SIZE; k++)
          sum += A[j*lda + k] * B[k*lda + i];
        C[i + j*lda] = sum;
      }
    barrier(ncores);
    for(i = 0; i < DIM_SIZE/2; i++)
      for(j = DIM_SIZE/2; j < DIM_SIZE; j++){
        data_t sum = 0;
        for(k = 0; k < DIM_SIZE; k++)
          sum += A[j*lda + k] * B[k*lda + i];
        C[i + j*lda] = sum;
      }
  }
}

void matmul_opt3(const size_t coreid, const size_t ncores, const size_t lda,  const data_t A[], const data_t B[], data_t C[])
{
  size_t i, j, k, l, u = coreid*lda/ncores, d = (coreid+1)*lda/ncores;
  size_t temp[STEP][DIM_SIZE];

  for (i = 0; i < lda; i+=STEP) {
    for(j = 0; j < lda; j++) {
      for(k = 0; k < STEP; k++) {
        temp[k][j] = B[j*lda + i + k];
      }
    }
    for (j = u; j < d; j+=1) {
      data_t sum[STEP];
      for (k = 0; k < STEP; k++)
        sum[k] = 0;
      for (k = 0; k < lda; k++)
        for (l = 0; l < STEP; l++) {
          sum[l] += A[j*lda + k] * temp[l][k];
        }
      for (k = 0; k < STEP; k++)
        C[j*lda + i + k] = sum[k];
    }
  }
}

void matmul_opt4(const size_t coreid, const size_t ncores, const size_t lda,  const data_t A[], const data_t B[], data_t C[])
{
  size_t i, j, k, l, u = coreid*lda/ncores, d = (coreid+1)*lda/ncores;
  size_t temp[STEP][DIM_SIZE];

  
  if(coreid == 0)
  {
    for (i = 0; i < DIM_SIZE/2; i+=STEP) {
      for(j = 0; j < lda; j++) {
        for(k = 0; k < STEP; k++) {
          temp[k][j] = B[j*lda + i + k];
        }
      }
      for (j = 0; j < DIM_SIZE/2; j++) {
        data_t sum[STEP];
        for (k = 0; k < STEP; k++)
          sum[k] = 0;
        for (k = 0; k < lda; k++)
          for (l = 0; l < STEP; l++) {
            sum[l] += A[j*lda + k] * temp[l][k];
          }
        for (k = 0; k < STEP; k++)
          C[j*lda + i + k] = sum[k];
      }
    }
    barrier(ncores);
    for (i = DIM_SIZE/2; i < DIM_SIZE; i+=STEP) {
      for(j = 0; j < lda; j++) {
        for(k = 0; k < STEP; k++) {
          temp[k][j] = B[j*lda + i + k];
        }
      }
      for (j = 0; j < DIM_SIZE/2; j++) {
        data_t sum[STEP];
        for (k = 0; k < STEP; k++)
          sum[k] = 0;
        for (k = 0; k < lda; k++)
          for (l = 0; l < STEP; l++) {
            sum[l] += A[j*lda + k] * temp[l][k];
          }
        for (k = 0; k < STEP; k++)
          C[j*lda + i + k] = sum[k];
      }
    }
  }
  if(coreid==1)
  {
    for (i = DIM_SIZE/2; i < DIM_SIZE; i+=STEP) {
      for(j = 0; j < lda; j++) {
        for(k = 0; k < STEP; k++) {
          temp[k][j] = B[j*lda + i + k];
        }
      }
      for (j = DIM_SIZE/2; j < DIM_SIZE; j++) {
        data_t sum[STEP];
        for (k = 0; k < STEP; k++)
          sum[k] = 0;
        for (k = 0; k < lda; k++)
          for (l = 0; l < STEP; l++) {
            sum[l] += A[j*lda + k] * temp[l][k];
          }
        for (k = 0; k < STEP; k++)
          C[j*lda + i + k] = sum[k];
      }
    }
    barrier(ncores);
    for (i = 0; i < DIM_SIZE/2; i+=STEP) {
      for(j = 0; j < lda; j++) {
        for(k = 0; k < STEP; k++) {
          temp[k][j] = B[j*lda + i + k];
        }
      }
      for (j = DIM_SIZE/2; j < DIM_SIZE; j++) {
        data_t sum[STEP];
        for (k = 0; k < STEP; k++)
          sum[k] = 0;
        for (k = 0; k < lda; k++)
          for (l = 0; l < STEP; l++) {
            sum[l] += A[j*lda + k] * temp[l][k];
          }
        for (k = 0; k < STEP; k++)
          C[j*lda + i + k] = sum[k];
      }
    }
  }
}

