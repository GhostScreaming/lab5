// See LICENSE for license details.

#include "stdlib.h"
#include "dataset.h"

//--------------------------------------------------------------------------
// vvadd function

void __attribute__((noinline)) vvadd(int coreid, int ncores, size_t n, const data_t* x, const data_t* y, data_t* z)
{
   size_t i;

   // interleave accesses
   for (i = coreid; i < n; i+=ncores)
   {
      z[i] = x[i] + y[i];
   }
}

void __attribute__((noinline)) vvadd_opt1(int coreid, int ncores, size_t n, const data_t* x, const data_t* y, data_t* z)
{
  // TODO: Your code here
   size_t i, l = coreid*n/ncores, r = (coreid+1)*n/ncores;

   // interleave accesses
   for (i = l; i < r; i++)
   {
      z[i] = x[i] + y[i];/*
      z[i+1] = x[i+1] + y[i+1];
      z[i+2] = x[i+2] + y[i+2];
      z[i+3] = x[i+3] + y[i+3];*/
   }
}

void __attribute__((noinline)) vvadd_opt2(int coreid, int ncores, size_t n, const data_t* x, const data_t* y, data_t* z)
{
  // TODO: Your code here
   size_t i, l = coreid*n/ncores, r = (coreid+1)*n/ncores;

   // interleave accesses
   for (i = l; i < r; i+=4)
   {
      z[i] = x[i] + y[i];
      z[i+1] = x[i+1] + y[i+1];
      z[i+2] = x[i+2] + y[i+2];
      z[i+3] = x[i+3] + y[i+3];
   }
}
