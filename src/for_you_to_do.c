#include "../include/for_you_to_do.h"

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 128;
  
}

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/


int mydgetrf(double *A, int *ipiv, int n) 
{
    /* add your code here */
     int i, j, k, t, maxind;
  double max;

  for (i = 0; i < (n - 1); i++) {

    maxind = i;
    max = fabs(A[i*n + i]);

    for (t = i + 1; t < n; t++) {
      if (fabs(A[t*n + i]) > max) {
        maxind = t;
        max = fabs(A[t*n + i]);
      }
    }

    if (max == 0) {
      printf("LUfactoration failed: coefficient matrix is singular\n");
      return -1;
    }
    else {
      if (maxind != i) {
        int temps = ipiv[i];
	ipiv[i] = ipiv[maxind];
        ipiv[maxind] = temps;

        for (j = 0; j < n; j++) {
          double tempv = A[i*n + j];
          A[i*n + j] = A[maxind*n + j];
          A[maxind*n +j] = tempv;
        }
      }
    }
    for (j = i + 1; j < n; j++) {
      A[j*n + i] = A[j*n + i] / A[i*n + i];
      for (k = i + 1; k < n; k++)
        A[j*n + k] = A[j*n + k] - A[j*n + i] * A[i*n + k];
    }
  }
  return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    /* add your code here */
     int i, j;
  double *y = (double*)calloc(n, sizeof(double));
  double *x = (double*)calloc(n, sizeof(double));

  if (UPLO == 'L') {
    y[0] = B[ipiv[0]];
    for (i = 1; i < n; i++) {

      double sum = 0.0;

      for (j = 0; j < (i); j++) {
        sum += y[j] * A[i*n + j];
      }
      y[i] = B[ipiv[i]] - sum;
    }
    memcpy(B, y, n * sizeof(double));
  }
  else if (UPLO == 'U') {
    memcpy(y, B, n * sizeof(double));
    x[n - 1] = y[n - 1]/A[(n - 1)*n + n - 1];
    for (i = n - 2; i >= 0; i--) {
      double sum = 0.0;

      for (j = i + 1; j < n; j++) {
        sum += x[j] * A[i*n + j];
      }

      x[i] = (y[i] - sum) / A[i*n + i];
    }
    memcpy(B, x, n * sizeof(double));
  }
  free(y);
  free(x);
  return;
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
    /* add your code here */
    /* please just copy from your lab1 function optimal( ... ) */
     int i1, j1, k1;

  for (i1 = i; i1 < i + b && i1 < n; i1 += 3) {
    for (j1 = j; j1 < j + b && j1 < n; j1 += 3) {

      register int t = i1*n + j1;
      register int tt = t + n;
      register int ttt = tt + n;
      register double c00 = C[t];
      register double c01 = C[t + 1];
      register double c02 = C[t + 2];
      register double c10 = C[tt];
      register double c11 = C[tt + 1];
      register double c12 = C[tt + 2];
      register double c20 = C[ttt];
      register double c21 = C[ttt + 1];
      register double c22 = C[ttt + 2];
      for (k1 = k; k1 < k + b && k1 < n; k1 += 3) {

	register int ta = i1*n + k1;
	register int tta = ta + n;
	register int ttta = tta + n;
	register int tb = k1*n + j1;
	register int ttb = tb + n;
	register int tttb = ttb + n;
	register double a0 = A[ta];
	register double a1 = A[tta];
	register double a2 = A[ttta];
	register double b0 = B[tb];
	register double b1 = B[tb+1];
	register double b2 = B[tb+2];

	c00 -= a0 * b0;
	c01 -= a0 * b1;
	c02 -= a0 * b2;
	c10 -= a1 * b0;
	c11 -= a1 * b1;
	c12 -= a1 * b2;
	c20 -= a2 * b0;
	c21 -= a2 * b1;
	c22 -= a2 * b2;
	a0 = A[ta + 1];
	a1 = A[tta + 1];
	a2 = A[ttta + 1];
	b0 = B[ttb];
	b1 = B[ttb + 1];
	b2 = B[ttb + 2];
	c00 -= a0 * b0;
	c01 -= a0 * b1;
	c02 -= a0 * b2;
	c10 -= a1 * b0;
	c11 -= a1 * b1;
	c12 -= a1 * b2;
	c20 -= a2 * b0;
	c21 -= a2 * b1;
	c22 -= a2 * b2;
	a0 = A[ta + 2];
	a1 = A[tta + 2];
	a2 = A[ttta + 2];
	b0 = B[tttb];
	b1 = B[tttb + 1];
	b2 = B[tttb + 2];
	c00 -= a0 * b0;
	c01 -= a0 * b1;
	c02 -= a0 * b2;
	c10 -= a1 * b0;
	c11 -= a1 * b1;
	c12 -= a1 * b2;
	c20 -= a2 * b0;
	c21 -= a2 * b1;
	c22 -= a2 * b2;
      }
      C[t] = c00;
      C[t + 1] = c01;
      C[t + 2] = c02;
      C[tt] = c10;
      C[tt + 1] = c11;
      C[tt + 2] = c12;
      C[ttt] = c20;
      C[ttt + 1] = c21;
      C[ttt + 2] = c22;
    }
  }
    return;
}

/**
 * 
 * this function computes LU decomposition
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *    
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *     
 *      b                , block size   
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
     int i, j, k, t, maxind, ib;
  double max;

  for (ib = 0; ib < n; ib += b) {

    for (i = ib; i < ib + b && i < n; i++) {

      maxind = i;
      max = fabs(A[i*n + i]);
    
      for (t = i + 1; t < n; t++) {
	if (fabs(A[t*n + i]) > max) {
	  maxind = t;
	  max = fabs(A[t*n + i]);
	}
      }
      
      if (max == 0) {
	printf("LUfactoration failed: coefficient matrix is singular\n");
	return -1;
      }
      else {
	if (maxind != i) {
	  int temps = ipiv[i];
	  ipiv[i] = ipiv[maxind];
	  ipiv[maxind] = temps;

	  for (j = 0; j < n; j++) {
	    double tempv = A[i*n + j];
	    A[i*n + j] = A[maxind*n + j];
	    A[maxind*n +j] = tempv;
	  }
	}
      }

      for (j = i + 1; j < n; j++) {
	A[j*n + i] = A[j*n + i] / A[i*n + i];
	for (k = i + 1; k < ib + b && k < n; k++)
	  A[j*n + k] = A[j*n + k] - A[j*n + i] * A[i*n + k];
      }
    }
    for (i = ib; i < ib + b && i < n; i++) {
      for (j = ib + b; j < n; j++) {
	double sum = 0.0;
	for (k = ib; k < i; k++) {
	  sum += A[i*n + k] * A[k*n + j];
	}
	A[i*n + j] -= sum;
      }
    }
    for (i = ib + b; i < n; i += b) {
      for (j = ib + b; j < n; j += b) {
	mydgemm(A, A, A, n, i, j, ib, b);
      }
    }
  }
  return 0;
}

