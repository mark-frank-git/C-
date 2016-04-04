#ifndef _MATRIXROUTINES_H
#define _MATRIXROUTINES_H
/************************************************************************
 *                                                                      *
 * This file contains some general purpose matrix routines.             *
 *                                                                      *
 * File:matrix_routines.h                                               *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 07/06/04  - Extracted from DoubleMatrix.                         *
 *  2. 05/17/05  - Add Complex versions and svd() routines.             *
 ************************************************************************/

#ifndef LOGICAL
#define LOGICAL                 char
#endif

//
// Double Matrix Routines
//
LOGICAL lu_decomp(const double *inputMatrix, int size, double *luMatrix, int *pivot);
LOGICAL linear_solve(const double *luMatrix, const int *pivot, double *bVector, int size);
LOGICAL invert_using_LU(const double *luMatrix, const int *pivot, int size, double *inverseMatrix);
void    toeplitz_solve(double *r, double *b, int size, double *x, double *conda, int *flag);
int     svd(const double *a, int aRows, int aCols, double *singularValues, int uFlag, double *u, int vFlag, double *v);
int     diagon(double *u, double *v, double *rv1, double *w, double anorm, int coln, int colu, int n, int m);
double  dsign(double to, double from);
//
// Complex Matrix Routines
//
LOGICAL lu_decomp(const Complex *inputMatrix, int size, Complex *luMatrix, int *pivot);
LOGICAL linear_solve(const Complex *luMatrix, const int *pivot, Complex *bVector, int size);
LOGICAL invert_using_LU(const Complex *luMatrix, const int *pivot, int size, Complex *inverseMatrix);
void    toeplitz_solve(Complex *r, Complex *b, int size, Complex *x, double *conda, int *flag);
int     svd(const Complex *aMatrix, int aRows, int aCols, double *singularValues, Complex *uMatrix, Complex *qMatrix);
double gsl_hypot (const double x, const double y);
Complex complex_dot_columns(const Complex *aMatrix, int aRows, int aCols, int column1, int column2);
double complex_norm_column(const Complex *aMatrix, int aRows, int aCols, int column);

#endif
