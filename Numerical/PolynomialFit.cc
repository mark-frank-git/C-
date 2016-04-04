/************************************************************************
 *                                                                      *
 * This subclass of object fits a polynomial to a set of data.          *
 *                                                                      *
 *                                                                      *
 * File:PolynomialFit.cc                                                *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 11/17/00 - Started.                                              *
 *  2. 07/14/04 - Used DoubleSquareMatrix.                              *
 ************************************************************************/

#include "PolynomialFit.h"
#include <stdio.h>
#include <Polynomials/DoublePoly.h>
#include <MatrixAlgebra/DoubleMatrix.h>
#include <MatrixAlgebra/DoubleSquareMatrix.h>

#define ABS(a)          ((a) >= 0 ? (a) : (-a))
#define SIGN(a,b)       ( ((b)>=0 ) ? ABS((a)) : -ABS((a)) )
#define MAX(a,b)        ((a)>(b)?(a):(b))
#define EPS             0.0000000001
// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the PolynomialFit class.
// Input:       order:          Order of polynomial to fit to data
//          
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
PolynomialFit::PolynomialFit(int order)
{
  setPolyOrder(order);
  setRegularization(0.);
  _fitPoly      = new DoublePoly(_polyOrder);
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the PolynomialFit class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
PolynomialFit::~PolynomialFit()
{
  delete _fitPoly;
  return;
}

// ############################# Public Method ###############################
// fitToData -- Fit a polynomial to the data
// Input:       xData:          Array of x data points
//              yData:          Array of y data points
//              numberPoints:   Size of arrays
//          
// Output:                      Returns YES if fit went OK, NO otherwise.
//
// Notes:
// 1. Call fitPoly() to get polynomial
// ############################# Public Method ###############################
LOGICAL PolynomialFit::fitToData(const float *xData, const float *yData, int numberPoints)
{
  int                   i, j, k, rows, columns;
  LOGICAL               flag;
  DoubleMatrix          A, b;
  DoubleMatrix          C, D;
  DoubleSquareMatrix    A_transp_A;
  DoubleSquareMatrix    *regularizer;
  double                product;
  double                *a_array, *b_array;
  const double          *x_array;

  if( (xData==NULL) || (yData==NULL) || (numberPoints<(_polyOrder+1)))
    return NO;
//
// Set up matrix sizes:
//
  rows          = numberPoints;
  columns       = _polyOrder + 1;
  a_array       = new double[rows*columns];
  b_array       = new double[rows];
//
// Calculate the elements of the A matrix:
//
  k                     = 0;
  for(i=0; i<rows; i++)
  {
    product             = 1.;
    for(j=0; j<columns; j++)
    {
      a_array[k++]      = product;
      product           *= xData[i];
    }
    b_array[i]          = yData[i];
  }
//
// Initialize the A and b matrices
//
  A.assign(rows, columns, a_array);
  b.assign(rows, 1, b_array);
//
// Get two copies of  A transpose:
//
  C     = A;
  C.transpose();
  D     = C;
//
// Get A'A, add lambda*I
//
  A_transp_A            = D*A;
  if(_lambda > EPS)
  {
    regularizer         = new DoubleSquareMatrix(_lambda, columns);
    A_transp_A          += *regularizer;
    delete regularizer;
  }
  _matrixDeterminant    = A_transp_A.findDeterminant();
//
// Find inverse:
//
  A_transp_A.invert(&flag);
  if(flag == NO)
    return NO;
//
// Find (A'A)A' (pseudo-inverse)
//
  D     = A_transp_A*C;         // Post multiply D by A'
  D     *= b;
//
// Get elements, and assign to output polynomial
//
  x_array       = D.getMatrix();
  _fitPoly->assign(_polyOrder, x_array);
  _fitPoly->reverse();
  
  delete [] a_array;
  delete [] b_array;

  return YES;
}
// ############################# Public Method ###############################
// fitCoefficients -- Returns the coefficients of the fit polynomial
// Input:                       None
//          
// Output:                      Polynomial coefficients
// ############################# Public Method ###############################
const double *PolynomialFit::fitCoefficients()
{
  return _fitPoly->getCoefficients();

}

// ############################# Public Method ###############################
// setPolyOrder -- Sets a new polynomial order
// Input:       order:          New polynomial order
//          
// Output:                      None
// ############################# Public Method ###############################
void PolynomialFit::setPolyOrder(int order)
{
  _polyOrder    = MAX(1, order);
  return;
}



