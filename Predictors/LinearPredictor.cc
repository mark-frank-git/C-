/************************************************************************
 *                                                                      *
 * This class solves the linear prediction equations for the linear     *
 * predictor filter coefficients. The Levinson Durbin algorithm is an   *
 * efficient means of solving the equations.                            *
 *                                                                      *
 * File:LinearPredictor.h                                               *
 *                                                                      *
 * The autocorrelation or linear prediction equations take the form:    *
 *                                                                      *
 *    Rxx * A = R                                                       *
 * where                                                                *
 * Rxx = |rxx(0) rxx(1) .. rxx(p-1)|                                    *
 *       |rxx(1) rxx(0) .. rxx(p-2)|                                    *
 *       |                         |                                    *
 *       |rxx(p-1)..        rxx(0) |                                    *
 * is the symmetric autocorrelation matrix.                             *
 * A = [a(0) a(1) .. a(p-1)]T   is the coefficient vector.              *
 * R = [rxx(1) rxx(2) .. rxx(p)] is the one step predictor correlation  *
 *                               vector.                                *
 * The purpose of this class is to solve for the coefficient vector, A  *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 04/04/00  - Started.                                             *
 ************************************************************************/

#include <math.h>
#include <stdio.h>
#include "LinearPredictor.h"                                    // Object prototypes
#include <GNU/Complex.h>
#include <MatrixAlgebra/ComplexMatrix.h>
#include <MatrixAlgebra/DoubleMatrix.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))


// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the LinearPredictor class.
//
// Input:       number:         number of filter coefficients/size of autocorrelations
//              lambda:         matrix conditioning, i.e., A + lambda*I
//              row:            first row of autocorrelation matrix
//              vector:         vector of autocorrelation
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
LinearPredictor::LinearPredictor(int number, double lambda, Complex *row, Complex *vector)
{
  _predictorMatrix              = NULL;
  _predictorCoefficients        = NULL;
  _oldNumber                    = 0;
  setNumberCoeffs(number);
  setLambda(lambda);
  setToeplitzRow(row);
  setAutocorrelationVector(vector);


  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the LinearPredictor class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
LinearPredictor::~LinearPredictor()
{
  delete [] _predictorCoefficients;
  delete _predictorMatrix;

  return;
}

// ############################# Public Method ###############################
// setNumberCoeffs -- Sets new number of filter coefficients to calculate
// Input:   number:             new number of coefficients
//          
// Output:                      None
//
// Notes:
// 1. The number of coefficients must match the size of the autocorrelation
//    matrix.
// ############################# Public Method ###############################
void LinearPredictor::setNumberCoeffs(int number)
{
  _numberCoeffs = MAX(1, number);
  _numberCoeffs = MIN(_numberCoeffs, MAX_PRED_COEFFS);
  return;
}

// ############################# Public Method ###############################
// calculatePredictorCoefficients -- Calculate and return the predictor
//                                   coefficients from the autocorrelation matrix
//
// Input:               None
//          
// Output:              array of predictor coefficients
//
// Notes:
//  1. This needs to be overridden in subclasses.
// ############################# Public Method ###############################
const Complex *LinearPredictor::calculatePredictorCoefficients()
{
  int           i, flag;
  double        conda;
  const Complex *output_coeffs;
  Complex       *toeplitz_plus_lambda;
//
// Error check:
//
  if((_toeplitzRow==NULL) || (_autocorrelationVector==NULL))
    return NULL;
  toeplitz_plus_lambda  = new Complex[_numberCoeffs];
  for(i=0; i<_numberCoeffs; i++)
    toeplitz_plus_lambda[i]     = _toeplitzRow[i];
  toeplitz_plus_lambda[0]       += _lambda;
//
// Allocate output array:
//
  if(_oldNumber < _numberCoeffs)
  {
    _oldNumber  = _numberCoeffs;
    delete [] _predictorCoefficients;
    _predictorCoefficients      = new Complex[_numberCoeffs];
  }
//
// Allocate Matrix object for solving equations:
//
  if(_predictorMatrix == NULL)
  {
    _predictorMatrix    = new ComplexMatrix();
  }
//
// Call the Toeplitz solver:
//
  _predictorMatrix->toeplitzSolve(toeplitz_plus_lambda, _autocorrelationVector, _numberCoeffs, &conda, &flag);
  output_coeffs = _predictorMatrix->getMatrix();
//
  for(i=0; i<_numberCoeffs; i++)
    _predictorCoefficients[i]   = output_coeffs[i];

  delete [] toeplitz_plus_lambda;
  
  return _predictorCoefficients;
}
