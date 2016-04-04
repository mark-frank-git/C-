#ifndef _LINEARPREDICTOR_H
#define _LINEARPREDICTOR_H    1
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

#ifndef YES
#define YES     1
#define NO      0
#endif

#define MAX_PRED_COEFFS 256

class Complex;                                  // class prototype
class ComplexMatrix;

class LinearPredictor
{
protected:
  int           _numberCoeffs;                  // Number of predictor coefficients
  int           _oldNumber;                     // Old number of coefficients

  double        _lambda;                        // constant added to diagonal of autocorrelation matrix

  Complex       *_toeplitzRow;                  // First row of autocorrelation matrix
  Complex       *_autocorrelationVector;        // Autocorrelation vector
  Complex       *_predictorCoefficients;        // Calculated predictor coefficients

  ComplexMatrix *_predictorMatrix;              // Matrix object used to calculate coefficients

//
// Private functions:
//

// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  LinearPredictor(int number=0, double lambda=0.0, Complex *row=NULL, Complex *vector=NULL);
  ~LinearPredictor();

/**********************
 * Set parameters:    *
 **********************/
  void  setNumberCoeffs(int number);
  void  setLambda(float lambda)                         {_lambda        = lambda;               return;}
  void  setToeplitzRow(Complex *row)                    {_toeplitzRow   = row;                  return;}
  void  setAutocorrelationVector(Complex *vector)       {_autocorrelationVector = vector;       return;}

/**********************
 * Get parameters:    *
 **********************/
  int           numberCoeffs()          {return _numberCoeffs;}
  double        lambda()                {return _lambda;}
  Complex       *toeplitzRow()          {return _toeplitzRow;}
  Complex       *autocorrelationVector(){return _autocorrelationVector;}

/**********************
 * Getting Outputs:   *
 **********************/
  const Complex *calculatePredictorCoefficients();

};
#endif
