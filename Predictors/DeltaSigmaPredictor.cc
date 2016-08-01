/************************************************************************
 *                                                                      *
 * This class calculates and solves the linear prediction equations     *
 * for a bandpass delta sigma converter.                                *
 *                                                                      *
 * File:DeltaSigmaPredictor.h                                           *
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
 * The autocorrelation vector is given by:                              *
 * rxx(n) = exp(-j*IF*n)*sin(BW*n)/(2*PI*n).                            *
 * Where the IF and BW frequencies are given in rad/s.  See the January *
 * 2000 issue of the IEEE Signal Processing Magazine for more details.  *
 *                                                                      *
 ************************************************************************/

#include <math.h>
#include <stdio.h>
#include <C_Libraries/constants.h>
#include "DeltaSigmaPredictor.h"
#include "LinearPredictor.h"
#include <GNU/Complex.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))


// ############################# Private Function ###############################
// calculateAutocorrelationVector: Calculate the autocorrelation vector for a band
//                                 pass predictor.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// 1. The autocorrelation vector is stored in _autocorrelationVector.
// 2. _numberCoeffs + 1 coefficients are calculated for 1 step predictor
// ############################# Class Constructor ###############################
void DeltaSigmaPredictor::calculateAutocorrelationVector()
{
  int           n;
  Complex       carg;
  
  if(_oldNumber < _numberCoeffs)
  {
    delete [] _autocorrelationVector;
    _oldNumber                  = _numberCoeffs;
    _autocorrelationVector      = new Complex[_numberCoeffs+1];
  }
//
// The autocorrelation vector is obtained by inverse DTFT'ing a bandlimited
// white noise process centered at the IF frequency:
//
  _autocorrelationVector[0]     = Complex(_bwFrequency/TWOPI, 0.);      // sinc function at 0
  for(n=1; n<=_numberCoeffs; n++)
  {
    carg                        = Complex(0., _ifFrequency*n);
    _autocorrelationVector[n]   = exp(carg)*sin(_bwFrequency*n/2.)/TWOPI;
  }
  return;
}
    
// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the DeltaSigmaPredictor class.
//
// Input:       number:         number of filter coefficients/size of autocorrelations
//              ifFreq:         center of passband in radians [0, PI]
//              bwFreq:         passband width in radians [0, PI]
//              lambda:         matrix conditioning, i.e., A + lambda*I
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
DeltaSigmaPredictor::DeltaSigmaPredictor(int number, double ifFreq, double bwFreq, double lambda)
{
  _linearPredictor              = new LinearPredictor(number, lambda);
  _oldNumber                    = 0;
  _autocorrelationVector        = NULL;
  setNumberCoeffs(number);
  setLambda(lambda);
  setIFFrequency(ifFreq);
  setBWFrequency(bwFreq);

  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the DeltaSigmaPredictor class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
DeltaSigmaPredictor::~DeltaSigmaPredictor()
{
  delete [] _autocorrelationVector;
  delete _linearPredictor;

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
void DeltaSigmaPredictor::setNumberCoeffs(int number)
{
  _numberCoeffs = MAX(1, number);
  _numberCoeffs = MIN(_numberCoeffs, MAX_PRED_COEFFS);
  _linearPredictor->setNumberCoeffs(_numberCoeffs);
  return;
}

// ############################# Public Method ###############################
// setIFFrequency -- Sets new center frequency in radians
// Input:   number:             new center frequency [0, PI]
//          
// Output:                      None
//
// Notes:
// ############################# Public Method ###############################
void DeltaSigmaPredictor::setIFFrequency(double ifFreq)
{
  _ifFrequency  = MAX(0., ifFreq);
  _ifFrequency  = MIN(_ifFrequency, PI);
  return;
}

// ############################# Public Method ###############################
// setBWFrequency -- Sets new bandwidth in radians
// Input:   bwFreq:             bandwidth in radians [0, PI]
//          
// Output:                      None
//
// Notes:
// ############################# Public Method ###############################
void DeltaSigmaPredictor::setBWFrequency(double bwFreq)
{
  _bwFrequency  = MAX(0., bwFreq);
  _bwFrequency  = MIN(_bwFrequency, PI);
  return;
}

// ############################# Public Method ###############################
// setLambda -- Sets matrix conditioning number, lambda
// Input:   bwFreq:             bandwidth in radians [0, PI]
//          
// Output:                      None
//
// Notes:
// ############################# Public Method ###############################
void DeltaSigmaPredictor::setLambda(double lambda)
{
  _lambda       = MAX(0.,lambda);
  _linearPredictor->setLambda(_lambda);
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
const Complex *DeltaSigmaPredictor::calculatePredictorCoefficients()
{
//
// First calculate the autocorrelation vector:
//
  calculateAutocorrelationVector();
  _linearPredictor->setToeplitzRow(_autocorrelationVector);
  _linearPredictor->setAutocorrelationVector(&_autocorrelationVector[1]);       // 1 step predictor
//
// Calculate and return coefficients:
//  
  return _linearPredictor->calculatePredictorCoefficients();
}
