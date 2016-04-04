/************************************************************************
 *                                                                      *
 * This class solves for the feedforward and/or feedback coefficients   *
 * for a general delta sigma modulator.  It is based on the block       *
 * diagram in 5.6.4 in Delta-Sigma Data Converters: Theory, Design,     *
 * and Simulation.                                                      *
 *                                                                      *
 * File:FeedforwardFeedback.h                                           *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 04/23/00  - Started.                                             *
 ************************************************************************/

#include <math.h>
#include <stdio.h>
#include "FeedforwardFeedback.h"                                        // Object prototypes
#include <GNU/Complex.h>
#include <Polynomials/ComplexPoly.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))


// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the FeedforwardFeedback class.
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
FeedforwardFeedback::FeedforwardFeedback(int order, Complex *ntfZeros, Complex *stfZeros, Complex *poles)
{
  _oldFFCoefficients            = 0;
  _oldFBCoefficients            = 0;
  _ffCoefficients               = NULL;
  _fbCoefficients               = NULL;
  setFilterOrder(order);
  setNTFZeros(ntfZeros);
  setSTFZeros(stfZeros);
  setTFPoles(poles);


  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the FeedforwardFeedback class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
FeedforwardFeedback::~FeedforwardFeedback()
{
  delete [] _ffCoefficients;
  delete [] _fbCoefficients;

  return;
}

// ############################# Public Method ###############################
// setFilterOrder -- Sets new delta sigma filter order
// Input:   order:              new filter order
//          
// Output:                      None
//
// Notes:
// 1. The number of NTF zeros, STF zeros and TF poles must match the filter order.
// ############################# Public Method ###############################
void FeedforwardFeedback::setFilterOrder(int order)
{
  _filterOrder  = MAX(1, order);
  _filterOrder  = MIN(_filterOrder, MAX_DELTA_ORDER);
  return;
}

// ############################# Public Method ###############################
// calculateFFCoefficients -- Calculate and return the delta sigma feedforward
//                            coefficients given the desired poles and zeros
//                            of the noise and signal transfer functions.
//
// Input:               None
//          
// Output:              array of feedforward coefficients
//
// Notes:
// ############################# Public Method ###############################
const Complex *FeedforwardFeedback::calculateFFCoefficients()
{
//
// Check for error:
//
  if(_ntfZeros==NULL || _stfZeros==NULL)
    return NULL;
//
// Allocate output array
//
  if(_oldFFCoefficients < _filterOrder)
  {
    delete [] _ffCoefficients;
    delete [] _fbCoefficients;
    _ffCoefficients     = new Complex[_filterOrder+1];
    _fbCoefficients     = new Complex[_filterOrder+1];
    _oldFFCoefficients  = _filterOrder;
  }
//
// Calculate the coefficients based on the filter order, and poles and zeros:
//
  switch(_filterOrder)
  {
    case DELTA_FIRST_ORDER:
      _ffCoefficients[0]        = Complex(1.,0.);               // b1
      break;
    case DELTA_SECOND_ORDER:
      _ffCoefficients[0]        = _ntfZeros[0] - _stfZeros[0];  // b1
      _ffCoefficients[1]        = Complex(1., 0.);              // b2
      break;
    case DELTA_THIRD_ORDER:
      _ffCoefficients[2]        = Complex(1., 0.);              // b3
      _ffCoefficients[1]        = _ntfZeros[0]+_ntfZeros[1]-_stfZeros[0]-_stfZeros[1];
                                                                // b2
      _ffCoefficients[0]        = _stfZeros[0]*_stfZeros[1]+_ntfZeros[0]*(_ffCoefficients[1]-_ntfZeros[1]);
      break;
    default:
      break;
  }
  return _ffCoefficients;
}

// ############################# Public Method ###############################
// calculateFBCoefficients -- Calculate and return the delta sigma feedback
//                            coefficients given the desired poles and zeros
//                            of the noise and signal transfer functions.
//
// Input:               None
//          
// Output:              array of feedback coefficients
//
// Notes:
// ############################# Public Method ###############################
const Complex *FeedforwardFeedback::calculateFBCoefficients()
{
  ComplexPoly   denominator_poly;
  const Complex *denominator_coeff;
//
// Check for error:
//
  if(_ntfZeros==NULL || _stfZeros==NULL || _tfPoles==NULL)
    return NULL;
//
// Allocate output array
//
  if(_oldFBCoefficients < _filterOrder)
  {
    delete [] _fbCoefficients;
    _fbCoefficients     = new Complex[_filterOrder+1];
    _oldFBCoefficients  = _filterOrder;
  }
  switch(_filterOrder)
  {
    case DELTA_FIRST_ORDER:
      _fbCoefficients[0]        = _ntfZeros[0]-_tfPoles[0];     // a1
      break;
    case DELTA_SECOND_ORDER:
      denominator_poly.setRoots(2, _tfPoles);                   // Find coefficients of poles
      printf("poles poly = %s\n", denominator_poly.getPolyString());
      denominator_coeff         = denominator_poly.getCoefficients();
      _fbCoefficients[1]        = denominator_coeff[1]+_ntfZeros[0]+_ntfZeros[1];       // a2
      _fbCoefficients[0]        = denominator_coeff[2]+
                                  _ntfZeros[0]*(_fbCoefficients[1]-_ntfZeros[1]);       // a1
      break;
    case DELTA_THIRD_ORDER:
      denominator_poly.setRoots(3, _tfPoles);                   // Find coefficients of poles
      printf("poles poly = %s\n", denominator_poly.getPolyString());
      denominator_coeff         = denominator_poly.getCoefficients();
      _fbCoefficients[2]        = denominator_coeff[1]+_ntfZeros[0]+_ntfZeros[1]+_ntfZeros[2];
                                                                                        // a3
      _fbCoefficients[1]        = denominator_coeff[2]+_fbCoefficients[2]*(_ntfZeros[0]+_ntfZeros[1])
                                  -_ntfZeros[0]*_ntfZeros[1]-_ntfZeros[1]*_ntfZeros[2]-_ntfZeros[0]*_ntfZeros[2];
                                                                                        // a2
      _fbCoefficients[0]        = denominator_coeff[3]+_ntfZeros[0]*(_fbCoefficients[1]-
                                  _fbCoefficients[2]*_ntfZeros[1]+_ntfZeros[1]*_ntfZeros[2]);
      break;
  }
  return _fbCoefficients;
}
