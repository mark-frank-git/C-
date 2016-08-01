/************************************************************************
 *                                                                      *
 * This subclass of AbstractFIR adds functionality for calculating      *
 * the filter coefficients for implmenting delay FIR filters.           *
 * The code is based on the article, "Splitting the unit delay,"        *
 * T. Laakso, V. Valimaki, M. Karjalainen, and U. Laine, IEEE Sig Proc  *
 * Magazine, Jan. 1996.                                                 *
 *                                                                      *
 * File:DelayFIR.h                                                      *
 *                                                                      *
 * The filter is stored in the forms:                                   *
 *                                                                      *
 *    H(z) =  b[n] + b[n-1]z + ... + b[0]z^n                            *
 *         =  b[0] + b[1]z^(-1) + ... + b[n]z^-(n)                      *
 *                                                                      *
 * Where b polynomials are stored in the polynomial objects.            *
 * NOTE: The order of the filter is the # taps - 1.                     *
 *                                                                      *
 ************************************************************************/

#include "DelayFIR.h"                                   // Object prototypes
#include <Polynomials/DoublePoly.h>
#include <math.h>
#include <stdio.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))

// ############################# Private Method ###############################
// transferForLagrange -- This routine finds the coefficients of the filter
//                        using the Lagrange interpolation method from p. 41 of
//                        the reference.
//
// Input:                   None
//          
// Output:                  Return YES if OK, return NO if polynomials don't exist
//
// Notes:
//  1. The instanceVariable, _bPolyObject is modified.
//  2. The number of taps is equal to order + 1.
// ############################# Private Method ###############################
LOGICAL DelayFIR::transferForLagrange()
{
  int       n, k;
  double    *b, total_delay;
  
  if(_bPolyObject == NULL)
     return NO;

  b             = new double[_numberTaps];
//
// First, calculate total delay per eqn. (21) with
// M = 0.
//
  total_delay   = findTotalDelay();
//
// Calculate the coefficients
//
  for(n=0; n<_numberTaps; n++)
  {
    b[n]        = 1.;
    for(k=0; k<_numberTaps; k++)
    {
      if(k!=n)
      {
        b[n]    *= (total_delay-k)/(n-k);
      }
    }
  }
  _bPolyObject->assign(_filterOrder, b);
  delete [] b;
//
// Set unity pass band gain:
//
  setPassBandGain();

  return YES;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the DelayFIR class.
//
// Input:       taps:                   Number of filter coefficients.
//              samplingFreq:           Digital sampling frequency in Hertz
//              delay:                  normalized delay in units of sampling time [0, 0.5]
//              type:                   Delay type (Lagrange, etc.)
//
// Output:                              None
//
// Notes:
// ############################# Class Constructor ###############################
DelayFIR::DelayFIR(int taps, float samplingFreq, float delay)
              :AbstractFIR(LOW_PASS, 0., 100., samplingFreq, taps-1)
{
  setDelay(delay);
  
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the DelayFIR class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
DelayFIR::~DelayFIR()
{
  return;
}

// ############################# Public Method ###############################
// setDelay -- Sets the normalized filter delay from 0 to 0.5 * sampling time
// Input:       delay:          new normalized delay
//          
// Output:                      None
// Notes:
// ############################# Public Method ###############################
void DelayFIR::setDelay(float delay)
{
  _normalizedDelay      = MAX(MIN_DELAY, delay);
  _normalizedDelay      = MIN(MAX_DELAY, _normalizedDelay);
  return;
}

// ############################# Public Method ###############################
// findTotalDelay -- Find the filter's total delay based on Eqn. (21) in the
//                   reference
//
// Input:               None
//          
// Output:              Integer plus fractional delay
//
// Notes:
// ############################# Public Method ###############################
float DelayFIR::findTotalDelay()
{
  float total_delay;
//
// First, calculate total delay per eqn. (21) with
// M = 0.
//
  if(_filterOrder%2==0)                 // even order
  {
    if(_normalizedDelay >= 0.)          // positive delay
    {
      if(_normalizedDelay >= 0.5)
        total_delay     = _filterOrder/2 - 1 + _normalizedDelay;
      else
        total_delay     = _filterOrder/2 + _normalizedDelay;
    }
    else                                // negative delay
    {
      if(_normalizedDelay <= -0.5)
        total_delay     = _filterOrder/2 + 1 + _normalizedDelay;
      else
        total_delay     = _filterOrder/2 + _normalizedDelay;
    }
  }
  else                                  // odd order
  {
    if(_normalizedDelay >= 0.)          // positive delay
      total_delay       = (_filterOrder-1)/2. + _normalizedDelay;
    else                                // negative delay
      total_delay       = (_filterOrder-1)/2. + 1. + _normalizedDelay;
  }
  return        total_delay;
}

// ############################# Public Method ###############################
// findTransferFunction -- Find the filter's transfer function
//
// Input:                   None
//          
// Output:                  None
//
// Notes:
//  1. The instanceVariables, _bPolyObject are modified.
// ############################# Public Method ###############################
void DelayFIR::findTransferFunction()
{
  transferForLagrange();
  zeroOutTaps();
  return;
}
