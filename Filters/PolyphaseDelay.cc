/************************************************************************
 *                                                                      *
 * This subclass of PolyphaseInterpolator implements a polyphase filter *
 * for delaying an input signal.                                        *
 * Based on Crochiere and Rabiner, pp. 79-85.                           *
 *                                                                      *
 * File:PolyphaseDelay.cc                                               *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 04/05/02  - Started.                                             *
 ************************************************************************/

#include "PolyphaseDelay.h"                                     // Object prototypes
#include "AbstractFIR.h"
#include <math.h>
#include <stdio.h>

#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define ABS(a)          ((a) >= 0 ? (a) : (-a))

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the PolyphaseDelay class.
//
// Input:       numberPhases:   Number of polyphase branches
//              tapsPerPhase:   Size of filter for each phase
//              cutoffFactor:   Cutoff factor for each polyphase filter for dec/interp
//              delay:          normalized delay referenced to 1 sampling period
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
PolyphaseDelay::PolyphaseDelay(int numberPhases, int tapsPerPhase, float cutoffFactor, float delay)
                        :PolyphaseInterpolator(numberPhases, tapsPerPhase, cutoffFactor)
{
  setDelay(delay);
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the PolyphaseDelay class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
PolyphaseDelay::~PolyphaseDelay()
{
 return;
}

// ############################# Public Method ###############################
// setNumberPhases -- Sets a new number of phases/number of filter branches
//
// Input:       phases:         number of phases/polyphase filters
//          
// Output:                      None
// ############################# Public Method ###############################
void PolyphaseDelay::setNumberPhases(int phases)
{
//
// Perform super's function:
//
  PolyphaseInterpolator::setNumberPhases(phases);
//
// Find the filter having the closest delay, see p. 82 in Rabiner and Schaefer
//
  setDelay(_normalizedDelay);
  return;
}

// ############################# Public Method ###############################
// setTapsPerPhase -- Sets a new number of taps per each polyphase branch
//
// Input:       taps:           number of taps for each polyphase filter
//          
// Output:                      None
// ############################# Public Method ###############################
void PolyphaseDelay::setTapsPerPhase(int taps)
{
//
// Perform super's function:
//
  PolyphaseFilter::setTapsPerPhase(taps);
//
// Add this:
//
  setDelay(_normalizedDelay);

  return;
}

// ############################# Public Method ###############################
// setDelay -- Sets the filter delay.
//
// Input:       delay:          new normalized delay
//          
// Output:                      None
// ############################# Public Method ###############################
void PolyphaseDelay::setDelay(float delay)
{
  int   i, n;
  float error;
  float min_error, polyphase_delay, fractional_delay;
//
// Error check input:
//
  _normalizedDelay      = MAX(0., delay);
  _normalizedDelay      = MIN(1., _normalizedDelay);
//
// Find the filter having the closest delay, see p. 82 in Rabiner and Schaefer
//
  min_error             = 10.;
  _selectedPolyphase    = 0;
  n                     = _tapsPerPhase*_interpolationRatio;    // prototype filter taps
  for(i=0; i<_interpolationRatio; i++)
  {
    polyphase_delay     = (n - 2*i - 1)/2./_interpolationRatio;
    fractional_delay    = polyphase_delay - (int)polyphase_delay;
    error               = (fractional_delay - _normalizedDelay);
    if(ABS(error) < min_error)
    {
      min_error                 = ABS(error);
      _selectedPolyphase        = i;
    }
  }

  return;
}

// ############################# Public Method ###############################
// filterFloatArray -- Filter the input array and overwrite it with output.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   x:              An array of floating point data, modified on output
//          numberPts:      Number of points in the input array.
//          
// Output:                  YES if went OK, No if _bPolyObject doesn't exist.
// ############################# Public Method ###############################
LOGICAL PolyphaseDelay::filterFloatArray(float * input, int numberPts)
{
//
// Filter the input data using the selected delay filter:
//
  return _polyphaseFilters[_selectedPolyphase]->filterFloatArray(input, numberPts);
}

// ############################# Public Method ###############################
// filterFloatData -- Filter an input data sample, and return the sampled value.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   input:          An input data sample
//          
// Output:                  The filtered data sample
// ############################# Public Method ###############################
double PolyphaseDelay::filterFloatData(float input)
{
//
// Filter the input data using the selected delay filter:
//
  return _polyphaseFilters[_selectedPolyphase]->filterFloatData(input);
}

// ############################# Public Method ###############################
// filterDoubleData -- Filter an input data sample, and return the sampled value.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   input:          An input data sample
//          
// Output:                  The filtered data sample
// ############################# Public Method ###############################
double PolyphaseDelay::filterDoubleData(double input)
{
//
// Filter the input data using the selected delay filter:
//
  return _polyphaseFilters[_selectedPolyphase]->filterDoubleData(input);
}

// ############################# Public Method ###############################
// filterComplexArray -- Filter the input array and overwrite it with output.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   x:              An array of floating point data, modified on output
//          numberPts:      Number of points in the input array.
//          
// Output:                  YES if went OK, No if _bPolyObject doesn't exist.
// ############################# Public Method ###############################
LOGICAL PolyphaseDelay::filterComplexArray(Complex *x, int numberPts)
{
//
// Filter the input data using the selected delay filter:
//
  return _polyphaseFilters[_selectedPolyphase]->filterComplexArray(x, numberPts);
}

// ############################# Public Method ###############################
// filterComplexData -- Filter a complex data sample, and return the sampled value.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   xin:            An input data sample
//          
// Output:                  The filtered data sample
// ############################# Public Method ###############################
Complex  PolyphaseDelay::filterComplexData(Complex &input)
{
//
// Filter the input data using the selected delay filter:
//
  return _polyphaseFilters[_selectedPolyphase]->filterComplexData(input);
}