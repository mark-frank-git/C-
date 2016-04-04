/************************************************************************
 *                                                                      *
 * This class implements generates coefficients for wavelet filters     *
 *                                                                      *
 * File:WaveletFilterCoefficients.h                                     *
 *                                                                      *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 08/09/05 - Started.                                              *
 *                                                                      *
 ************************************************************************/
#include "WaveletFilterCoefficients.h"                                  // Object prototypes
#if defined(WIN32)
#include <GNU/Complex.h>
#include <Filters/DataWindow.h>
#else
#include "Complex.h"
#include "DataWindow.h"
#endif

#include <math.h>
#include <stdio.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))
const double EPS        = 1.e-20;
const double PI         = 3.141592654;

const int       DEFAULT_FILTER_TYPE     = 0;
const int       DEFAULT_WINDOW_TYPE     = 3;            // Hamming
const int       DEFAULT_FILTER_LENGTH   = 32;
const double    DEFAULT_COMPRESSION     = 1.9937587238;
const double    DEFAULT_SCALING         = 1.0061848868;

const int       MIN_TAPS                = 2;
const int       MAX_TAPS                = 50000;


// ############################# Private Method ###############################
// sinc -- Calculates the sinc() function.
// Input:       x:              argument to sinc function
//          
// Output:                      none
// Notes:
// ############################# Private Method ###############################
double WaveletFilterCoefficients::sinc(double x)
{
  if(ABS(x) < EPS)
    return 1.;
  else
    return sin(PI*x)/PI/x;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the WaveletFilterCoefficients class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
WaveletFilterCoefficients::WaveletFilterCoefficients()
{
//
// Initialize instance variables:
//
  _lowPassCoefficients  = NULL;                 //!< The low pass filter coefficients
  _highPassCoefficients = NULL;                 //!< The high pass filter coefficients
  setFilterType(DEFAULT_FILTER_TYPE);
  setWindowType(DEFAULT_WINDOW_TYPE);
  setFilterLength(DEFAULT_FILTER_LENGTH);
  setCompressionVariable(DEFAULT_COMPRESSION);
  setScalingVariable(DEFAULT_SCALING);

  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the WaveletFilterCoefficients class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
WaveletFilterCoefficients::~WaveletFilterCoefficients()
{
  delete [] _lowPassCoefficients;
  delete [] _highPassCoefficients;
  return;
}

// ############################# Public Method ###############################
// setFilterType -- Sets a new wavelet filter type.
// Input:       type:           new wavelet filter type
//          
// Output:                      none
// Notes:
// ############################# Public Method ###############################
void WaveletFilterCoefficients::setFilterType(int type)
{
  _filterType   = type;

  return;
}

// ############################# Public Method ###############################
// setWindowType -- Sets a new wavelet window type for Sinc filter.
// Input:       type:           new window type
//          
// Output:                      none
// Notes:
// 1. See DataWindow.h for valid window types
// ############################# Public Method ###############################
void WaveletFilterCoefficients::setWindowType(int type)
{
  _windowType   = type;

  return;
}

// ############################# Public Method ###############################
// setFilterLength -- Sets a new wavelet filter coefficient size/length.
// Input:       length:         new coefficient length
//          
// Output:                      none
// Notes:
// ############################# Public Method ###############################
void WaveletFilterCoefficients::setFilterLength(int length)
{
  _filterLength = MAX(MIN_TAPS, length);
  _filterLength = MIN(_filterLength, MAX_TAPS);

  return;
}

// ############################# Public Method ###############################
// setCompressionVariable -- Sets a new sinc filter parameter.
// Input:       compression:    new compression variable
//          
// Output:                      none
// Notes:
// ############################# Public Method ###############################
void WaveletFilterCoefficients::setCompressionVariable(double compression)
{
  _compressionVariable  = compression;

  return;
}

// ############################# Public Method ###############################
// setScalingVariable -- Sets a new sinc filter parameter.
// Input:       scaling:        new scaling variable
//          
// Output:                      none
// Notes:
// ############################# Public Method ###############################
void WaveletFilterCoefficients::setScalingVariable(double scaling)
{
  _scalingVariable      = scaling;

  return;
}

// ############################# Public Method ###############################
// calculateLowPassCoefficients -- Calculates and returns the low pass filter
//                                 coefficients.
// Input:                       None
//          
// Output:                      Pointer to low pass coefficients
// Notes:
// ############################# Public Method ###############################
const float *WaveletFilterCoefficients::calculateLowPassCoefficients()
{
//
// Allocate the coefficients:
//
  delete [] _lowPassCoefficients;
  _lowPassCoefficients  = new float [_filterLength];

//
// Calculate coefficients:
//
  DataWindow *window    = new DataWindow(_windowType);
  const double *w       = window->windowFunction(_filterLength);

  int index             = 0;
  for(int n=-_filterLength/2; n<=(_filterLength-2)/2; n++)
  {
    _lowPassCoefficients[index] = sqrt(_scalingVariable/2.)*sinc((n+0.5)/_compressionVariable)*w[index];
    index++;
  }
  delete window;
  return _lowPassCoefficients;
}

// ############################# Public Method ###############################
// calculateHighPassCoefficients -- Calculates and returns the high pass filter
//                                 coefficients.
// Input:                       None
//          
// Output:                      Pointer to high pass coefficients
// Notes:
// ############################# Public Method ###############################
const float *WaveletFilterCoefficients::calculateHighPassCoefficients()
{
//
// Allocate the coefficients:
//
  delete [] _highPassCoefficients;
  _highPassCoefficients = new float [_filterLength];

//
// Calculate coefficients:
//
  DataWindow *window    = new DataWindow(_windowType);
  const double *w       = window->windowFunction(_filterLength);

  int coeff_index       = _filterLength-1;
  int window_index      = 0;
  LOGICAL even          = NO;
  for(int n=-_filterLength/2; n<=(_filterLength-2)/2; n++)
  {
    double      arg     = sqrt(_scalingVariable)/2.*sinc((n+0.5)/_compressionVariable)*w[window_index];
    if(even)
      _highPassCoefficients[coeff_index]        = arg;
    else
      _highPassCoefficients[coeff_index]        = -arg;
    even                = !even;
    window_index++;
    coeff_index--;
  }
  delete window;
  
  return _highPassCoefficients;
}
