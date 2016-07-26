/************************************************************************
 *                                                                      *
 * This class implements a digital FIR filter having real coefficients  *
 *                                                                      *
 * File:RealFIR.cc                                                      *
 *                                                                      *
 * The filter is stored in the forms:                                   *
 *                                                                      *
 *    H(z) =  b[0] + b[1]z^(-1) + ... + b[n]z^-(n)                      *
 *                                                                      *
 * Where the b coefficients are calculated outside of this class.       *
 *                                                                      *
 *                                                                      *
 ************************************************************************/
#include <stdio.h>
#include "RealFIR.h"                                    // Object prototypes
#if defined(WIN32)
#include <GNU/Complex.h>
#else
#include "Complex.h"
#endif

#include <math.h>
#include <stdio.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the RealFIR class.
//
// Input:       coeffs:         Filter coefficients
//              length:         Number of coefficients
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
RealFIR::RealFIR(const float *coeffs, int length)
{
//
// Initialize instance variables:
//
  _numberTaps           = 0;
  _shiftPtr             = 0;
  _numberRealOutput     = _oldRealOutput        = 0;
  _numberComplexOutput  = _oldComplexOutput     = 0;

  _filterCoefficients   = NULL;
  _realShiftBuffer      = NULL;
  _complexShiftBuffer   = NULL;
  _realOutputData       = NULL;
  _complexOutputData    = NULL;

//
// Load coefficients, if available:
//
  if(coeffs != NULL)
    setFilterCoefficients(coeffs, length);

  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the RealFIR class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
RealFIR::~RealFIR()
{
  delete [] _filterCoefficients;
  delete [] _realShiftBuffer;
  delete [] _complexShiftBuffer;
  delete [] _realOutputData;
  delete [] _complexOutputData;
  
  return;
}

// ############################# Public Method ###############################
// setFilterCoefficients -- Sets a new set of coefficients for the filter.
// Input:       coeffs:         filter coefficients
//              length:         number of taps
//          
// Output:                      none
// Notes:
// ############################# Public Method ###############################
void RealFIR::setFilterCoefficients(const float *coeffs, int length)
{
  _numberTaps           = MAX(MIN_TAPS, length);
  _numberTaps           = MIN(MAX_TAPS, _numberTaps);
  _shiftSize            = _numberTaps - 1;
//
// Allocate arrays based on new length
//
  delete []     _filterCoefficients;
  delete []     _realShiftBuffer;
  delete []     _complexShiftBuffer;

  _filterCoefficients   = new float[_numberTaps];
  _realShiftBuffer      = new float[_numberTaps-1];
  _complexShiftBuffer   = new Complex[_numberTaps-1];
//
// load the coefficients
//
  for(int i=0; i<_numberTaps; i++)
  {
    _filterCoefficients[i]      = coeffs[i];
  }

  return;
}

// ############################# Public Method ###############################
// reset -- Resets the filter's shift register.
// Input:                       None
//          
// Output:                      None
// ############################# Public Method ###############################
void RealFIR::reset()
{
//
// Perform super's function
//
  _shiftPtr     = MAX((_shiftSize - 1), 0);     // The pointer decrements, so set to end
  for(int i=0; i<_shiftSize; i++)
  {
    _realShiftBuffer[i]         = 0.;
    _complexShiftBuffer[i]      = Complex(0., 0.);
  }
  
  return;
}

// ############################# Public Method ###############################
// filterFloatArray -- Filter the input array and save it in _realOutputData.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   x:              An array of floating point data.
//          numberPts:      Number of points in the input array.
//          
// Output:                  YES if went OK, No if _bPolyObject doesn't exist.
// Notes:
// 1. Call realOutputData() to get results
// ############################# Public Method ###############################
void RealFIR::filterFloatArray(const float * input, int numberPts)
{
  int       shift_ptr;
  double    yout;
  const     float *b;

//
//  allocate output array
//
  _numberRealOutput     = numberPts;
  if(_oldRealOutput < _numberRealOutput)
  {
    _oldRealOutput      = _numberRealOutput;
    delete [] _realOutputData;
    _realOutputData     = new float[_numberRealOutput];
  }
//
// Filter the data:
//
  for(int j=0; j<numberPts; j++)
  {
    b           = _filterCoefficients;
    yout        = (*b) * input[j];
    b++;
    shift_ptr   = _shiftPtr;
    for(int i=0; i<_shiftSize; i++)
    {
      yout      += (*b) * _realShiftBuffer[shift_ptr];
      b++;
      shift_ptr++;
      if(shift_ptr == _shiftSize)
        shift_ptr       = 0;
    }
//
// Add the new data to the shift register
//
    _shiftPtr--;
    if(_shiftPtr < 0)
      _shiftPtr                 = MAX((_shiftSize - 1), 0);
    _realShiftBuffer[_shiftPtr] = input[j];
//
// Save the output:
//
    _realOutputData[j]          = (float)yout;
  }

  return;
}

// ############################# Public Method ###############################
// filterFloatData -- Filter an input data sample, and return the sampled value.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   input:          An input data sample
//          
// Output:                  The filtered data sample
// ############################# Public Method ###############################
float RealFIR::filterFloatData(float input)
{
  int           shift_ptr;
  double        yout;
  const         float *b;

//
// Now, get the coefficients
//
  b         = _filterCoefficients;
//
// Filter the data:
//
  yout          = (*b) * input;
  b++;
  shift_ptr     = _shiftPtr;
  for(int i=0; i<_shiftSize; i++)
  {
    yout        += (*b) * _realShiftBuffer[shift_ptr];
    b++;
    shift_ptr++;
    if(shift_ptr == _shiftSize)
      shift_ptr = 0;
  }
//
// Add the new data to the shift register
//
  _shiftPtr--;
  if(_shiftPtr < 0)
    _shiftPtr                   = MAX((_shiftSize - 1), 0);
  _realShiftBuffer[_shiftPtr]   = input;

  return (float)yout;
}

// ############################# Public Method ###############################
// filterComplexArray -- Filter the input array and save in _complexOutputData
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   x:              An array of Complex data
//          numberPts:      Number of points in the input array.
//          
// Output:                  None.
// ############################# Public Method ###############################
void RealFIR::filterComplexArray(const Complex *input, int numberPts)
{
  int       shift_ptr;
  Complex    yout;
  const     float *b;

//
//  allocate output array
//
  _numberComplexOutput  = numberPts;
  if(_oldComplexOutput < _numberComplexOutput)
  {
    _oldComplexOutput   = _numberComplexOutput;
    delete [] _complexOutputData;
    _complexOutputData  = new Complex[_numberComplexOutput];
  }
//
// Filter the data:
//
  for(int j=0; j<numberPts; j++)
  {
    b           = _filterCoefficients;
    yout        = (*b) * input[j];
    b++;
    shift_ptr   = _shiftPtr;
    for(int i=0; i<_shiftSize; i++)
    {
      yout      += (*b) * _complexShiftBuffer[shift_ptr];
      b++;
      shift_ptr++;
      if(shift_ptr == _shiftSize)
        shift_ptr       = 0;
    }
//
// Add the new data to the shift register
//
    _shiftPtr--;
    if(_shiftPtr < 0)
      _shiftPtr                         = MAX((_shiftSize - 1), 0);
    _complexShiftBuffer[_shiftPtr]      = input[j];
//
// Save the output:
//
    _complexOutputData[j]               = yout;
  }

  return;
}

// ############################# Public Method ###############################
// filterComplexData -- Filter a complex data sample, and return the sampled value.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   xin:            An input data sample
//          
// Output:                  The filtered data sample
// ############################# Public Method ###############################
Complex  RealFIR::filterComplexData(const Complex &input)
{
  int           shift_ptr;
  Complex       yout;
  const         float *b;

//
// Now, get the coefficients
//
  b         = _filterCoefficients;
//
// Filter the data:
//
  yout          = (*b) * input;
  b++;
  shift_ptr     = _shiftPtr;
  for(int i=0; i<_shiftSize; i++)
  {
    yout        += (*b) * _complexShiftBuffer[shift_ptr];
    b++;
    shift_ptr++;
    if(shift_ptr == _shiftSize)
      shift_ptr = 0;
  }
//
// Add the new data to the shift register
//
  _shiftPtr--;
  if(_shiftPtr < 0)
    _shiftPtr                   = MAX((_shiftSize - 1), 0);
  _complexShiftBuffer[_shiftPtr]        = input;

  return yout;
}

