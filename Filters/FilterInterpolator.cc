/************************************************************************
 *                                                                      *
 * This subclass of object implements a class for interpolating an      *
 * input sequence using filter coefficients (e.g, sinc filter).         *
 * See Theory of Bandlimited interpolation by J.O. Smith.               *
 * http://www-ccrma.stanford.edu.                                       *
 *                                                                      *
 * File:FilterInterpolator.h                                            *
 *                                                                      *
 * The interpolator operates as follows:                                *
 * 1. For each output point center the filter at the output point.      *
 * 2. With the filter centered at the output point, 
 *                                                                      *
 * Revision history:                                                    *
 *  1. 08/23/99  - Started.                                             *
 ************************************************************************/

#include "FilterInterpolator.h"                                 // Object prototypes
#include <math.h>
#include <stdio.h>

#define MAX(a,b)    ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)    ( ((a)>(b)) ? (b) : (a) )


// ############################# Private Method ###############################
// interpolateFilter -- Interpolates and returns filter coefficient value between
//                      samples (uses 2 point interpolation)
//
// Input:       index:          index to start interpolation
//              p:              increment from index
//          
// Output:                      interpolated value
//
//   |     |     |     |     |     |     |    |     <- filter sample times
//               ^
//               | index
//
//         ----> |  | <---- p
// ############################# Private Method ###############################
float FilterInterpolator::interpolateFilter(int index, float p)
{
  float f0, f1, f_interp;

  f0            = _filterCoefficients[index];
  index         = MIN((index+1), (_filterLength-1));
  f1            = _filterCoefficients[index];
  f_interp      = (1.-p)*f0 + p*f1;
  return f_interp;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the FilterInterpolator class.
//
// Input:       coefficients:   Set of filter coefficients
//              length:         # of coefficients
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
FilterInterpolator::FilterInterpolator(const float *coefficients, int length)
{
  if( (coefficients!=NULL) && (length>0) )
    setFilterCoefficients(coefficients, length);
  else
  {
    _filterLength       = 0;
    _filterCoefficients = NULL;
  }
  _outputSize           = 0;
  _interpolatorOutput   = NULL;

  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the FilterInterpolator class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
FilterInterpolator::~FilterInterpolator()
{
  delete [] _interpolatorOutput;
  delete [] _filterCoefficients;
  return;
}

// ############################# Public Method ###############################
// setFilterCoefficients -- Sets a new set of filter coefficients
//
// Input:       points:         new number of points
//          
// Output:                      None
// ############################# Public Method ###############################
void FilterInterpolator::setFilterCoefficients(const float *coefficients, int length)
{
  int   i;

  delete [] _filterCoefficients;
  _filterLength                 = MAX(1, length);
  _filterCoefficients           = new float[_filterLength];
  for(i=0; i<_filterLength; i++)
    _filterCoefficients[i]      = coefficients[i];
  if( (_filterLength%2) == 1)
    _filterOdd                  = YES;
  else
    _filterOdd                  = NO;
  return;
}

// ############################# Public Method ###############################
// interpolateArray -- Interpolates an input array using filter interpolation
//
// Input:       input:          array of input data
//              inputSize:      total size of array
//              startIndex:     index of input array to start interpolating
//              outputStepSize: step size of interpolated array relative to input
//                              step size == 1
//              firstStep       initialization of first step
//              outputSize:     number of points in output array
//          
// Output:                      None
//
// Notes:
// 1. Call interpolatorOutput() to get output from this interpolation.
// 2. We assume the spacing of the filter coefficients is the same as
//    the input array, i.e., delta_t=1.
// 3. The interpolator operates by placing the center of the filter at
//    each output point, and then summing the (interpolated) filter
//    coefficients times the input signal points at adjacent points.
// 3. The start index is used as follows:
//
//  |    |    |    |    |    |    |    |    |    |    |    |    |     <--- input array points
//                 ^
//              startIndex
//                    ^
//                first interpolated output point
//
//          x    x    x    x    x    x    x    <- filter coefficient locations
//
//  where we have assumed that outputStepSize ~ 0.75, and first step = 0
// ############################# Public Method ###############################
void FilterInterpolator::interpolateArray(float *input, int inputSize, int startIndex, float outputStepSize,
                                            float firstStep, int outputSize)
{
  int           i, j, input_index, half_filter;
  float         spacing, p;
  double        sum;
//
// Allocate output array:
//
  if(_outputSize < outputSize)
  {
    delete [] _interpolatorOutput;
    _interpolatorOutput         = new float [outputSize];
    _outputSize                 = outputSize;
  }
//
// Initialize the increment.
//
  p     = firstStep;
//
// Handle odd and even length filters differently:
//
  if(_filterOdd)
  {
    half_filter = (_filterLength-1)/2;
//
// Loop through number of output points:
//
    for(i=0; i<_outputSize; i++)
    {
      while(p > 1.)
      {
        startIndex++;
        p--;
      }
      if(startIndex >= inputSize)                       // check for end condition
        break;
//
// sum over the left half of filter:
//
      sum               = 0.;
      input_index       = startIndex;
      for(j=0; j<half_filter; j++)
      {
        if( YES )
          spacing       = 1.-p;
        else
          spacing       = 1.-p;
        sum             += input[input_index]*interpolateFilter(half_filter-j-1, spacing);
        input_index--;
        if(input_index < 0)
          break;
      }
      _interpolatorOutput[i]    = sum;
//
// sum over the right half of filter:
//
      sum               = 0.;
      input_index       = startIndex+1;
      for(j=0; j<half_filter; j++)
      {
        if( YES )
          spacing       = 1.-p;
        else
          spacing       = 1.-p;
        sum             += input[input_index]*interpolateFilter(half_filter+j, spacing);
        input_index++;
        if(input_index >= inputSize)
          break;
      }
      _interpolatorOutput[i]    += sum;
      p                         += outputStepSize;
    }           // end loop over i
  }
  else
  {
    half_filter = (_filterLength)/2;                    // even length filter
//
// Loop through number of output points:
//
    for(i=0; i<_outputSize; i++)
    {
      while(p > 1.)
      {
        startIndex++;
        p--;
      }
      if(startIndex >= inputSize)                       // check for end condition
        break;
//
// sum over the left half of filter:
//
      sum               = 0.;
      input_index       = startIndex;
      for(j=0; j<half_filter; j++)
      {
        if( (j%2) == 0 )
          spacing       = p;
        else
          spacing       = 1.-p;
        sum             += input[input_index]*interpolateFilter(half_filter-j, spacing);
        input_index--;
        if(input_index < 0)
          break;
      }
//
// sum over the right half of filter:
//
      input_index       = startIndex;
      for(j=0; j<half_filter; j++)
      {
        if( (j%2) == 0 )
          spacing       = p;
        else
          spacing       = 1.-p;
        sum             += input[input_index]*interpolateFilter(half_filter+j+1, spacing);
        input_index++;
        if(input_index >= outputSize)
          break;
      }
    }
    p   += outputStepSize;
  }
  return;
}


// ############################# Public Method ###############################
// filterArray -- Interpolates an input array using Lagrange interpolation
//
// Input:       input:          array of input data
//              inputSize:      total size of array
//              startIndex:     index of input array to start interpolating
//              outputStepSize: step size of interpolated array relative to input
//                              step size == 1
//              firstStep       initialization of first step
//              outputSize:     number of points in output array
//          
// Output:                      None
//
// Notes:
// 1. Call interpolatorOutput() to get output from this interpolation.
// 2. We assume the spacing of the filter coefficients is the same as
//    the input array, i.e., delta_t=1.
// 3. The interpolator operates by placing the center of the filter at
//    each output point, and then summing the (interpolated) filter
//    coefficients times the input signal points at adjacent points.
// 3. The start index is used as follows:
//
//  |    |    |    |    |    |    |    |    |    |    |    |    |     <--- input array points
//                 ^
//              startIndex
//                    ^
//                first interpolated output point
//
//          x    x    x    x    x    x    x    <- filter coefficient locations
//
//  where we have assumed that outputStepSize ~ 0.75, and first step = 0
// ############################# Public Method ###############################
void FilterInterpolator::filterArray(float *input, int inputSize, int startIndex, float outputStepSize,
                                            float firstStep, int outputSize)
{
  int           i, j, input_index, half_filter;
  float         spacing, p;
  double        sum;
//
// Allocate output array:
//
  if(_outputSize < outputSize)
  {
    delete [] _interpolatorOutput;
    _interpolatorOutput         = new float [outputSize];
    _outputSize                 = outputSize;
  }
//
// Initialize the increment.
//
  p     = firstStep;
//
// Handle odd and even length filters differently:
//
  if(_filterOdd)
  {
    half_filter = (_filterLength-1)/2;
//
// Loop through number of output points:
//
    for(i=0; i<_outputSize; i++)
    {
      _interpolatorOutput[i]    = 0.;
      startIndex++;
      if(startIndex < inputSize)                        // check for end condition
      {
//
// sum over the left half of filter:
//
        sum             = 0.;
        input_index     = startIndex;
        for(j=0; j<half_filter; j++)
        {
          sum           += input[input_index]*_filterCoefficients[half_filter-j-1];
          input_index--;
          if(input_index < 0)
            break;
        }
//
// sum over the right half of filter:
//
        input_index     = startIndex+1;
        for(j=0; j<half_filter; j++)
        {
          sum           += input[input_index]*_filterCoefficients[half_filter+j];
          input_index++;
          if(input_index >= inputSize)
            break;
        }
        _interpolatorOutput[i]  = sum;
      }
    }           // end loop over i
  }
  else
  {
    half_filter = (_filterLength)/2;                    // even length filter
//
// Loop through number of output points:
//
    for(i=0; i<_outputSize; i++)
    {
      while(p > 1.)
      {
        startIndex++;
        p--;
      }
      if(startIndex >= inputSize)                       // check for end condition
        break;
//
// sum over the left half of filter:
//
      sum               = 0.;
      input_index       = startIndex;
      for(j=0; j<half_filter; j++)
      {
        if( (j%2) == 0 )
          spacing       = p;
        else
          spacing       = 1.-p;
        sum             += input[input_index]*interpolateFilter(half_filter-j, spacing);
        input_index--;
        if(input_index < 0)
          break;
      }
//
// sum over the right half of filter:
//
      input_index       = startIndex;
      for(j=0; j<half_filter; j++)
      {
        if( (j%2) == 0 )
          spacing       = p;
        else
          spacing       = 1.-p;
        sum             += input[input_index]*interpolateFilter(half_filter+j+1, spacing);
        input_index++;
        if(input_index >= outputSize)
          break;
      }
    }
    p   += outputStepSize;
  }
  return;
}
