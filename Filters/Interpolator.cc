/************************************************************************
 *                                                                      *
 * This class implements an interpolator (no filtering).                *
 *                                                                      *
 * File:Interpolator.cc                                                 *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 06/08/02  - Started.                                             *
 ************************************************************************/
#include        "Interpolator.h"
#include        <stdio.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the Interpolator class.
//
// Input:       interpolation:          interpolation factor
//
// Output:                              None
//
// Notes:
// ############################# Class Constructor ###############################
Interpolator::Interpolator(int interpolation)
{
  _interpolateOutput    = NULL;
  setInterpolateFactor(interpolation);
  return;
}

// ############################# Class Destructor ###############################
// Class Constructor -- Destructor for the Interpolator class.
//
// Input:                               None
//
// Output:                              None
//
// Notes:
// ############################# Class Destructor ###############################
Interpolator::~Interpolator()
{
  delete [] _interpolateOutput;
  return;
}

// ############################# Public Method ###############################
// interpolateInput -- Interpolate an input array by selecting the ith input.
//                     Return the interpolated array.
// Input:       input:                  The floating point input array
//                      numberPts:              Size of the array
//                      
// Output:                                      The interpolated array
//
// Notes:
// ############################# Public Method ###############################
float *Interpolator::interpolateInput(float *input, int numberPts)
{
  int   i, j, k, output_pts;
  float  temp;
  
  output_pts = numberPts*_interpolateFactor;
  if(_oldInterpolatePoints < output_pts)
  {
    if(_interpolateOutput != NULL)
      delete [] _interpolateOutput;
    _interpolateOutput  = new float[output_pts];
    _oldInterpolatePoints       = output_pts;
  }

  k             = 0;
  temp  = (float)_interpolateFactor;
  for(i=0; i<numberPts; i++)
  {
    _interpolateOutput[k++] = temp*input[i];
    for(j=1; j<_interpolateFactor; j++)
      _interpolateOutput[k++] = 0.;
  }
  return _interpolateOutput;
}
