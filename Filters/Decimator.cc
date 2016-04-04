/************************************************************************
 *                                                                      *
 * This class implements a decimator (no filtering).                    *
 *                                                                      *
 * File:Decimator.h                                                     *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 06/08/99  - Started.                                             *
 ************************************************************************/
#include        "Decimator.h"
#include        <GNU/Complex.h>
#include        <stdio.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the Decimator class.
//
// Input:       decimation:             decimation factor
//
// Output:                              None
//
// Notes:
// ############################# Class Constructor ###############################
Decimator::Decimator(int decimation)
{
  _decimateOutput       = NULL;
  _complexOutput        = NULL;
  _oldDecimatePoints    = 0;
  _oldComplexPoints     = 0;
  setDecimateFactor(decimation);
  return;
}

// ############################# Class Destructor ###############################
// Class Constructor -- Destructor for the Decimator class.
//
// Input:                               None
//
// Output:                              None
//
// Notes:
// ############################# Class Destructor ###############################
Decimator::~Decimator()
{
  delete [] _decimateOutput;
  delete [] _complexOutput;
  return;
}

// ############################# Public Method ###############################
// decimateInput -- Decimate an input array by selecting the ith input.
//                     Return the decimated array.
// Input:       input:                  The floating point input array
//              numberPts:              Size of the array
//                      
// Output:                              The decimated array
//
// Notes:
// ############################# Public Method ###############################
float *Decimator::decimateInput(float *input, int numberPts)
{
  int   i, k, output_pts;
  
  if(_decimateFactor)
  {
    output_pts = numberPts/_decimateFactor;
    if(_oldDecimatePoints < output_pts)
    {
      delete [] _decimateOutput;
      _decimateOutput           = new float[output_pts];
      _oldDecimatePoints        = output_pts;
    }
    k = 0;
    for(i=0; i<output_pts; i++)
    {
      _decimateOutput[i]        = input[k];
      k                         += _decimateFactor;
    }
  }
  return _decimateOutput;
}


// ############################# Public Method ###############################
// decimateComplexData -- Decimate an input array by selecting the ith input.
//                     Return the decimated array.
// Input:       input:                  The complex input array
//              numberPts:              Size of the array
//                      
// Output:                              The decimated array
//
// Notes:
// ############################# Public Method ###############################
Complex *Decimator::decimateComplexData(Complex *input, int numberPts)
{
  int   i, k, output_pts;
  
  if(_decimateFactor)
  {
    output_pts = numberPts/_decimateFactor;
    if(_oldComplexPoints < output_pts)
    {
      delete [] _complexOutput;
      _complexOutput            = new Complex[output_pts];
      _oldComplexPoints         = output_pts;
    }
    k = 0;
    for(i=0; i<output_pts; i++)
    {
      _complexOutput[i] = input[k];
      k                 += _decimateFactor;
    }
  }
  return _complexOutput;
}

