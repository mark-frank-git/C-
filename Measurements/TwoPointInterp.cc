/************************************************************************
 *                                                                      *
 * This is used for finding the peak location from two points.  The     *
 * algorithm assumes the underlying function is triangle shaped, like   *
 * a PN autocorrelation.                                                *
 *                                                                      *
 * File: /User/frank/C++/Measurements/TwoPointInterp.h                  *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 07/17/99 - Started.                                              *
 ************************************************************************/

#include "TwoPointInterp.h"
#include <stdio.h>
#include <math.h>

#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )

#ifndef YES
#define YES 1
#define NO  0
#endif


// ############################# Class Constructor #################################
// TwoPointInterp -- Constructor for the TwoPointInterp class
//
// Input:       fs:     Sampling frequency
//
// Output:              None
// ############################# Class Constructor #################################
TwoPointInterp::TwoPointInterp()
{
// 
// Initialize instance variables:
//
  setTriangleWidth(1.5);
  return;
}


// ############################# Class Destructor ###############################
// TwoPointInterp -- Destructor for the TwoPointInterp class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
TwoPointInterp::~TwoPointInterp()
{
  return;
}

// ############################ Public Function ###################################
// setTriangleWidth - Sets the width of the triangle base.
//
// Input:       width:          New width of triangle base
//
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void TwoPointInterp::setTriangleWidth(float width)
{
  _triangleWidth        = MAX(MIN_TRIANGLE_WIDTH, width);
  return;
}

// ############################ Public Function ###################################
// peakPosition - Find and return the peak position from a set of autocorrelation
//                data.
//
// Input:       autoData:       Array of real parts of FFT data
//              lowIndex:       Low index of array to start search
//              highIndex:      High index of array to end search
//
// Output:                      Interpolated index of peak position
//
// Notes:
// 1. This function uses the algorithm described in the Magna Simulation document.
//    It assumes a triangular shaped autocorrelation function with height offset
//    = noise level.
// ############################ Public Function ###################################
float TwoPointInterp::peakPosition(float *autoData, int lowIndex, int highIndex)
{
  int   i, max_index, noise_count;
  int   int_width;
  int   x1, x2;
  float noise_level, slope;
  float b1, b2, y1, y2;
  float max_mag, mag, denominator, x_peak;
//
// Error check indices:
//
  if(highIndex <= lowIndex)
    return highIndex;
//
//  Now, find maximum bin:
//
  max_index     = lowIndex;
  max_mag       = autoData[max_index];
  for(i=lowIndex+1; i<=highIndex; i++)
  {
    mag         = autoData[i];
    if(mag > max_mag)
    {
      max_mag   = mag;
      max_index = i;
    }
  }
//
// Now, define two high points for interpolation:
//
  if(max_index == lowIndex)
  {
    x1  = lowIndex;
    x2  = lowIndex+1;
  }
  else if (max_index == highIndex)
  {
    x1  = highIndex - 1;
    x2  = highIndex;
  }
  else
  {
    if(autoData[max_index+1] > autoData[max_index-1])
    {
      x1        = max_index;
      x2        = max_index+1;
    }
    else
    {
      x1        = max_index-1;
      x2        = max_index;
    }
  }
  y1    = autoData[x1];
  y2    = autoData[x2];
//
// Now find noise level from points outside of signal autocorrelation:
//
  
  int_width     = (int) (_triangleWidth/2.+1.);
  noise_level   = 0.;
  noise_count   = 0;
  for(i=lowIndex; i<max_index-int_width; i++)
  {
    noise_level += autoData[i];
    noise_count++;
  }
  for(i=max_index+int_width; i<highIndex; i++)
  {
    noise_level += autoData[i];
    noise_count++;
  }
  if(noise_count > 0)
    noise_level /= noise_count;
//
// Now, find the slope of the lines through the two points:
//
  denominator   = x2 - x1 - _triangleWidth;
  slope         = 100.;
  if(denominator != 0.)
    slope       = (2.*noise_level - y1 - y2)/denominator;
//
// Next, find the y intercepts of the two lines:
//
  b1            = y1 - slope*x1;
  b2            = y2 + slope*x2;
//
// Finally, find the peak location:
//
  x_peak        = max_index;
  if(slope != 0.)
    x_peak      = (b2 - b1)/2./slope;
  return x_peak;
}