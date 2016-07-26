/************************************************************************
 *                                                                      *
 * This is used for finding the peak position of a set of FFT bins.     *
 * The methods used here are from DSP GURU.                             *
 *                                                                      *
 * File: /User/frank/C++/Measurements/PeakDetect.cc                     *
 *                                                                      *
 ************************************************************************/

#include "PeakDetect.h"
#include <stdio.h>
#include <math.h>

#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define MAGNITUDE_SQ(a, b)      ( (a)*(a) + (b)*(b) )
#ifndef YES
#define YES 1
#define NO  0
#endif


// ############################ Private Function ###################################
// apFN - Local function for evaluating Quinn's methods.
//
// Input:       realData1:      real data point
//              reaData0:       real data point at max_index
//              imagData1:      imag data point
//              imagData0:      imag data point at max_index
//
// Output:                      function output
//
// Notes:
// ############################ Private Function ###################################
double PeakDetect::apFN(float realData1, float realData0, float imagData1, float imagData0)
{
  double ap_fn, den;
//
  ap_fn         = realData1*realData0 + imagData1*imagData0;
  den           = realData0*realData0 + imagData0*imagData0;
  if(den > 0.)
    ap_fn       /= den;
  return ap_fn;
}

// ############################ Private Function ###################################
// tauFN - Local function for evaluating Quinn's methods.
//
// Input:       x:              real data point
//
// Output:                      function output
//
// Notes:
// ############################ Private Function ###################################
double PeakDetect::tauFN(float x)
{
  double tau_fn, second_term, arg, den;
//
/***************************** old mistake (works better?) ****************
  tau_fn        = second_term   = 0.;
  arg           = 3*x*x + 6*x + 1;
  if(arg > 0.)
    tau_fn      = log(arg)*0.25;
  arg           = x + 1. - sqrt(2./3.);
  if(arg > 0.)
    second_term = log(arg);
  second_term   *= sqrt(6.)/24.;
  arg           = x + 1. + sqrt(2./3.);
  if(arg > 0.)
    second_term /= arg;
  tau_fn        -= second_term;
*****************************************************************************/

  tau_fn        = second_term   = 0.;
  arg           = 3*x*x + 6*x + 1;
  if(arg > 0.)
    tau_fn      = log(arg)*0.25;
  arg           = x + 1. - sqrt(2./3.);
  den           = x + 1. + sqrt(2./3.);
  if(den != 0.)
    second_term = log(arg/den);
  second_term   *= sqrt(6.)/24.;
  tau_fn        -= second_term;

  return        tau_fn;
}

// ############################# Class Constructor #################################
// PeakDetect -- Constructor for the PeakDetect class
//
// Input:       fs:     Sampling frequency
//
// Output:              None
// ############################# Class Constructor #################################
PeakDetect::PeakDetect()
{
// 
// Initialize instance variables:
//
  setPeakDetectType(NO_INTERPOLATE);
  return;
}


// ############################# Class Destructor ###############################
// PeakDetect -- Destructor for the PeakDetect class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
PeakDetect::~PeakDetect()
{
  return;
}

// ############################ Public Function ###################################
// peakPosition - Find and return the peak position from a set of FFT data.
//
// Input:       realFFTData:    Array of real parts of FFT data
//              imagFFTData:    Array of imag parts of FFT data
//              lowIndex:       Low index of array to start search
//              highIndex:      High index of array to end search
//
// Output:                      Interpolated index of peak position
//
// Notes:
// ############################ Public Function ###################################
float PeakDetect::peakPosition(float *realFFTData, float *imagFFTData, int lowIndex, int highIndex)
{
  int   i, max_index;
  int   k_plus_1, k_minus_1;
  float ap, dp, am, dm, d;
  float mag, peak_posn, max_mag, den;
  float y1, y2, y3, a;
//
// Error check indices:
//
  if(highIndex <= lowIndex)
    return highIndex;
//
//  Now, find maximum bin:
//
  max_index     = lowIndex;
  max_mag       = MAGNITUDE_SQ(realFFTData[max_index], imagFFTData[max_index]);
  for(i=lowIndex+1; i<=highIndex; i++)
  {
    mag         = MAGNITUDE_SQ(realFFTData[i],imagFFTData[i]);
    if(mag > max_mag)
    {
      max_mag   = mag;
      max_index = i;
    }
  }
//
// Now, interpolate using one of the methods from DSP GURU:
//
  k_minus_1     = MAX(lowIndex, (max_index-1));
  k_plus_1      = MIN(highIndex, (max_index+1));
  d             = 0.;
  switch(m_detectType)
  {
    case NO_INTERPOLATE:
    default:
      break;
    case QUADRATIC_METHOD:
      y1                = MAGNITUDE_SQ(realFFTData[k_minus_1], imagFFTData[k_minus_1]);
      y1                = sqrt(y1);
      y2                = sqrt(max_mag);
      y3                = MAGNITUDE_SQ(realFFTData[k_plus_1], imagFFTData[k_plus_1]);
      y3                = sqrt(y3);
      den               = 2.*(y2+y2-y1-y3);
      if(den != 0.)
        d               = (y3-y1)/den;
      break;
    case BARYCENTRIC_METHOD:
      y1                = MAGNITUDE_SQ(realFFTData[k_minus_1], imagFFTData[k_minus_1]);
      y1                = sqrt(y1);
      y2                = sqrt(max_mag);
      y3                = MAGNITUDE_SQ(realFFTData[k_plus_1], imagFFTData[k_plus_1]);
      y3                = sqrt(y3);
      den               = (y1+y2+y3);
      if(den != 0.)
        d               = (y3-y1)/den;
      break;
    case QUINNS_FIRST:
      dp                = dm    = 0.;
      ap                = apFN(realFFTData[k_plus_1], realFFTData[max_index], imagFFTData[k_plus_1],
                                 imagFFTData[max_index]);
      den               = 1. - ap;
      if(den != 0.)
        dp              = -ap/den;
      am                = apFN(realFFTData[k_minus_1], realFFTData[max_index], imagFFTData[k_minus_1],
                                 imagFFTData[max_index]);
      den               = 1. - am;
      if(den != 0.)
        dm              = am/den;
      if( (dp>0.) && (dm>0.) )
        d               = dp;
      else
        d               = dm;
      break;
    case QUINNS_SECOND:
      dp                = dm    = 0.;
      ap                = apFN(realFFTData[k_plus_1], realFFTData[max_index], imagFFTData[k_plus_1],
                                 imagFFTData[max_index]);
      den               = 1. - ap;
      if(den != 0.)
        dp              = -ap/den;
      am                = apFN(realFFTData[k_minus_1], realFFTData[max_index], imagFFTData[k_minus_1],
                                 imagFFTData[max_index]);
      den               = 1. - am;
      if(den != 0.)
        dm              = am/den;
      d                 = (dp+dm)/2. + tauFN(dp*dp) - tauFN(dm*dm);
      break;
    case JAINS_METHOD:
      y1                = MAGNITUDE_SQ(realFFTData[k_minus_1], imagFFTData[k_minus_1]);
      y1                = sqrt(y1);
      y2                = sqrt(max_mag);
      y3                = MAGNITUDE_SQ(realFFTData[k_plus_1], imagFFTData[k_plus_1]);
      y3                = sqrt(y3);
      if( (y1>y3) && (y1>0.) )
      {
        a               = y2/y1;
        den             = a+1.;
        if(den > 0.)
          d             = a/den;
        d               -= 1.;
      }
      else if( y2>0. )
      {
        a               = y3/y2;
        den             = a+1.;
        if(den > 0.)
          d             = a/den;
      }
      break;
  }
  peak_posn             = max_index + d;
  return peak_posn;
}