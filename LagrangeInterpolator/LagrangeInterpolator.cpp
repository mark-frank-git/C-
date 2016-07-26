/************************************************************************
 *                                                                      *
 * This subclass of object implements a class for interpolating an      *
 * input sequence using Lagrange's formula.  There are also methods     *
 * for real time processing of input data samples.                      *
 * See Abramowitz & Stegun, pp. 878-879.                                *
 *                                                                      *
 * File:LagrangeInterpolator.cc                                         *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 08/23/00  - Started.                                             *
 ************************************************************************/

#include "LagrangeInterpolator.h"                                       // Object prototypes

#define MAX(a,b)    ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)    ( ((a)>(b)) ? (b) : (a) )
// ############################# Private Method ###############################
// interpolate -- Interpolates an input set of points using Lagrange interpolation
//
// Input:       f_m2:           function evaluated at -2 of index
//              f_m1:           function evaluated at -1 of index
//              f_0:            function evaluated at 0 of index
//              f_p1:           function evaluated at +1 of index
//              f_p2:           function evaluated at +2 of index
//              p:              x distance from 0 index
//          
// Output:                      The interpolated output
//
// Notes:
// ############################# Public Method ###############################
double CLagrangeInterpolator::interpolate(double f_m2, double f_m1, double f_0, double f_p1, double f_p2, double p)
{
  double output, p_squared;
//
//  Now do the interpolation based on numberInterpPoints requested:
//
  p_squared     = p*p;
  switch(m_iNumberInterpPoints)
  {
    case 2:
    default:
      output    = f_0 + p*(f_p1-f_0);
      break;
    case 3:
      output    = 0.5*p*(f_p1-f_m1) + p_squared*(0.5*f_m1 - f_0 + 0.5*f_p1) + f_0;
      break;
    case 4:
      output    = p*(1.-p)*(p-2.)*f_m1/6. + 0.5*(p_squared-1.)*(p-2.)*f_0 +
                                  0.5*p*(p+1.)*(2.-p)*f_p1 + p*(p_squared-1.)*f_p2/6.;
      break;
    case 5:
      output    = (p_squared-1.)*p*(p-2.)*f_m2/24. + (1.-p)*p*(p_squared-4.)*f_m1/6. +
                                  0.25*(p_squared-1.)*(p_squared-4.)*f_0 + (p+1.)*p*(4-p_squared)*f_p1/6. +
                                  (p_squared-1.)*p*(p+2.)*f_p2/24.;
      break;
  }
  return output;
}

// ############################# Private Method ###############################
// interpolateDouble -- Interpolates an input set of points using Lagrange interpolation
//
// Input:       f[]:            function evaluated at -/+ input points
//              p:              x distance from 0 index
//          
// Output:                      The interpolated output
//
// Notes:
// ############################# Public Method ###############################
double CLagrangeInterpolator::interpolateDouble(double *f, double p)
{
  double output, p_squared;
  double        f_m2, f_m1, f_0, f_p1, f_p2, f_p3;
//
//  Now do the interpolation based on numberInterpPoints requested:
//
  p_squared     = p*p;
  switch(m_iNumberInterpPoints)
  {
    case 2:
    default:
      f_0       = f[0];
      f_p1      = f[1];
      output    = f_0 + p*(f_p1-f_0);
      break;
    case 3:
      f_m1      = f[0];
      f_0       = f[1];
      f_p1      = f[2];
      output    = 0.5*p*(f_p1-f_m1) + p_squared*(0.5*f_m1 - f_0 + 0.5*f_p1) + f_0;
      break;
    case 4:
      f_m1      = f[0];
      f_0       = f[1];
      f_p1      = f[2];
      f_p2      = f[3];
      output    = p*(1.-p)*(p-2.)*f_m1/6. + 0.5*(p_squared-1.)*(p-2.)*f_0 +
                                  0.5*p*(p+1.)*(2.-p)*f_p1 + p*(p_squared-1.)*f_p2/6.;
      break;
    case 5:
      f_m2      = f[0];
      f_m1      = f[1];
      f_0       = f[2];
      f_p1      = f[3];
      f_p2      = f[4];
      output    = (p_squared-1.)*p*(p-2.)*f_m2/24. + (1.-p)*p*(p_squared-4.)*f_m1/6. +
                                  0.25*(p_squared-1.)*(p_squared-4.)*f_0 + (p+1.)*p*(4-p_squared)*f_p1/6. +
                                  (p_squared-1.)*p*(p+2.)*f_p2/24.;
      break;
    case 6:
      f_m2      = f[0];
      f_m1      = f[1];
      f_0       = f[2];
      f_p1      = f[3];
      f_p2      = f[4];
      f_p3      = f[5];
      output    = p*(1.-p_squared)*(p-2.)*(p-3.)*f_m2/120. + p*(p-1.)*(p_squared-4.)*(p-3.)*f_m1/24. +
                                  (3.-p)*(p_squared-1.)*(p_squared-4.)*f_0/12. + p*(p+1.)*(p_squared-4.)*(p-3.)*f_p1/12. +
                                  p*(p_squared-1.)*(p+2.)*(3.-p)*f_p2/24.+p*(p_squared-1.)*(p_squared-4.)*f_p3/120.;
      break;
  }
  return output;
}



// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the CLagrangeInterpolator class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
CLagrangeInterpolator::CLagrangeInterpolator(int numberPoints)
{
  setNumberInterpPoints(numberPoints);

  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the CLagrangeInterpolator class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
CLagrangeInterpolator::~CLagrangeInterpolator()
{
   return;
}

// ############################# Public Method ###############################
// setNumberInterpPoints -- Sets a new number points for interpolation
//
// Input:       points:         new number of points
//          
// Output:                      None
// ############################# Public Method ###############################
void CLagrangeInterpolator::setNumberInterpPoints(int points)
{
  m_iNumberInterpPoints   = MAX(MIN_INTERP_POINTS, points);
  m_iNumberInterpPoints   = MIN(m_iNumberInterpPoints, MAX_INTERP_POINTS);
  return;
}

// ############################# Public Method ###############################
// interpolateArray -- Interpolates an input array using Lagrange interpolation
//
// Input:       input:          array of input data
//              inputSize:      total size of array
//              startIndex:     index of input array to start interpolating
//              outputStepSize: step size of interpolated array relative to input
//                              step size == 1
//              firstStep       initialization of first step output step (set to 0 or outputStepSize)
//              outputSize:     number of points in output array
//          
// Output:      output:         interpolated array
//
// Notes:
// 1. Call interpolatorOutput() to get output from this interpolation.
// 2. The start index is used as follows:
//
//  |    |    |    |    |    |    |    |    |    |    |    |    |     <--- input array points
//                 ^
//              startIndex
//                    ^
//                first interpolated output point
//  where we have assumed that outputStepSize ~ 0.75
// ############################# Public Method ###############################
void CLagrangeInterpolator::interpolateArray(const float *input, int inputSize, int startIndex, double outputStepSize,
                                            double firstStep, double *output, int outputSize)
{
  int   i, index;
  double f_m2, f_m1, f_0, f_p1, f_p2;
  double p;
//
// Initialize the increment.
//
  p     = firstStep;
//
// Loop through number of output points:
//
  for(i=0; i<outputSize; i++)
  {
    while(p > 1.)
    {
      startIndex++;
      p--;
    }
    if(startIndex >= inputSize)                         // check for end condition
      break;
//
// Set the history and future function values, start index = x0.
//
    index       = MAX((startIndex - 2), 0);
    f_m2        = input[index];
    index       = MAX((startIndex - 1), 0);
    f_m1        = input[index];
    f_0         = input[startIndex];
    index       = MIN((startIndex+1), (inputSize-1));
    f_p1        = input[index];
    index       = MIN((startIndex+2), (inputSize-1));
    f_p2        = input[index];
//
//  Now do the interpolation based on numberInterpPoints requested:
//
    output[i]      = interpolate(f_m2, f_m1, f_0, f_p1, f_p2, p);
    p   += outputStepSize;
  }
  return;
}

// ############################# Public Method ###############################
// interpolateArray -- Interpolates an unevenly spaced input array using Lagrange interpolation
//
// Input:       inputYData:     array of input data
//              inputXData:     array of input sample points (possibly uneven spacing)
//              outputXData:    array of output samples points to perform interpolation
//              inputSize:      total size of array
//              outputSize:     number of points in output array
//          
// Output:                      None
//
// Notes:
// 1. Call interpolatorOutput() to get output from this interpolation.
// ############################# Public Method ###############################
void CLagrangeInterpolator::interpolateUnevenArray(const float *inputYData, const float *inputXData, const float *outputXData, 
                                                  int inputSize, double *output, int outputSize)
{
  int   i, j, last_j, index;
  float f_m2, f_m1, f_0, f_p1, f_p2;
  float p, delta;
//
// Allocate output array:
//
//
// Loop through number of output points:
//
  last_j        = 0;
  for(i=0; i<outputSize; i++)
  {
//
// Search for first input x data point less than output x data point
    for(j=last_j; j<inputSize; j++)
    {
      if(inputXData[j] > outputXData[i])
        break;
    }
    j           = MAX((j-1), 0);
    j           = MIN(j, (inputSize-1));
    last_j      = j;
    p           = outputXData[i] - inputXData[j];
    if(j<inputSize-1)
    {
      delta     = inputXData[j+1] - inputXData[j];
      if(delta > 0.)
        p       /= delta;
    }
//
// Set the history and future function values, start index = x0.
//
    index       = MAX((j - 2), 0);
    f_m2        = inputYData[index];
    index       = MAX((j - 1), 0);
    f_m1        = inputYData[index];
    f_0         = inputYData[j];
    index       = MIN((j+1), (inputSize-1));
    f_p1        = inputYData[index];
    index       = MIN((j+2), (inputSize-1));
    f_p2        = inputYData[index];
//
//  Now do the interpolation based on numberInterpPoints requested:
//
    output[i]      = interpolate(f_m2, f_m1, f_0, f_p1, f_p2, p);
  }
  return;
}
