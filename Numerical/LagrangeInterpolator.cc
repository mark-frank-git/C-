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
 ************************************************************************/

#include "LagrangeInterpolator.h"                                       // Object prototypes
#include <math.h>
#include <stdio.h>
#include <Buffers/DoubleBuffer.h>

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
double LagrangeInterpolator::interpolate(float f_m2, float f_m1, float f_0, float f_p1, float f_p2, float p)
{
  double output, p_squared;
//
//  Now do the interpolation based on numberInterpPoints requested:
//
  p_squared     = p*p;
  switch(_numberInterpPoints)
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
double LagrangeInterpolator::interpolateDouble(double *f, double p)
{
  double output, p_squared;
  double        f_m2, f_m1, f_0, f_p1, f_p2, f_p3;
//
//  Now do the interpolation based on numberInterpPoints requested:
//
  p_squared     = p*p;
  switch(_numberInterpPoints)
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
// Class Constructor -- Constructor for the LagrangeInterpolator class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
LagrangeInterpolator::LagrangeInterpolator(int numberPoints)
{
  _outputSize           = 0;
  _interpolatorOutput   = NULL;
  _inputBuffer          = NULL;
  setNumberInterpPoints(numberPoints);
  setInputTimeStep(DEFAULT_TIME_STEP);
  setOutputTimeStep(DEFAULT_TIME_STEP);
  reset();

  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the LagrangeInterpolator class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
LagrangeInterpolator::~LagrangeInterpolator()
{
  delete [] _interpolatorOutput;
  delete _inputBuffer;
  return;
}

// ############################# Public Method ###############################
// reset -- Resets the interpolator for real time processing
//
// Input:                       None
//          
// Output:                      None
// ############################# Public Method ###############################
void LagrangeInterpolator::reset()
{
  _inputSample  = _outputSample = 0;
  _inputBuffer->resetBuffer();
  return;
}

// ############################# Public Method ###############################
// setNumberInterpPoints -- Sets a new number points for interpolation
//
// Input:       points:         new number of points
//          
// Output:                      None
// ############################# Public Method ###############################
void LagrangeInterpolator::setNumberInterpPoints(int points)
{
  _numberInterpPoints   = MAX(MIN_INTERP_POINTS, points);
  _numberInterpPoints   = MIN(_numberInterpPoints, MAX_INTERP_POINTS);
  delete _inputBuffer;
  _inputBuffer  = new DoubleBuffer(_numberInterpPoints);
  return;
}

// ############################# Public Method ###############################
// setInputSampleTime -- Sets a new input sampling time for real time processing
//
// Input:       time:           New input sample time in seconds
//          
// Output:                      None
// ############################# Public Method ###############################
void LagrangeInterpolator::setInputTimeStep(float time)
{
  _inputTimeStep        = time;
  return;
}

// ############################# Public Method ###############################
// setOutputSampleTime -- Sets a new output sampling time for real time processing
//
// Input:       time:           New output sample time in seconds
//          
// Output:                      None
// ############################# Public Method ###############################
void LagrangeInterpolator::setOutputTimeStep(float time)
{
  _outputTimeStep       = time;
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
//              firstStep       initialization of first step
//              outputSize:     number of points in output array
//          
// Output:                      None
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
void LagrangeInterpolator::interpolateArray(const float *input, int inputSize, int startIndex, float outputStepSize,
                                            float firstStep, int outputSize)
{
  int   i, index;
  float f_m2, f_m1, f_0, f_p1, f_p2;
  float p;
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
// Loop through number of output points:
//
  for(i=0; i<_outputSize; i++)
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
    _interpolatorOutput[i]      = interpolate(f_m2, f_m1, f_0, f_p1, f_p2, p);
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
void LagrangeInterpolator::interpolateArray(const float *inputYData, const float *inputXData, const float *outputXData, int inputSize,
                                            int outputSize)
{
  int   i, j, last_j, index;
  float f_m2, f_m1, f_0, f_p1, f_p2;
  float p, delta;
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
// Loop through number of output points:
//
  last_j        = 0;
  for(i=0; i<_outputSize; i++)
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
    _interpolatorOutput[i]      = interpolate(f_m2, f_m1, f_0, f_p1, f_p2, p);
  }
  return;
}

// ############################# Public Method ###############################
// numberOfInputPointsNeeded -- Returns the number of input points needed to
//                              produce next output sample for real time proc.
//
// Input:                       None
//          
// Output:                      # of samples needed
//
// Notes:
// 1. After this routine is called, supply the required number of samples using
//    bufferNextInput().
// 2. The variable, _inputSample tracks the number of input samples received.
// ############################# Public Method ###############################
int LagrangeInterpolator::numberOfInputPointsNeeded()
{
  int           input_samples;
  double        output_time;
//
// Check for boundary condition:
//
  if(_outputSample == 0)
    return _numberInterpPoints/2 + 1;
//
// The number of input samples depends on the output time, and the number of
// interpolation points:
//
  output_time   = _outputSample * _outputTimeStep;
//
// Calculate total number of input samples needed
//
  input_samples = (int)(output_time/_inputTimeStep);
  input_samples += 1 + _numberInterpPoints/2;
//
// Number of samples needed = total number - number received
//
  return        MAX(0, (input_samples - _inputSample));
}

// ############################# Public Method ###############################
// bufferNextInput -- This function receives a new input sample, and puts it
//                      into the buffer
//
// Input:       inputData:      Next input sample
//          
// Output:                      None
//
// Notes:
// ############################# Public Method ###############################
void    LagrangeInterpolator::bufferNextInput(double inputData)
{
  int   i;
  _inputSample++;
  _inputBuffer->writeToBuffer(&inputData);
  if(_inputSample==1)                           // Fill up the buffer if first sample
    for(i=0; i<_numberInterpPoints-1; i++)
      _inputBuffer->writeToBuffer(&inputData);
  return;
}

// ############################# Public Method ###############################
// getNextOutput -- This function outputs the next interpolated output
//
// Input:                       None
//          
// Output:                      Next interpolated data point
//
// Notes:
// 1. The functions, numberOfInputPointsNeeded() and bufferNextInput() need to
//    be called prior to calling this function.
// ############################# Public Method ###############################
double LagrangeInterpolator::getNextOutput()
{
  int           i;
  const         double  *buffer_output;
  double        f[MAX_INTERP_POINTS];
  double        p;
  double        input_time, output_time;
//
// We first need to find, p, the percentage step of the ouput point between
// input samples:
//
  output_time   = _outputSample * _outputTimeStep;
  input_time    = (_inputSample - _numberInterpPoints/2)*_inputTimeStep;
  p             = (output_time - input_time)/_inputTimeStep;
//
// Update output counter:
//
  _outputSample++;
//
// Get the output from the input buffer.
//
  _inputBuffer->setReadPointerToWritePointer();
  for(i=0; i<_numberInterpPoints; i++)
  {
    buffer_output       = _inputBuffer->readFromBuffer();
    f[i]                = buffer_output[0];
  }
//
// Return the interpolated value:
//
  return interpolateDouble(f, p);
}
