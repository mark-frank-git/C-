/************************************************************************
 *                                                                      *
 * This subclass of object implements a class for interpolating an      *
 * input sequence using a windowed sinc pulse.                          *
 *                                                                      *
 * File:SincInterpolator.cc                                             *
 *                                                                      *
 *                                                                      *
 ************************************************************************/

#include "SincInterpolator.h"                                   // Object prototypes
#include <math.h>
#include <stdio.h>
#include <Buffers/DoubleBuffer.h>
#include <Filters/FIRFilter.h>

#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define ROUND(a)        ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )

// ############################# Private Function ###############################
// calculateSincFunction -- Calculates a new sinc function for the class parameters.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// 1. The sinc function is based on the input sampling rate times the oversampling
//    factor.  The number of points in the sinc is based on the desired number
//    of crossings and the sampling rate.
// 2. The calculated sinc function is stored in the instance variable, _sincFunction.
// ############################# Private Function ###############################
void SincInterpolator::calculateSincFunction()
{
  int           half_taps, i;
  float         f_s, f_c, f_s_over, first_cross;
  double        fo;
  const         double *b_coeff;
//
// Calculate the input (over) sampling rate:
//
  f_s                   = 1.;
  if(_inputTimeStep > 0.)
    f_s                 = 1./_inputTimeStep;
  f_s_over              = f_s*_oversamplingFactor;
//
// Now calculate the BW of the ideal LPF (-> sinc in time domain):
// And the number of filter taps based on this bandwidth
//
  f_c                   = (1. - _transitionWidth)*f_s/2.;
  first_cross           = 0.5/f_c;
  half_taps             = ROUND((_numberCrossings*first_cross*f_s_over));
  _numberFilterTaps     = 2*half_taps+1;                        // Make it odd
//
// Now, use the FIR filter to calculate sinc function (= filter coefficients)
//
  fo                    = 0.;
  _firFilter->setFIRFilterType(FIR_WINDOW);
  _firFilter->setNumberTaps(_numberFilterTaps);
  _firFilter->setFilterParameters(fo, f_c, f_s_over);
  _firFilter->findTransferFunction();
//
// Get the coefficients from the filter, and store into _sincFunction
//
  delete [] _sincFunction;
   b_coeff              = _firFilter->bCoeffs();
  _sincFunction         = new double [_numberFilterTaps];
  for(i=0; i<_numberFilterTaps; i++)
    _sincFunction[i]    = b_coeff[i];
  return;
}

// ############################# Private Method ###############################
// interpolateDouble -- Interpolates an input set of points using a sinc function
//                      centered at the output sample point.  The sinc function is
//                      linearly interpreted to get its values at the input sample
//                      times.
//
// Input:       p:      the fractional step of the ouput point from the next lower
//                      input sample          
// Output:              The interpolated output
//
// Notes:
// ############################# Public Method ###############################
double SincInterpolator::interpolateDouble(double p)
{
  int           i, n, p_int;
  double        eta, p_norm, one_minus_eta, sinc_interp;
  double        sinc_sum, output_sum;
  const         double *signal_sample;

//
// Set the ring buffer to point to oldest signal sample
//
  _inputBuffer->setReadPointerToWritePointer();
//
// Loop over input samples, # of input samples and filter taps should be odd:
//
  p_norm                = _oversamplingFactor*(1.-p);
  p_int                 = (int)(p_norm);                // index of filter for first point
  eta                   = p_norm - p_int;               // filter tap interpolation factor
  one_minus_eta         = 1. - eta;
  n                     = p_int;
  sinc_sum              = 0.;
  output_sum            = 0.;
  for(i=0; i<_numberBufferedSamples; i++)
  {
    sinc_interp         = _sincFunction[n]*(one_minus_eta) + _sincFunction[n+1]*eta;
    n                   += _oversamplingFactor;
    signal_sample       = _inputBuffer->readFromBuffer();
    sinc_sum            += sinc_interp;
    output_sum          += sinc_interp* (*signal_sample);
    if( (n+3) > _numberFilterTaps)                      // Error check this shouldn't happen
    {
      printf("Array index out of range in SincInterpolator::interpolateDouble()\n");
      printf("_inputSample = %d, _outputSample = %d, _firstOutput = %d, p_int = %d, p = %g\n",
             _inputSample, _outputSample, _firstOutputSample, p_int, p);
      printf("_buffered = %d, filter taps = %d_inputStep = %g, _outputStep = %g\n", 
             _numberBufferedSamples, _numberFilterTaps, _inputTimeStep, _outputTimeStep);
      break;
    }
  }
  if(sinc_sum != 0.)
    output_sum          /= sinc_sum;
  return output_sum;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the SincInterpolator class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
SincInterpolator::SincInterpolator(int numberCrossings)
{
  _inputBuffer          = NULL;
  _sincFunction         = NULL;
  _firFilter            = new FIRFilter();

  setInputTimeStep(DEFAULT_TIME_STEP);
  setOutputTimeStep(DEFAULT_TIME_STEP);
  setTransitionWidth(DEFAULT_TRANS_WIDTH);
  setNumberCrossings(numberCrossings);
  setOversamplingFactor(DEFAULT_OVER_SAMPLING);
  setWindowType(DEFAULT_WINDOW);
  reset();

  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the SincInterpolator class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
SincInterpolator::~SincInterpolator()
{
  delete [] _sincFunction;
  delete _firFilter;
  delete _inputBuffer;
  return;
}

// ############################# Public Method ###############################
// reset -- Resets the interpolator for real time processing
//
// Input:                       None
//          
// Output:                      None
//
// NOTES:
// 1. This function must be called before performing interpolations, since
//    it calls calculateSincFunction().
// ############################# Public Method ###############################
void SincInterpolator::reset()
{
  int           half_input;
  double        output_time;

  _inputSample  = _outputSample = 0;
  calculateSincFunction();
//
// Calculate the number of input points needed to buffer:
//
  delete _inputBuffer;
  _numberBufferedSamples        = _numberFilterTaps/_oversamplingFactor - 2;
  _numberBufferedSamples        = MAX(1, _numberBufferedSamples);
  _inputBuffer                  = new DoubleBuffer(_numberBufferedSamples);
  _inputBuffer->resetBuffer();
//
//  Also, calculate the first non-zero output sample:
//
  half_input                    = _numberBufferedSamples/2;
  output_time                   = _inputTimeStep*half_input;
  _firstOutputSample            = 0;
  if(_outputTimeStep > 0.)
    _firstOutputSample          = (int)(output_time/_outputTimeStep);
  if( (_firstOutputSample*_outputTimeStep) < output_time)
    _firstOutputSample          += 2;                           // Non-coincident input-output sampling
  else
    _firstOutputSample++;                                       // Coincident input-output sampling
  return;
}

// ############################# Public Method ###############################
// setOversamplingFactor -- Sets a new oversampling factor for the sinc function
//
// Input:       factor:         New oversampling factor
//          
// Output:                      None
//
// NOTES:
// 1. The oversampling factor determines the amount of oversampling of the sinc
//    function.  When using the sinc interpolation, the values of the sinc
//    function are determined by linear interpolation, therefore, an oversampling
//    factor of at least 10 needs to be used for adequate performance
// ############################# Public Method ###############################
void SincInterpolator::setOversamplingFactor(int factor)
{
  _oversamplingFactor   = MAX(1, factor);
  return;
}

// ############################# Public Method ###############################
// setNumberCrossings -- Sets a new number of sinc crossings
//
// Input:       crossings:      Number of one sided sinc crossings
//          
// Output:                      None
// ############################# Public Method ###############################
void SincInterpolator::setNumberCrossings(int crossings)
{
  _numberCrossings      = MAX(1, crossings);
  return;
}

// ############################# Public Method ###############################
// setWindowType -- Sets a new window function for the sinc calculation
//
// Input:       type:           New type of window (see DataWindow.h)
//          
// Output:                      None
//
// ############################# Public Method ###############################
void SincInterpolator::setWindowType(int type)
{
  _firFilter->setWindowType(type);
  _windowType           = _firFilter->windowType();
  return;
}

// ############################# Public Method ###############################
// setInputSampleTime -- Sets a new input sampling time for real time processing
//
// Input:       time:           New input sample time in seconds
//          
// Output:                      None
// ############################# Public Method ###############################
void SincInterpolator::setInputTimeStep(float time)
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
void SincInterpolator::setOutputTimeStep(float time)
{
  _outputTimeStep       = time;
  return;
}

// ############################# Public Method ###############################
// setTransitionWidth -- Sets a new transition width as a fraction
//
// Input:       width:          Transition width as a fraction of fs/2
//          
// Output:                      None
// ############################# Public Method ###############################
void SincInterpolator::setTransitionWidth(float width)
{
  _transitionWidth      = MAX(MIN_TRANSITION, width);
  _transitionWidth      = MIN(MAX_TRANSITION, _transitionWidth);
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
int SincInterpolator::numberOfInputPointsNeeded()
{
  int           input_samples;
  double        output_time;
//
// For the boundary condition, we set the number of
// points needed equal to the buffer size:
//
  if(_outputSample == 0)
    return _numberBufferedSamples;
//
// The number of input samples depends on the output time, and the filter
// length:
//
  output_time   = _outputSample * _outputTimeStep;
//
// Calculate total number of input samples needed.  Since the output sample
// is centered on the sinc function, the number of input samples is equal to
// ((output_time/_inputTimeStep) + _bufferSize/2, where the plus 1 is due
// to the number samples = number sections + 1
//
  input_samples = (int)(output_time/_inputTimeStep) + _numberBufferedSamples/2;
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
void    SincInterpolator::bufferNextInput(double inputData)
{
  int   i;
  _inputSample++;
  _inputBuffer->writeToBuffer(&inputData);
  if(_inputSample==1)                           // Fill up the buffer if first sample
    for(i=0; i<_numberFilterTaps-1; i++)
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
double SincInterpolator::getNextOutput()
{
  int           half_input;
  double        p;
  double        input_time, output_time;
//
// We first need to find, p, the percentage step of the ouput point from
// the next lower input sample
//
  half_input    = _numberBufferedSamples/2;
  output_time   = _outputSample * _outputTimeStep;
  input_time    = (_inputSample - half_input)*_inputTimeStep;
  p             = (output_time - input_time)/_inputTimeStep;
  p             = MAX(0., p);
  p             = MIN(p, 1.);                                   // Error check
//
// Update output counter:
//
  _outputSample++;
  if(_outputSample < _firstOutputSample)
    return 0.;
  else
  {
//
// Return the interpolated value:
//
    return interpolateDouble(p);
  }
}
