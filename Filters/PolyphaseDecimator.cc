/************************************************************************
 *                                                                      *
 * This subclass of PolyphaseFilter implements a polyphase filter for   *
 * decimating an input signal.                                          *
 * Based on Crochiere and Rabiner, pp. 79-85.                           *
 *                                                                      *
 * File:PolyphaseDecimator.cc                                           *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 04/05/02  - Started.                                             *
 ************************************************************************/

#include "PolyphaseDecimator.h"                                 // Object prototypes
#include "AbstractFIR.h"
#include <math.h>
#include <stdio.h>

#define MAX(a,b)    ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)    ( ((a)>(b)) ? (b) : (a) )

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the PolyphaseDecimator class.
//
// Input:       numberPhases:   Number of polyphase branches
//              tapsPerPhase:   Size of filter for each phase
//              cutoffFactor:   Cutoff factor for each polyphase filter for dec/interp
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
PolyphaseDecimator::PolyphaseDecimator(int numberPhases, int tapsPerPhase, float cutoffFactor)
                        :PolyphaseFilter(numberPhases, tapsPerPhase, cutoffFactor)
{
  setNumberPhases(numberPhases);

  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the PolyphaseDecimator class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
PolyphaseDecimator::~PolyphaseDecimator()
{
 return;
}

// ############################# Public Method ###############################
// zeroOutTaps -- Resets the polyphase filter
//
// Input:                       None
//          
// Output:                      None
// ############################# Public Method ###############################
void PolyphaseDecimator::zeroOutTaps()
{
//
// Perform super's function:
//
  PolyphaseFilter::zeroOutTaps();
//
// Add this:
//
  _outputSumReal        = _outputSumImag        = 0.;
  _outputSumComplex     = Complex(0., 0.);
  return;
}

// ############################# Public Method ###############################
// setNumberPhases -- Sets a new number of phases/number of filter branches
//
// Input:       phases:         number of phases/polyphase filters
//          
// Output:                      None
// ############################# Public Method ###############################
void PolyphaseDecimator::setNumberPhases(int phases)
{
//
// Perform super's function:
//
  PolyphaseFilter::setNumberPhases(phases);
//
// Add this:
//
  _decimationRatio      = _numberPhases;
  return;
}

// ############################# Public Method ###############################
// loadInputReal -- Loads a new floating point number into the decimator
//
// Input:       inputSample:    new input sample to filter/decimate
//          
// Output:                      None
//
// Notes:
// 1. Call readyForOutput(), and getFloatOutput() to access output data.
// ############################# Public Method ###############################
void PolyphaseDecimator::loadInputReal(float inputSample)
{
  double        filter_output;
//
// Update the commutator:
//
  _commutatorLocation--;
  if(_commutatorLocation < 0)
    _commutatorLocation = _decimationRatio - 1;
//
// Filter the new input using the appropriate polyphase branch:
//
  filter_output         = _polyphaseFilters[_commutatorLocation]->filterFloatData(inputSample);
  _outputSumReal        += filter_output;
  return;
}

// ############################# Public Method ###############################
// loadInputRealImag -- Loads a new set of real and imaginary data into filter
//
// Input:       inputReal:      real part of new input sample to filter/decimate
//              inputImag:      imag part of new input sample to filter/decimate
//          
// Output:                      None
//
// Notes:
// 1. Call readyForOutput(), and getFloatOutput() to access output data.
// ############################# Public Method ###############################
void PolyphaseDecimator::loadInputRealImag(float inputReal, float inputImag)
{
  Complex       filter_input, filter_output;
//
// Update the commutator:
//
  _commutatorLocation--;
  if(_commutatorLocation < 0)
    _commutatorLocation = _decimationRatio - 1;
//
// Filter the new input using the appropriate polyphase branch:
//
  filter_input          = Complex(inputReal, inputImag);
  filter_output         = _polyphaseFilters[_commutatorLocation]->filterComplexData(filter_input);
  _outputSumReal        += filter_output.real();
  _outputSumImag        += filter_output.imag();
  return;
}

// ############################# Public Method ###############################
// loadInputComplex -- Loads a new complex number into the decimator
//
// Input:       inputSample:    new input sample to filter/decimate
//          
// Output:                      None
//
// Notes:
// 1. Call readyForOutput(), and getComplexOutput() to access output data.
// ############################# Public Method ###############################
void PolyphaseDecimator::loadInputComplex(Complex &inputSample)
{
  Complex       filter_output;
//
// Update the commutator:
//
  _commutatorLocation--;
  if(_commutatorLocation < 0)
    _commutatorLocation = _decimationRatio - 1;
//
// Filter the new input using the appropriate polyphase branch:
//
  filter_output         = _polyphaseFilters[_commutatorLocation]->filterComplexData(inputSample);
  _outputSumComplex     += filter_output;
  return;
}

// ############################# Public Method ###############################
// readyForOutput -- Checks whether the decimator is ready to output data
//
// Input:                       None
//          
// Output:                      None
//
// Notes:
// 1. If this function returns YES, call getFloatOutput or getComplexOutput()
//    to access output data.
// ############################# Public Method ###############################
LOGICAL PolyphaseDecimator::readyForOutput()
{
//
// The following assumes data is input to branch M-1, M-2, ..., 1, 0
  return (_commutatorLocation == 0);
}

// ############################# Public Method ###############################
// getRealOutput -- Returns the decimated floating point output after M inputs
//
// Input:                       None
//          
// Output:                      Filtered/decimated output
//
// Notes:
// ############################# Public Method ###############################
float PolyphaseDecimator::getRealOutput()
{
  float output;

  output                = _outputSumReal;
  _outputSumReal        = 0.;                   // reset accumulator
  return        output;
}

// ############################# Public Method ###############################
// getImagOutput -- Returns the imag part of the decimated floating point
//                  output after M inputs
//
// Input:                       None
//          
// Output:                      Filtered/decimated output
//
// Notes:
// ############################# Public Method ###############################
float PolyphaseDecimator::getImagOutput()
{
  float output;

  output                = _outputSumImag;
  _outputSumImag        = 0.;                   // reset accumulator
  return        output;
}

// ############################# Public Method ###############################
// getComplexOutput -- Returns the decimated complex output after M inputs
//
// Input:                       None
//          
// Output:                      Filtered/decimated output
//
// Notes:
// ############################# Public Method ###############################
Complex &PolyphaseDecimator::getComplexOutput()
{

  _complexOutput        = _outputSumComplex;
  _outputSumComplex     = Complex(0.,0.);                       // reset accumulator
  return        _complexOutput;
}

// ############################# Public Method ###############################
// filterFloatArray -- interpolates an input array, and returns pointer to output
//
// Input:       inputData:      Array of input data
//              samples:        size of above array
//          
// Output:                      Filtered/interpolated output
//
// Notes:
// 1. Call outputSamples() to get size of output array
// ############################# Public Method ###############################
float *PolyphaseDecimator::filterFloatArray(const float *inputData, int samples)
{
  int   i, j, k;

//
// Allocate output array:
//
  _outputSamples        = samples/_decimationRatio;
  if(_outputSamples > _oldOutputSamples)
  {
    _oldOutputSamples   = _outputSamples;
    delete [] _filterOutput;
    _filterOutput       = new float[_outputSamples];
  }
//
// FIlter the input data:
//
  k     = 0;
  for(i=0; i<_outputSamples; i++)
  {
    for(j=0; j<_decimationRatio; j++)
    {
      loadInputReal(inputData[k]);
      k++;
    }
    _filterOutput[i]    = getRealOutput();
  }
  return _filterOutput;
}