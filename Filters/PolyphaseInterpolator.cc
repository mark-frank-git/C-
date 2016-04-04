/************************************************************************
 *                                                                      *
 * This subclass of PolyphaseFilter implements a polyphase filter for   *
 * decimating an input signal.                                          *
 * Based on Crochiere and Rabiner, pp. 79-85.                           *
 *                                                                      *
 * File:PolyphaseInterpolator.cc                                        *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 04/05/99  - Started.                                             *
 ************************************************************************/

#include "PolyphaseInterpolator.h"                                      // Object prototypes
#include "AbstractFIR.h"
#include <math.h>
#include <stdio.h>

#define MAX(a,b)    ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)    ( ((a)>(b)) ? (b) : (a) )

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the PolyphaseInterpolator class.
//
// Input:       numberPhases:   Number of polyphase branches
//              tapsPerPhase:   Size of filter for each phase
//              cutoffFactor:   Cutoff factor for each polyphase filter for dec/interp
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
PolyphaseInterpolator::PolyphaseInterpolator(int numberPhases, int tapsPerPhase, float cutoffFactor)
                        :PolyphaseFilter(numberPhases, tapsPerPhase, cutoffFactor)
{
  setNumberPhases(numberPhases);
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the PolyphaseInterpolator class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
PolyphaseInterpolator::~PolyphaseInterpolator()
{
 return;
}

// ############################# Public Method ###############################
// setNumberPhases -- Sets a new number of phases/number of filter branches
//
// Input:       phases:         number of phases/polyphase filters
//          
// Output:                      None
// ############################# Public Method ###############################
void PolyphaseInterpolator::setNumberPhases(int phases)
{
//
// Perform super's function:
//
  PolyphaseFilter::setNumberPhases(phases);
//
// Add this:
//
  _interpolationRatio   = _numberPhases;
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
// 1. Call readyForOutput(), and getRealOutput() to access output data.
// ############################# Public Method ###############################
void PolyphaseInterpolator::loadInputReal(float inputSample)
{
//
// Store the input data for processing by each of the polyphase branches:
//
  _inputRealData        = inputSample;
  return;
}

// ############################# Public Method ###############################
// loadInputFloat -- Loads a new floating point number into the decimator
//
// Input:       inputReal:      real part input sample to filter/decimate
//              inputImag:      imag part of input
//          
// Output:                      None
//
// Notes:
// 1. Call readyForOutput(), and getComplexOutput() to access output data.
// ############################# Public Method ###############################
void PolyphaseInterpolator::loadInputRealImag(float inputReal, float inputImag)
{
//
// Store the input data for processing by each of the polyphase branches:
//
  _inputComplexData     = Complex(inputReal, inputImag);
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
void PolyphaseInterpolator::loadInputComplex(Complex &inputSample)
{
//
// Store the input data for processing by each of the polyphase branches:
//
  _inputComplexData     = inputSample;
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
LOGICAL PolyphaseInterpolator::readyForOutput()
{
//
// The following assumes data is output from branches 0, 1, 2, ..., L-1
  return (_commutatorLocation < _interpolationRatio);
}

// ############################# Public Method ###############################
// getFloatOutput -- Returns the decimated floating point output after M inputs
//
// Input:                       None
//          
// Output:                      Filtered/decimated output
//
// Notes:
// ############################# Public Method ###############################
float PolyphaseInterpolator::getRealOutput()
{
  float output;

  if(_commutatorLocation == _interpolationRatio)
    _commutatorLocation = 0;
  output                = _polyphaseFilters[_commutatorLocation]->filterFloatData(_inputRealData);
  _commutatorLocation++;
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
Complex &PolyphaseInterpolator::getComplexOutput()
{

  if(_commutatorLocation == _interpolationRatio)
    _commutatorLocation = 0;
  _complexOutput        = _polyphaseFilters[_commutatorLocation]->filterComplexData(_inputComplexData);
  _commutatorLocation++;
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
float *PolyphaseInterpolator::filterFloatArray(const float *inputData, int samples)
{
  int   i, j, k;

//
// Allocate output array:
//
  _outputSamples        = samples*_interpolationRatio;
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
  for(i=0; i<samples; i++)
  {
    loadInputReal(inputData[i]);
    for(j=0; j<_interpolationRatio; j++)
    {
      _filterOutput[k]  = getRealOutput();
      k++;
    }
  }
  return _filterOutput;
}