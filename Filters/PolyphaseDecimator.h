#ifndef _POLYPHASE_DECIMATOR_H
#define _POLYPHASE_DECIMATOR_H 1
/************************************************************************
 *                                                                      *
 * This subclass of PolyphaseFilter implements a polyphase filter for   *
 * decimating an input signal.                                          *
 * Based on Crochiere and Rabiner, pp. 79-85.                           *
 *                                                                      *
 * File:PolyphaseDecimator.h                                            *
 *                                                                      *
 ************************************************************************/
#include        "PolyphaseFilter.h"
#include        <GNU/Complex.h>


#ifndef LOGICAL
#define LOGICAL                 char
#endif
#ifndef YES
#define YES                     1
#define NO                      0
#endif

class                           Complex;        // Class prototype

class PolyphaseDecimator : public PolyphaseFilter
{
protected:
  int           _decimationRatio;               // Same as _numberPhases

  double        _outputSumReal;                 // output floating point real data
  double        _outputSumImag;                 // output floating point imag data

  Complex       _outputSumComplex;              // complex output accumulator
  Complex       _complexOutput;                 // complex output data
//
// Private methods:
//
//
// Public methods:
//
public:
/********************************
 * Constructors, destructor     *
 *******************************/
  PolyphaseDecimator(int numberPhases, int tapsPerPhase, float cutoffFactor=DEFAULT_CUTOFF_FACTOR);
  virtual  ~PolyphaseDecimator();

/********************************
 * Resetting the filter         *
 ********************************/
  virtual void  zeroOutTaps();

/**************************
 * Setting parameters:    *
 **************************/
  void  setNumberPhases(int phases);

/************************
 * Get parameters:      *
 ***********************/

/************************
 * Inputting new data:  *
 ************************/
  void  loadInputReal(float inputSample);
  void  loadInputRealImag(float inputReal, float inputImag);
  void  loadInputComplex(Complex &inputSample);

/********************************
 * Checking for output data     *
 ********************************/
  LOGICAL       readyForOutput();
  
/************************
 * Getting output data: *
 ************************/
  float         getRealOutput();
  float         getImagOutput();
  Complex       &getComplexOutput();

/************************
 * Filtering arrays:    *
 ************************/
  float         *filterFloatArray(const float *inputData, int samples);
  
};

#endif
