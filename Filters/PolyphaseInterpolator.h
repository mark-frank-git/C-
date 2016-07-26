#ifndef _POLYPHASE_INTERPOLATOR_H
#define _POLYPHASE_INTERPOLATOR_H 1
/************************************************************************
 *                                                                      *
 * This subclass of PolyphaseFilter implements a polyphase filter for   *
 * interpolating an input signal.                                       *
 * Based on Crochiere and Rabiner, pp. 79-85.                           *
 *                                                                      *
 * File:PolyphaseInterpolator.h                                         *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 04/05/02  - Started.                                             *
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

class PolyphaseInterpolator : public PolyphaseFilter
{
protected:
  int           _interpolationRatio;            // Same as _numberPhases

  float         _inputRealData;                 // input floating point data
  double        _inputImagData;                 // input floating point imag data

  Complex       _inputComplexData;              // input complex data
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
  PolyphaseInterpolator(int numberPhases, int tapsPerPhase, float cutoffFactor=DEFAULT_CUTOFF_FACTOR);
  virtual  ~PolyphaseInterpolator();

/**************************
 * Setting parameters:    *
 **************************/
  virtual void  setNumberPhases(int phases);

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
  Complex       &getComplexOutput();
  
/************************
 * Filtering arrays:    *
 ************************/
  float         *filterFloatArray(const float *inputData, int samples);

};

#endif
