#ifndef _POLYPHASE_DELAY_H
#define _POLYPHASE_DELAY_H      1
/************************************************************************
 *                                                                      *
 * This subclass of PolyphaseInterpolator implements a polyphase filter *
 * for delaying an input signal.                                        *
 * Based on Crochiere and Rabiner, pp. 79-85.                           *
 *                                                                      *
 * File:PolyphaseDelay.h                                                *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 04/05/99  - Started.                                             *
 ************************************************************************/
#include        "PolyphaseInterpolator.h"
#include        <GNU/Complex.h>


#ifndef LOGICAL
#define LOGICAL                 char
#endif
#ifndef YES
#define YES                     1
#define NO                      0
#endif


class PolyphaseDelay : public PolyphaseInterpolator
{
protected:
  int           _selectedPolyphase;                     // Filter selected on delay
  float         _normalizedDelay;                       // Delay relative to 1 sampling period
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
  PolyphaseDelay(int numberPhases, int tapsPerPhase, float cutoffFactor=DEFAULT_CUTOFF_FACTOR, float delay=0.0);
  virtual  ~PolyphaseDelay();

/**************************
 * Setting parameters:    *
 **************************/
  void  setNumberPhases(int phases);
  void  setTapsPerPhase(int taps);
  void  setDelay(float delay);

/************************
 * Get parameters:      *
 ***********************/
  float delay()                                 {return _normalizedDelay;}

/****************************************
 * Filtering float and Complex data     *
 * assuming a polyphase structure.      *
 ****************************************/
  LOGICAL       filterFloatArray(float * input, int numberPts); // Filter float array
  double        filterFloatData(float input);                   // Filter float data
  double        filterDoubleData(double input);                 // Filter double data
  LOGICAL       filterComplexArray(Complex *x, int numberPts);  // Filter complex array
  Complex       filterComplexData(Complex &input);              // Filter complex data

};

#endif
