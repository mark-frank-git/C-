#ifndef _POLYPHASE_FILTER_H
#define _POLYPHASE_FILTER_H 1
/************************************************************************
 *                                                                      *
 * This subclass of object implements an abstract polyphase filter.     *
 * For actual implementations, see the subclasses, PolyphaseDecimator,  *
 * and PolyphaseInterpolator.                                           *
 * Based on Crochiere and Rabiner, pp. 79-85.                           *
 *                                                                      *
 * File:PolyphaseFilter.h                                               *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 04/02/99  - Started.                                             *
 ************************************************************************/

#define DEFAULT_CUTOFF_FACTOR   0.95
#define MAX_PHASES              100

#ifndef LOGICAL
#define LOGICAL                 char
#endif
#ifndef YES
#define YES                     1
#define NO                      0
#endif

class   FIRFilter;                                      // Class prototypes
class   AbstractFIR;

class PolyphaseFilter
{
protected:
  int           _numberPhases;                  // This will also be equal to the
                                                // decimate or interpolate factor and number of filters
  int           _tapsPerPhase;                  // # of taps for each of the polyphase filters
  int           _commutatorLocation;            // Location of commutator switch
  int           _outputSamples;                 // # of samples in output array
  int           _oldOutputSamples;              // Old # of output array samples
  
  float         _cutoffFactor;                  // Factor to select cutoff frequency
  float         _samplingFrequency;             // Always set to 1

  float         *_filterOutput;                 // output filtered array

  FIRFilter     *_prototypeFilter;              // The prototype filter for finding polyphase branches
  AbstractFIR   **_polyphaseFilters;            // The filters for each phase
  LOGICAL       _useLoadedCoefficients;         // If yes, use pre-calculated coefficients

//
// Private methods:
//
  void          findFilterCoefficients();       // Find the coefficients for the filters
//
// Public methods:
//
public:
/********************************
 * Constructors, destructor     *
 *******************************/
  PolyphaseFilter(int numberPhases, int tapsPerPhase, float cutoffFactor=DEFAULT_CUTOFF_FACTOR);
  virtual  ~PolyphaseFilter();

/********************************
 * Resetting the filter         *
 ********************************/
  virtual void  zeroOutTaps();
  
/**************************
 * Setting parameters:    *
 **************************/
  virtual void  setPrototypeCoeffs(double *b, int size);
  virtual void  setNumberPhases(int phases);
  virtual void  setTapsPerPhase(int taps);
  void          setCutoffFactor(float factor);

/**********************
 * Get parameters:    *
 **********************/
  int           numberPhases()          {return _numberPhases;}
  int           tapsPerPhase()          {return _tapsPerPhase;}
  int           outputSamples()         {return _outputSamples;}
  float         cutoffFactor()          {return _cutoffFactor;}

};

#endif
