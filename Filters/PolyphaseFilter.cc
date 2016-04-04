/************************************************************************
 *                                                                      *
 * This subclass of object implements an abstract polyphase filter.     *
 * For actual implementations, see the subclasses, PolyphaseDecimator,  *
 * and PolyphaseInterpolator.                                           *
 * Based on Crochiere and Rabiner, pp. 79-85.                           *
 *                                                                      *
 * File:PolyphaseFilter.cc                                              *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 04/02/99  - Started.                                             *
 ************************************************************************/

#include "PolyphaseFilter.h"                                    // Object prototypes
#include "AbstractFIR.h"
#include "FIRFilter.h"
#include "DataWindow.h"
#include <math.h>
#include <stdio.h>

#define MAX(a,b)    ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)    ( ((a)>(b)) ? (b) : (a) )

// ############################# Private Method ###############################
// findFilterCoefficients -- Finds a new set of filter coefficients for the
//                           polyphase filters.
//
// Input:                       None
//          
// Output:                      None
// ############################# Private Method ###############################
void PolyphaseFilter::findFilterCoefficients()
{
  int           i, j, k, proto_taps;
  double        cutoff, center;
  const         double *proto_coeffs;
  double        *poly_coeffs;
//
// Allocate the prototype filter:
//
  cutoff                = _cutoffFactor*_samplingFrequency/_numberPhases/2.;
  center                = 0.;                           // low pass filter
  proto_taps            = _tapsPerPhase*_numberPhases;
  if(_prototypeFilter == NULL)
    _prototypeFilter    = new FIRFilter(LOW_PASS, center, cutoff, _samplingFrequency, proto_taps-1);
  else
  {
    if(!_useLoadedCoefficients)
    {
      _prototypeFilter->setNumberTaps(proto_taps);
      _prototypeFilter->setFilterParameters(center, cutoff, _samplingFrequency);
    }
  }
  if(!_useLoadedCoefficients)
  {
    _prototypeFilter->setFIRFilterType(FIR_WINDOW);
    _prototypeFilter->setWindowType(BLACKMAN_HARRIS4);
    _prototypeFilter->findTransferFunction();
  }
//
// Set the polyphase coefficients for each phase/branch:
//
  poly_coeffs           = new double[_tapsPerPhase];
  proto_coeffs          = _prototypeFilter->bCoeffs();
  for(i=0; i<_numberPhases; i++)
  {
    k   = i;
    for(j=0; j<_tapsPerPhase; j++)
    {
      poly_coeffs[j]    = proto_coeffs[k];
      k                 += _numberPhases;
    }
    _polyphaseFilters[i]->setFilterBCoeffs(poly_coeffs);
  }
  
  delete [] poly_coeffs;
    
  return;
}


// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the PolyphaseFilter class.
//
// Input:       numberPhases:   Number of polyphase branches
//              tapsPerPhase:   Size of filter for each phase
//              cutoffFactor:   Cutoff factor for each polyphase filter for dec/interp
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
PolyphaseFilter::PolyphaseFilter(int numberPhases, int tapsPerPhase, float cutoffFactor)
{
  _polyphaseFilters             = NULL;
  _prototypeFilter              = NULL;
  _filterOutput                 = NULL;
  _numberPhases                 = 0;
  _cutoffFactor                 = cutoffFactor;
  _tapsPerPhase                 = tapsPerPhase;
  _samplingFrequency            = 1.;
  _commutatorLocation           = 0;
  _useLoadedCoefficients        = NO;
  _outputSamples                = _oldOutputSamples             = 0;
  setNumberPhases(numberPhases);

  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the PolyphaseFilter class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
PolyphaseFilter::~PolyphaseFilter()
{
  int   i;
  for(i=0; i<_numberPhases; i++)
    delete _polyphaseFilters[i];
  delete [] _polyphaseFilters;
  delete    _prototypeFilter;
  delete [] _filterOutput;
  return;
}

// ############################# Public Method ###############################
// zeroOutTaps -- Resets the polyphase filter
//
// Input:                       None
//          
// Output:                      None
// ############################# Public Method ###############################
void PolyphaseFilter::zeroOutTaps()
{
  int   i;
//
// Reset the polyphase branch filters:
//
  for(i=0; i<_numberPhases; i++)
    _polyphaseFilters[i]->zeroOutTaps();
  _commutatorLocation   = 0;
  return;
}

// ############################# Public Method ###############################
// setPrototypeCoeffs -- Loads a set of prototype coefficients.
//
// Input:       coeffs:         filter tap weights for prototype filter
//              size:           number of coefficients
//
// Output:                      None
// ############################# Public Method ###############################
void PolyphaseFilter::setPrototypeCoeffs(double *coeffs, int size)
{
  if(size != (_numberPhases*_tapsPerPhase))             // Check for error
  {
    printf("Number coeffs = %d != (_numberPhases*_tapsPerPhase) = %d in setPrototypeCoeffs\n",
           size, _numberPhases*_tapsPerPhase);
    return;
  }
  _useLoadedCoefficients        = YES;
  _prototypeFilter->setNumberTaps(size);
  _prototypeFilter->setFilterBCoeffs(coeffs);
  return;
}

// ############################# Public Method ###############################
// setNumberPhases -- Sets a new number of phases/number of filter branches
//
// Input:       phases:         number of phases/polyphase filters
//          
// Output:                      None
// ############################# Public Method ###############################
void PolyphaseFilter::setNumberPhases(int phases)
{
  int   i;
  float cutoff, center;
//
// Delete old filters:
//
  for(i=0; i<_numberPhases; i++)
    delete _polyphaseFilters[i];
  delete  [] _polyphaseFilters;
//
// Allocate new filters:
//
  _numberPhases                 = MAX(0, phases);
  _numberPhases                 = MIN(_numberPhases, MAX_PHASES);
  _polyphaseFilters             = new AbstractFIR* [_numberPhases];
  cutoff                        = _samplingFrequency/4.;        // this isn't used
  center                        = 0.;                           // also not used
  for(i=0; i<_numberPhases; i++)
    _polyphaseFilters[i]        = new AbstractFIR(LOW_PASS, center, cutoff, _samplingFrequency, _tapsPerPhase-1);
//
// Find the filter coefficients:
//
  findFilterCoefficients();
  return;
}

// ############################# Public Method ###############################
// setTapsPerPhase -- Sets a new number of taps per each polyphase branch
//
// Input:       taps:           number of taps for each polyphase filter
//          
// Output:                      None
// ############################# Public Method ###############################
void PolyphaseFilter::setTapsPerPhase(int taps)
{
  int   i;
//
// Set filter taps:
//
  _tapsPerPhase                 = taps;
  for(i=0; i<_numberPhases; i++)
    _polyphaseFilters[i]->setNumberTaps(_tapsPerPhase);
//
// Find the filter coefficients:
//
  findFilterCoefficients();

  return;
}

// ############################# Public Method ###############################
// setCutoffFactor -- Sets a new cutoff factor
//
// Input:       factor:         cutoff factor
//          
// Output:                      None
//
// Notes:
// 1. The cutoff factor is used to set the low pass filter cutoff frequency.
//    The cutoff frequency is equal to cutoffFactor*samplingFrequency/numberPhases.
//    Where the number of phases will also be equal to the decimation/interpolation
//    factor.
// ############################# Public Method ###############################
void PolyphaseFilter::setCutoffFactor(float factor)
{
//
// Set filter cutoff:
//
  _cutoffFactor                 = factor;
//
// Find the filter coefficients with new cutoff factor
//
  findFilterCoefficients();

  return;
}


