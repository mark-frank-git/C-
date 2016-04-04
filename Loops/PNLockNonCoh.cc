/************************************************************************
 *                                                                      *
 * This subclass of PNLock supports locking to a PN code, assuming      *
 * non-coherent operation (i.e., carrier frequency offset).             *
 *                                                                      *
 * File:PNLockNonCoh.cc                                                 *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 03/17/99  - Subclassed from PNLock.                              *
 ************************************************************************/

#include "PNLockNonCoh.h"                                       // Object prototypes
#include <Filters/IIRFilter.h>
#include <C_Libraries/constants.h>
#include <Generators/PNFiltered.h>
#include <stdio.h>
#include <math.h>

#ifndef EPS
#define EPS     1.e-5
#endif

#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define BIT_MAP(a)      ( ((a)==0 ) ? 1. : -1. )                // (0,1) -> (1, -1)
#define ROUND(a)        ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )
#define SGN(a)          ( ((a)>=0.) ? 1. : -1. )


// ############################# Private Function ###############################
// loopLocked -- Checks to see if the loop has locked
//
// Input:       input:          next complex input
//
// Output:                      YES if lock detector > threshold
//
// Notes:
// ############################# Private Function ###############################
LOGICAL PNLockNonCoh::loopLocked(Complex &input)
{
  float         despread_i;
  float         filter_out;
//
// Check for presence of pilot:
//
  despread_i                            = input.real()*_iPNOnTime;
  filter_out                            = _onTimeFilter->filterFloatData(despread_i);
  _lockDetectorSum                      += filter_out*filter_out;               // energy detector for non-coherent
//
// Update lock statistics:
//
  if(++_lockDetectCounter >= _lockDetectorCount)
  {
    _lockDetectorSum                    /= _lockDetectorCount;
    if(_lockDetectorSum >= _lockThreshold)
      return YES;
    _lockDetectorSum                    = 0.;
    _lockDetectCounter                  = 0;
  }
  return NO;
}

// ############################# Private Function ###############################
// findPDOutput -- Processes a new complex input through the phase detector
//
// Input:       input:          next complex input
//
// Output:                      None
//
// Notes:
// 1. The instance variable, _pdOutput, is modified.
// 2. This function implements the PD for a data directed PLL.
// 3. For the data directed PLL, we assume the complex down conversion is done
//    prior to this block, with the sliced data taken after the down conversion
// ############################# Private Function ###############################
void PNLockNonCoh::findPDOutput(Complex &input)
{
  Complex       filter_out;
//
// First, check for loop locked (this call should probably be located elsewhere):
//
  if(!_loopLocked  && (_lockDetectorType != NO_LOCK_DETECTOR) )
  {
    _loopLocked         = loopLocked(input);
    if(_loopLocked)
      initializeLoopFilterFor(_loopFilterType, _inputAmplitude, YES);   // Set new loop bandwidth
  }
      
//
// Get output of early and late arms after quadratic detector:
//
  switch(_phaseDetectorType)
  {
      case COMPLEX_PHASE_DETECTOR:
      default:
        _iEarlyOut      = _iPNEarly*input.real()        + _qPNEarly*input.imag();       // input * conj(pn)
        _qEarlyOut      = _iPNEarly*input.imag()        - _qPNEarly*input.real();
        _iLateOut       = _iPNLate*input.real()         + _qPNLate*input.imag();        // input * conj(pn)
        _qLateOut       = _iPNLate*input.imag()         - _qPNLate*input.real();
        filter_out      = Complex(_iEarlyOut, _qEarlyOut);
        filter_out      = _earlyFilter->filterComplexData(filter_out);
        _iEarlyOut      = norm(filter_out);
        filter_out      = Complex(_iLateOut, _qLateOut);
        filter_out      = _lateFilter->filterComplexData(filter_out);
        _iLateOut       = norm(filter_out);                                     // absolute value squared;
        break;
      case REAL_PHASE_DETECTOR:
        _iEarlyOut      = _iPNEarly*input.real();                               // input * pn
        _iLateOut       = _iPNLate*input.real();                                // input * pn
        _qEarlyOut      = _qLateOut     = 0.;
        filter_out      = Complex(_iEarlyOut, _qEarlyOut);
        filter_out      = _earlyFilter->filterComplexData(filter_out);
        _iEarlyOut      = norm(filter_out);
        filter_out      = Complex(_iLateOut, _qLateOut);
        filter_out      = _lateFilter->filterComplexData(filter_out);
        _iLateOut       = norm(filter_out);                                     // absolute value squared;
        break;
  }
//
// phase detector output = early - late
//
  _pdOutputReal         += _iEarlyOut - _iLateOut;
  _integrateSamples++;
  if(_integrateSamples >= _integrateSize)
  {
    _pdOutputReal               /= _integrateSize;
  }
//
// Select output to send:
//
  switch(_phaseDetectorType)
  {
     case COMPLEX_PHASE_DETECTOR:
     case REAL_PHASE_DETECTOR:
     default:
       _pdOutput        = Complex(_pdOutputReal, 0.);
       break;
     case EXTERNAL_PD:
       _pdOutput        = Complex(_iPNEarly, 0.);
       break;
     case 4:
      _pdOutput         = Complex(_iPNOnTime, 0.);
      break;                                                            // input to loop filter
     case 5:
      _pdOutput         = Complex(_qPNEarly, 0.);
      break;                                                            // input to loop filter
     case 6:
      _pdOutput         = Complex(_qPNLate, 0.);
      break;
     case 7:
       _pdOutput        = Complex(_iPNEarly, 0.);
       break;
     case 8:
       _pdOutput        = Complex(_iPNLate, 0.);
       break;
     case 9:
      _pdOutput         = Complex(_iEarlyOut, 0.);
      break;                                                            // input to loop filter
     case 10:
      _pdOutput         = Complex(_iLateOut, 0.);
      break;
     case 11:
       _pdOutput        = Complex(_qEarlyOut, 0.);
       break;
     case 12:
       _pdOutput        = Complex(_qLateOut, 0.);
       break;
  }

  return;
}

// ############################# Private Function ###############################
// findPDVCOGain -- Find the gain of the phase detector times the VCO gain
//
// Input:       inputAmplitude: Maximum amplitude of input signal
//
// Output:                      The phase detector * VCO gain
//
// Notes:
// ############################# Private Function ###############################
float PNLockNonCoh::findPDVCOGain(float inputAmplitude)
{
  float pd_gain, vco_gain, pd_vco_gain;
//
// Assume  VCO gain = fs rad/s/v
// Phase detector gain = A/PI v/rad
// -> pd_vco_gain is in 1/sec
//
  vco_gain              = _samplingFrequency;
  switch(_phaseDetectorType)
  {
      case COMPLEX_PHASE_DETECTOR:
        pd_gain         = 4*inputAmplitude*inputAmplitude/PI;   // (input_i + j input_q)*(pn_i + j pn_q) 
        break;
      case REAL_PHASE_DETECTOR:
      default:
        pd_gain         = inputAmplitude*inputAmplitude/PI;     // (input_i)*(pn_i)
        break;
  }
  pd_vco_gain           = pd_gain*vco_gain;
  return pd_vco_gain;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the PNLockNonCoh class.
//
// Input:               loopFilter:     type of loop filter
//                      output:         type of output
//                      input:          amplitude of input sine wave
//                      zeta:           dampling factor
//                      fn:             loop bandwidth in Hertz
//                      sampling:       sampling frequency in Hertz
//                      tune:           initial VCO offset frequency in Hertz
//                      pnGen:          passed in PN generator (for sharing)
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
PNLockNonCoh::PNLockNonCoh(int loopFilter, int output, float input, float zeta, float fn, float sampling,
                         float tune, int pnPolyI, IS95PNGen *gen)
            :PNLock(loopFilter, output, input, zeta, fn, sampling, tune, pnPolyI, gen)
{
  _earlyFilter          = _lateFilter   = _onTimeFilter = NULL;
  setSamplingAndTune(sampling, tune);
  initializePLL();
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the PNLockNonCoh class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
PNLockNonCoh::~PNLockNonCoh()
{
  delete        _earlyFilter;
  delete        _lateFilter;
  delete        _onTimeFilter;
  return;
}

// ############################# Public Function ###############################
// initializePLL -- Initializes the PLL parameters and loop filter for loop operation.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// 1. We assume a full quadrant multiply of a complex by a complex so that the
//    phase detector gain = Im[A(cos{wt}+jsin{wt})*(cos{wt+phi}-jsin{wt+phi})]
//                        ~ Asin(phi) ~ A*phi.
//    if this is not the case, set _inputAmplitude accordingly.
// ############################# Public Function ###############################
void PNLockNonCoh::initializePLL()
{
//
// Perform super's function
//
  PNLock::initializePLL();
//
// Also, initialize BPFs
//
  if(_earlyFilter == NULL)
  {
    _earlyFilter        = new IIRFilter(BAND_PASS, _bpfFrequency, _bpfFrequency/BW_FACTOR, _samplingFrequency, FILTER_ORDER);
    _lateFilter         = new IIRFilter(BAND_PASS, _bpfFrequency, _bpfFrequency/BW_FACTOR, _samplingFrequency, FILTER_ORDER);
    _onTimeFilter       = new IIRFilter(BAND_PASS, _bpfFrequency, _bpfFrequency/BW_FACTOR, _samplingFrequency, FILTER_ORDER);
  }
  _earlyFilter->zeroOutTaps();
  _lateFilter->zeroOutTaps();
  _onTimeFilter->zeroOutTaps();

  return;
}


// ############################# Public Function ###############################
// setSamplingAndTune -- Sets a new sampling frequency in Hertz.
//
// Input:       sampling:       new sampling frequency in Hertz
//              tune:           new tune frequency in Hertz
//
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void PNLockNonCoh::setSamplingAndTune(float sampling, float tune)
{
//
// Perform super's functions
//
  PNLock::setSamplingFrequency(sampling);
  PNLock::setTuneFrequency(tune);
//
// Also adjust the early and late band pass filter responses:
//
  if(_earlyFilter == NULL)
  {
    _earlyFilter        = new IIRFilter(BAND_PASS, _bpfFrequency, _bpfFrequency/BW_FACTOR, _samplingFrequency, FILTER_ORDER);
    _lateFilter         = new IIRFilter(BAND_PASS, _bpfFrequency, _bpfFrequency/BW_FACTOR, _samplingFrequency, FILTER_ORDER);
    _onTimeFilter       = new IIRFilter(BAND_PASS, _bpfFrequency, _bpfFrequency/BW_FACTOR, _samplingFrequency, FILTER_ORDER);
  }
  else
  {
    _earlyFilter->setFilterParameters(_bpfFrequency, _bpfFrequency/BW_FACTOR, _samplingFrequency);
    _lateFilter->setFilterParameters(_bpfFrequency, _bpfFrequency/BW_FACTOR, _samplingFrequency);
    _onTimeFilter->setFilterParameters(_bpfFrequency, _bpfFrequency/BW_FACTOR, _samplingFrequency);
  }
  _earlyFilter->findTransferFunction();
  _lateFilter->findTransferFunction();
  _onTimeFilter->findTransferFunction();

  return;
}

// ############################# Public Function ###############################
// setBPFFrequency -- Sets a new BPF center frequency for non-coherent loop.
//
// Input:       fo:             new center frequency in Hertz
//
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void PNLockNonCoh::setBPFFrequency(float fo)
{
  _bpfFrequency = MAX(0., fo);
//
// Adjust the early and late band pass filter responses:
//
  if(_earlyFilter == NULL)
  {
    _earlyFilter        = new IIRFilter(BAND_PASS, _bpfFrequency, _bpfFrequency/BW_FACTOR, _samplingFrequency, FILTER_ORDER);
    _lateFilter         = new IIRFilter(BAND_PASS, _bpfFrequency, _bpfFrequency/BW_FACTOR, _samplingFrequency, FILTER_ORDER);
    _onTimeFilter       = new IIRFilter(BAND_PASS, _bpfFrequency, _bpfFrequency/BW_FACTOR, _samplingFrequency, FILTER_ORDER);
  }
  else
  {
    _earlyFilter->setFilterParameters(_bpfFrequency, _bpfFrequency/BW_FACTOR, _samplingFrequency);
    _lateFilter->setFilterParameters(_bpfFrequency, _bpfFrequency/BW_FACTOR, _samplingFrequency);
    _onTimeFilter->setFilterParameters(_bpfFrequency, _bpfFrequency/BW_FACTOR, _samplingFrequency);
  }
  _earlyFilter->findTransferFunction();
  _lateFilter->findTransferFunction();
  _onTimeFilter->findTransferFunction();

  return;
}

