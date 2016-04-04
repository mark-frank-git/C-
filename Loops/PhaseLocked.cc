/************************************************************************
 *                                                                      *
 * This class implements a phase-locked loop.  It designs and imple-    *
 * ments different types of phase-locked loops.                         *
 *                                                                      *
 * File:PhaseLocked.cc                                                  *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 12/19/97  - Started.                                             *
 *  2. 01/27/98  - Went to phaseDetectorType, data directed PLL.        *
 *  3. 11/23/98  - Made this an abstract class, see CarrierLock, Code-  *
 *                 lock subclasses.                                     *
 *  4. 11/30/98  - Move discriminator, VCO gain calculations to sub-    *
 *                 classes.                                             *
 ************************************************************************/

#include "PhaseLocked.h"                                        // Object prototypes
#include <C_Libraries/constants.h>
#include <stdio.h>
#include <math.h>

#ifndef EPS
#define EPS     1.e-5
#endif

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )

// ############################# Private Function ###############################
// loopFilter -- Filters a double input through the loop filter.
//
// Input:       input:          next double input
//
// Output:                      None
//
// Notes:
// ############################# Private Function ###############################
double PhaseLocked::loopFilter(double input)
{
  int           i, i_minus_1;
  double        yout;

//
// First consider case when _numberNumCoeffs == _numberDenCoeffs
// This should always be the case:
//
//
// Check for boundary case:
//
    if( _numberDenCoeffs < 2 )
    {
       return input*_numCoeffs[0];
    }
//
// Else, implement denominator of IIR filter:
//
    yout        = _numCoeffs[1]*_shiftArray[0];
    for(i=(_numberDenCoeffs-1); i>1 ;i--)
    {
       i_minus_1                = i-1;
       yout                     += _numCoeffs[i]*_shiftArray[i_minus_1];
       input                    -= _denCoeffs[i]*_shiftArray[i_minus_1];
       _shiftArray[i_minus_1] = _shiftArray[i-2];
    }
    input                       -= _denCoeffs[1]*_shiftArray[0];
    _shiftArray[0]              = input;
    yout                        += _numCoeffs[0]*input;

  return yout;
}

// ############################# Private Function ###############################
// loopLocked -- Checks to see if the loop has locked
//
// Input:       input:          phase detector output
//
// Output:                      YES if lock detector meets threshold
//
// Notes:
// 1. This function should be overridden by subclasses needing special lock detector
//    behavior.
// ############################# Private Function ###############################
LOGICAL PhaseLocked::loopLocked(Complex &input)
{
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
// 1. This function needs to be overridden in subclasses, it is just a place holder.
// ############################# Private Function ###############################
void PhaseLocked::findPDOutput(Complex &input)
{
  return;
}

// ############################# Private Function ###############################
// findSamplingTime -- Returns the sampling time in seconds
//
// Input:                       None
//
// Output:                      sampling time in seconds
//
// Notes:
// 1. This should be overriden in any class that has a sampling time modifier,
//    for example an integrate and dump.
// ############################# Private Function ###############################
float PhaseLocked::findSamplingTime()
{
  return        1/_samplingFrequency;
}

// ############################# Private Function ###############################
// processPDOutput -- Processes the output of the phase detector
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// 1. This is a dummy function, see the subclasses for real implementations.
// ############################# Private Function ###############################
void PhaseLocked::processPDOutput()
{
  return;
}

// ############################# Private Function ###############################
// findPDVCOGain -- Find the gain of the phase detector times the VCO gain
//
// Input:       inputAmplitude: Maximum input amplitude
//
// Output:                      The phase detector * VCO gain
//
// Notes:
// 1. This function should be overridden in subclasses.  The code below works
//    for a carrier recovery loop.
// ############################# Private Function ###############################
float PhaseLocked::findPDVCOGain(float inputAmplitude)
{
  float pd_vco_gain;
//
// Assume  VCO gain = fs rad/s/V
// where loop_gain = K*PD_GAIN*VCO_GAIN
//
  switch(_phaseDetectorType)
  {
    case COMPLEX_PHASE_DETECTOR:
    default:
      pd_vco_gain       = _samplingFrequency*inputAmplitude;            // exp(-jphi)*exp(-jPhi)
      break;
    case REAL_PHASE_DETECTOR:
      pd_vco_gain       = _samplingFrequency*inputAmplitude/2.;         // cos(phi)*cos(Phi)
      break;
  }
  return pd_vco_gain;
}

#define TAU1_TAU2_RATIO         10
// ############################# Private Function ###############################
// initializeLoopFilterFor -- Initializes the loop filter coefficients and shifts.
//
// Input:       type:           type of loop filter
//              inputAmplitude: maximum amplitude of input signal
//              initializeTaps: YES = initialize taps to 0
//
// Output:                      None
//
// Notes:
// 1. The backwards difference approximation is used: s = (z-1)/zT.
// 2. See Ziemer and Peterson P. 265 for noise bandwidths
// ############################# Private Function ###############################
void PhaseLocked::initializeLoopFilterFor(int type, float inputAmplitude, LOGICAL initializeTaps)
{
  int           i;
  float         sampling_time, loop_gain, filter_gain;
  float         tau, tau1, tau2, omega_n, pd_vco_gain;
//
// Allocate array space:
//
  if(_numCoeffs == NULL)
  {
    _numCoeffs          = new double[MAX_ORDER];
    _denCoeffs          = new double[MAX_ORDER];
    _shiftArray         = new double[MAX_ORDER];
  }
  if(initializeTaps)
    for(i=0; i<MAX_ORDER; i++)
      _shiftArray[i]    = 0.;
//
// Calculate the phase detector * vco gain, the
// called function should be overridden in subclasses.
//
  pd_vco_gain           = findPDVCOGain(inputAmplitude);
//
// Calculate other loop constants:
//
  sampling_time         = findSamplingTime();
  if(_loopLocked)
    omega_n             = TWOPI*_lockedLoopBandwidth;
  else
    omega_n             = TWOPI*_loopBandwidth;
//
// Calculate filter coefficients based on loop filter type:
//
  switch(_loopFilterType)
  {
    default:
    case NO_FILTER:                                     // F(s) = 1, F(z) = 1
      loop_gain         = omega_n;                      // K = wn
      filter_gain       = loop_gain/pd_vco_gain;        // 2 sided noise BW = K/2 in Hz
      _numCoeffs[0]     = filter_gain;
      _numCoeffs[1]     = 0.;                           // Note: _denCoeffs aren't used
      _loopOrder        = 1;
      _loopType         = 1;
      _numberNumCoeffs  = 1;
      _numberDenCoeffs  = 1;
      break;
    case ONE_POLE:                                      // F(s) = 1/(1+sTau)
                                                        // F(z) = zT/[z(T+Tau) - Tau]
      loop_gain         = omega_n/2./_dampingFactor;    // K = wn/(2*zeta)
      filter_gain       = loop_gain/pd_vco_gain;        // 2 sided noise BW = K/2 in Hz
      tau               = 0.5/_dampingFactor/omega_n;   // Tau = 1./(2.*zeta*wn)
      _numCoeffs[0]     = filter_gain*sampling_time/(sampling_time + tau);
      _numCoeffs[1]     = 0.;
      _denCoeffs[0]     = 1.;
      _denCoeffs[1]     = -tau/(sampling_time+tau);
      _loopOrder        = 2;
      _loopType         = 1;
      _numberNumCoeffs  = 2;
      _numberDenCoeffs  = 2;
      break;
    case PHASE_LEAD:                                    // F(s) = (1+sTau)/s
                                                        // F(z) = (zT + [z-1]Tau)/(z-1)
      loop_gain         = omega_n*omega_n;              // K = wn^2
      filter_gain       = loop_gain/pd_vco_gain;        // 2 side Noise BW = (zeta+1/zeta)*wn
      tau               = 2.*_dampingFactor/omega_n;    // Tau = 2*zeta/wn
      _numCoeffs[0]     = (sampling_time + tau)*filter_gain;
      _numCoeffs[1]     = -tau*filter_gain;
      _denCoeffs[0]     = 1.;
      _denCoeffs[1]     = -1.;
      _loopOrder        = 2;
      _loopType         = 2;
      _numberNumCoeffs  = 2;
      _numberDenCoeffs  = 2;
      break;
    case ONE_POLE_PHASE_LEAD:                   // F(s) = (1+sTau2)/(1+sTau1)
                                                // F(z) = (z[{T+Tau2}/{T+Tau1}]-Tau2/{T+Tau1}) /
                                                //        (z - Tau1/{T+Tau1})
                                                // assume Tau1 = 10*Tau2
                                                // 2 sided Noise BW = K*Tau2*(1/Tau2^2+K/Tau1) /
                                                // (2*(K+1/Tau2))
      loop_gain         = omega_n*(TAU1_TAU2_RATIO*_dampingFactor+
                            sqrt(TAU1_TAU2_RATIO*(TAU1_TAU2_RATIO*_dampingFactor*_dampingFactor - 1.)));
      tau1              = loop_gain/omega_n/omega_n;
      tau2              = tau1/TAU1_TAU2_RATIO;
      filter_gain       = loop_gain/pd_vco_gain;
      _numCoeffs[0]     = filter_gain*(sampling_time + tau2)/(sampling_time + tau1);
      _numCoeffs[1]     = -filter_gain*tau2/(sampling_time+tau1);
      _denCoeffs[0]     = 1.;
      _denCoeffs[1]     = -tau1/(sampling_time + tau1);
       _loopOrder       = 2;
      _loopType         = 1;
      _numberNumCoeffs  = 2;
      _numberDenCoeffs  = 2;
      break;
  }
    
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the PhaseLocked class.
//
// Input:               loopFilter:     type of loop filter
//                      output:         type of output
//                      input:          amplitude of input sine wave
//                      zeta:           dampling factor
//                      fn:             loop bandwidth in Hertz
//                      sampling:       sampling frequency in Hertz
//                      tune:           initial VCO offset frequency in Hertz
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
PhaseLocked::PhaseLocked(int loopFilter, int output, float input, float zeta, float fn, float sampling,
                         float tune)
{
//
// Initialize instance variables:
//
  setLoopFilterType(loopFilter);
  setOutputType(output);
  setLockDetectorType(NO_LOCK_DETECTOR);
  setInputAmplitude(input);
  setDampingFactor(zeta);
  setLoopBandwidth(fn);
  setSamplingFrequency(sampling);
  setTuneFrequency(tune);
  setLockedLoopBandwidth(fn/LOCK_BW_SCALE);
  setLockThreshold(DEFAULT_LOCK_THRESH);
  setLockDetectorCount(DEFAULT_LOCK_COUNT);

  _numberNumCoeffs              = _numberDenCoeffs      = 0;
  _oldNumberOutput              = 0;
  _numCoeffs                    = _denCoeffs            = NULL;
  _outputArray                  = NULL;
  _shiftArray                   = NULL;
  _openLoopOperation            = NO;
  _phaseDetectorType            = COMPLEX_PHASE_DETECTOR;

//
// Initialize the loop:
//
  initializePLL();
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the PhaseLocked class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
PhaseLocked::~PhaseLocked()
{
  delete [] _numCoeffs;
  delete [] _denCoeffs;
  delete [] _outputArray;
  delete [] _shiftArray;
  return;
}

// ############################# Public Function ###############################
// setLoopFilterType -- Sets a new loop filter type.
//
// Input:       type:           new loop filter
//
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void PhaseLocked::setLoopFilterType(int type)
{
  _loopFilterType       = type;
  return;
}

// ############################# Public Function ###############################
// setOutputType -- Sets a new loop output.
//
// Input:       type:           new loop output
//
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void PhaseLocked::setOutputType(int type)
{
  _outputType   = MAX(0, type);
  return;
}

// ############################# Public Function ###############################
// setLockDetectorType -- Sets a lock detector type.
//
// Input:       type:           new lock detector type
//
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void PhaseLocked::setLockDetectorType(int type)
{
  _lockDetectorType     = MAX(0, type);
  _lockDetectorType     = MIN(MAX_DET_TYPE, _lockDetectorType);
  return;
}

// ############################# Public Function ###############################
// setInputAmplitude -- Sets a new input sine wave amplitude.
//
// Input:       amplitude:      new sine wave amplitude
//
// Output:                      None
//
// Notes:
// 1. The phase detector gain will be equal to the input amplitude divided by
//    two.
// ############################# Public Function ###############################
void PhaseLocked::setInputAmplitude(float amplitude)
{
  _inputAmplitude       = MAX(EPS, amplitude);
  return;
}

// ############################# Public Function ###############################
// setSamplingFrequency -- Sets a new sampling frequency in Hertz.
//
// Input:       sampling:       new sampling frequency in Hertz
//
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void PhaseLocked::setSamplingFrequency(float sampling)
{
  _samplingFrequency    = MAX(EPS, sampling);
  setTuneFrequency(_tuneFrequency);             // _vcoOffset depends on sampling freq
  return;
}

// ############################# Public Function ###############################
// setDampingFactor -- Sets a new loop damping factor.
//
// Input:       zeta:           new loop damping factor
//
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void PhaseLocked::setDampingFactor(float zeta)
{
  _dampingFactor        = MAX(MIN_ZETA, zeta);
  return;
}

// ############################# Public Function ###############################
// setLoopBandwidth -- Sets a new loop bandwidth in Hertz.
//
// Input:       fn:             new loop bandwidth in Hertz
//
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void PhaseLocked::setLoopBandwidth(float fn)
{
  _loopBandwidth        = MAX(EPS, fn);
  return;
}

// ############################# Public Function ###############################
// setLockedLoopBandwidth -- Sets a new locked loop bandwidth in Hertz.
//
// Input:       fn:             new loop bandwidth in Hertz
//
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void PhaseLocked::setLockedLoopBandwidth(float fn)
{
  _lockedLoopBandwidth  = MAX(EPS, fn);
  return;
}

// ############################# Public Function ###############################
// setLockDetectorCount -- Sets a new value for the lock detector counter.
//
// Input:       offset:         Offset in chips
//
// Output:                      None
//
// Notes:
// 1. If the offset is negative, then delay getting chips from the PN generator
//    until offset is achieved.  Otherwise, pre-fetch chips and throw away.
// ############################# Public Function ###############################
void PhaseLocked::setLockDetectorCount(int count)
{
  _lockDetectorCount    = MAX(1, count);
  return;
}

// ############################# Public Function ###############################
// setLockThreshold -- Sets a new lock threshold value.
//
// Input:       thresh:         New threshold value
//
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void PhaseLocked::setLockThreshold(float thresh)
{
  _lockThreshold        = thresh;
  return;
}

// ############################# Public Function ###############################
// setTuneFrequency -- Sets a new VCO offset frequency in Hertz
//
// Input:       tune:           new VCO offset frequency in Hertz
//
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void PhaseLocked::setTuneFrequency(float tune)
{
  _tuneFrequency        = tune;
  if(_samplingFrequency > 0.)
    _vcoOffset          = TWOPI*_tuneFrequency/_samplingFrequency;
  else
    _vcoOffset          = 0.;
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
void PhaseLocked::initializePLL()
{
  
  _vcoOutput                    = 0.;
  _unModuloVCOOutput            = 0.;
  _vcoSlope                     = 0.;
  _numberSlopePoints            = 0;
  _lockDetectCounter            = 0;
  _peakDetect                   = 0.;
  _loopFilterOutput             = 0.;
  _integratedFilterOutput       = 0.;
  _loopLocked                   = NO;

  initializeLoopFilterFor(_loopFilterType, _inputAmplitude);
  return;
}

// ############################# Private Function ###############################
// pllOutputFor -- Processes a new complex input through the phase locked loop
//
// Input:       input:          next complex input
//
// Output:                      None
//
// Notes:
// 1. This function must be overriden in subclasses, it is just a place holder.
// ############################# Private Function ###############################
Complex PhaseLocked::pllOutputFor(Complex &input)
{
  return _pllOutput;
}

// ############################# Private Function ###############################
// pllOutputFor -- Processes an array of complex input through the phase locked loop
//
// Input:       input:          array of complex inputs
//              number:         size of array
//
// Output:                      None
//
// Notes:
// 1. This function must be overriden in subclasses, it is just a place holder.
// ############################# Private Function ###############################
Complex *PhaseLocked::pllOutputFor(Complex *input, int number)
{
  return _outputArray;
}


// ############################# Private Function ###############################
// pllOutputFor -- Processes an array of real and imaginary array data through
//                 the phase locked loop
//
// Input:       realInput:      array of real inputs
//              imagInput:      array of imag inputs
//              number:         size of array
//
// Output:                      None
//
// Notes:
// 1. This function must be overriden in subclasses, it is just a place holder.
// ############################# Private Function ###############################
Complex *PhaseLocked::pllOutputFor(const float *realInput, const float *imagInput, int number)
{
  return _outputArray;
}
