/************************************************************************
 *                                                                      *
 * This subclass of PhaseLocked supports locking to a carrier.          *
 *                                                                      *
 * File:CarrierLock.cc                                                  *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 11/23/01  - Subclassed from PhaseLocked.                         *
 ************************************************************************/

#include "CarrierLock.h"                                        // Object prototypes
#include <C_Libraries/constants.h>
#include <stdio.h>
#include <math.h>

#ifndef EPS
#define EPS     1.e-5
#endif

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ABS(a)    ((a) >= 0. ? (a) : (-a))

// ############################# Private Function ###############################
// loopLocked -- Checks to see if the loop has locked
//
// Input:       input:          phase detector output
//
// Output:                      YES if lock detector meets threshold
//
// Notes:
// ############################# Private Function ###############################
LOGICAL CarrierLock::loopLocked(Complex &input)
{
  switch(_lockDetectorType)
  {
    case NO_LOCK_DETECTOR:
      return NO;
    case Q_LOCK_DETECTOR:
      _lockDetectorSum  += input.imag();
      if(++_lockDetectCounter >= _lockDetectorCount)
      {
        _lockDetectorSum        /= _lockDetectorCount;
        if(ABS(_lockDetectorSum) <= _lockThreshold)                     // Q should go to zero
          return YES;
        _lockDetectorSum        = 0.;
        _lockDetectCounter      = 0;
      }
      break;
    case I_LOCK_DETECTOR:
      _lockDetectorSum          += input.real();
      if(++_lockDetectCounter >= _lockDetectorCount)
      {
        _lockDetectorSum        /= _lockDetectorCount;
        if(ABS(_lockDetectorSum) >= _lockThreshold)                     // I should go high
          return YES;
        _lockDetectorSum        = 0.;
        _lockDetectCounter      = 0;
      }
      break;
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
// ############################# Private Function ###############################
void CarrierLock::findPDOutput(Complex &input)
{
  Complex       carg;
  double        pd_real;
//
// Phase detector output:
//
  switch(_phaseDetectorType)
  {
    case COMPLEX_PHASE_DETECTOR:
    default:
      carg              = Complex(0., -_vcoOutput);
      _oscillatorOutput = exp(carg);                                    // exp(-j*phi)
      _pdOutput         = _oscillatorOutput*input;
      break;
    case REAL_PHASE_DETECTOR:
      pd_real           = cos(_vcoOutput);
      _oscillatorOutput = Complex(pd_real, 0.);
      pd_real           *= input.real();                                // Real phase detector
      _pdOutput         = Complex(pd_real, pd_real);                    // Store in re and im for
      break;                                                            // input to loop filter
  }
  return;
}


// ############################# Private Function ###############################
// findPDOutput -- Processes a new complex input through the phase detector
//
// Input:       input:          next complex input
//              sliceData       sliced reference data
//
// Output:                      None
//
// Notes:
// 1. The instance variable, _pdOutput, is modified.
// 2. This function implements the PD for a data directed PLL.
// 3. For the data directed PLL, we assume the complex down conversion is done
//    prior to this block, with the sliced data taken after the down conversion
// ############################# Private Function ###############################
void CarrierLock::findPDOutput(Complex &input, Complex &sliceData)
{
  Complex       carg;
  double        slice_magnitude;

  carg                  = Complex(0., -_vcoOutput);
  _oscillatorOutput     = exp(carg);                            // exp(-j*phi)
  _zn                   = input;                                // down conversion performed 'off chip'
  slice_magnitude       = norm(sliceData);
//
// Phase detector output:
//
  if(slice_magnitude != 0.)
  {
    carg        = conj(sliceData);                              // see eqn 6.70 in Gitlin
    _pdOutput   = _zn*carg/slice_magnitude;
  }
  else
    _pdOutput   = _zn;
  return;
}

// ############################# Private Function ###############################
// processPDOutput -- Processes the output of the phase detector
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// 1. The instance variables, _loopFilterOutput, _vcoOutput, _pdOutput, _oscillatorOutput
//    contain the outputs from the various parts of the loop.
// 2. The PLL operates as follows: the input complex signal is multiplied by the
//    complex sine wave from the VCO, the imaginary component (Q) of the product
//    is then sent through the loop filter, and back around to the VCO.  For a
//    NTSC/PAL signal, the demodulated output will be the I component of the
//    VCO/input signal product.
// 3. You must first call findPDOutput().
// ############################# Private Function ###############################
void CarrierLock::processPDOutput()
{
//
// First, check for loop locked (this call should probably be located elsewhere):
//
  if(!_loopLocked)
  {
    _loopLocked                 = loopLocked(_pdOutput);
    if(_loopLocked)
    {
      _loopLockCounter          = 0;
      _integratedFilterOutput   = 0.;
      _loopRelocked             = NO;
    }
  }
  else if( (_loopLockCounter++ > LOCK_DETECT_COUNT) && !_loopRelocked)
  {
    _vcoOffset          += _integratedFilterOutput/LOCK_DETECT_COUNT;   // Adjust tuning frequency
    initializeLoopFilterFor(_loopFilterType, _inputAmplitude, NO);      // Set new loop bandwidth, don't init
                                                                        // taps
    _loopRelocked       = YES;
  }
//
// Loop filter output:
//
  if(!_openLoopOperation)
  {
    _loopFilterOutput   = loopFilter(_pdOutput.imag());
//
// VCO output, mod it 2PI to prevent overflows.
// NOTE: _vcoOffset is used to pre-tune the VCO
// NOTE: idle mode option is no longer functional
//
    _vcoOutput                  += _vcoOffset + _loopFilterOutput;
    _unModuloVCOOutput          += _vcoOffset + _loopFilterOutput;      // Not mod'd by TWOPI
    _integratedFilterOutput     += _loopFilterOutput;
  }
  else if(_numberSlopePoints > 0)
    _vcoOutput          += _vcoSlope/_numberSlopePoints;        // Use old slope for idle mode
  if(_vcoOutput > TWOPI)
    _vcoOutput          -= TWOPI;
  else if(_vcoOutput < -TWOPI)
    _vcoOutput          += TWOPI;
  _oldVCOOutput         = _vcoOutput;
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the CarrierLock class.
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
CarrierLock::CarrierLock(int loopFilter, int output, float input, float zeta, float fn, float sampling,
                         float tune)
            :PhaseLocked(loopFilter, output, input, zeta, fn, sampling, tune)
{
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the CarrierLock class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
CarrierLock::~CarrierLock()
{
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
// ############################# Private Function ###############################
Complex CarrierLock::pllOutputFor(Complex &input)
{

  findPDOutput(input);
  processPDOutput();                    // Processes input, saves output in instance variables
  switch (_outputType)
  {
    default:
    case PD_OUTPUT:
    case DEMOD_OUTPUT:
      return _pdOutput;
    case LOOP_FILTER_OUTPUT:
      _pllOutput        = Complex(_loopFilterOutput, 0.);
      break;
    case INTEGRATED_LOOP_FILTER:
      _pllOutput        = Complex(_integratedFilterOutput, _phaseRotationOutput);
      break;
    case VCO_OUTPUT:
      _pllOutput        = Complex(_unModuloVCOOutput, 0.);
      break;
    case OSCILLATOR_OUTPUT:
      return _oscillatorOutput;
    case VCO_SLOPE:
      if(_numberSlopePoints > 0)
        _pllOutput      = Complex(_vcoSlope/_numberSlopePoints, 0.);
      else
        _pllOutput      = Complex(0., 0.);
      break;
    case LOCK_DETECTOR:
      _pllOutput        = Complex((double)_lockDetectorSum, (double)_lockDetectorSum);
      break;
    case LOCK_STATUS:
      _pllOutput        = Complex((double)_loopLocked, 0.);
      break;
  }

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
// ############################# Private Function ###############################
Complex *CarrierLock::pllOutputFor(Complex *input, int number)
{
  int           i;
  Complex       output;
//
// Allocate output:
//
  if(number > _oldNumberOutput)
  {
    _oldNumberOutput    = number;
    delete [] _outputArray;
    _outputArray        = new Complex[number];
  }
//
// Loop over number of input samples:
//
  for(i=0; i<number; i++)
  {
    _outputArray[i]     = pllOutputFor(input[i]);
  }
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
// ############################# Private Function ###############################
Complex *CarrierLock::pllOutputFor(const float *realInput, const float *imagInput, int number)
{
  int           i;
  Complex       output, input;
//
// Allocate output:
//
  if(number > _oldNumberOutput)
  {
    _oldNumberOutput    = number;
    delete [] _outputArray;
    _outputArray        = new Complex[number];
  }
//
// Loop over number of input samples:
//
  for(i=0; i<number; i++)
  {
    input               = Complex(realInput[i], imagInput[i]);
    _outputArray[i]     = pllOutputFor(input);
  }
  return _outputArray;
}

// ############################# Public Function ###############################
// pllOutputForSlice -- Processes a new complex input through the phase locked loop
//
// Input:       input:          next complex input
//              sliceData:      complex sliced data
//
// Output:                      None
//
// Notes:
// 1. This function implements a data directed PLL.
// 2. For the data directed PLL, we assume the complex down conversion is done
//    prior to this block, with the sliced data taken after the down conversion
// ############################# Public Function ###############################
Complex CarrierLock::pllOutputForSlice(Complex &input, Complex &sliceData)
{

  findPDOutput(input, sliceData);
  processPDOutput();                    // Processes input, saves output in instance variables
  switch (_outputType)
  {
    default:
    case PD_OUTPUT:
      return    _pdOutput;
    case LOOP_FILTER_OUTPUT:
      _pllOutput        = Complex(_loopFilterOutput, 0.);
      break;
    case VCO_OUTPUT:
      _pllOutput        = Complex(_unModuloVCOOutput, 0.);
      break;
    case OSCILLATOR_OUTPUT:
      return    _oscillatorOutput;
    case DEMOD_OUTPUT:
      return    _zn;
    case VCO_SLOPE:
      if(_numberSlopePoints > 0)
        _pllOutput      = Complex(_vcoSlope/_numberSlopePoints, 0.);
      else
        _pllOutput      = Complex(0., 0.);
      break;
  }

  return _pllOutput;
}
