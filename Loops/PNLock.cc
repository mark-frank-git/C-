/************************************************************************
 *                                                                      *
 * This subclass of PhaseLocked supports locking to a PN code(s).       *
 *                                                                      *
 * File:PNLock.cc                                                       *
 *                                                                      *
 ************************************************************************/

#include "PNLock.h"                                     // Object prototypes
#include <C_Libraries/constants.h>
#include <Generators/PNFiltered.h>
#include <Generators/IS95PNGen.h>
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
// 1. We use a different lock detector sum for each sample per chip in case only
//    one of the chips has zero inter-chip interference (ICI).
// ############################# Private Function ###############################
LOGICAL PNLock::loopLocked(Complex &input)
{
  float         despread_i, despread_q;
  float         norm_threshold;

//
// Check for pilot signal:
//
  despread_i                            = input.real()*_iPNOnTime;
  if(_rotationInsensitive)                                      // Constellation may be rotated
  {
    despread_i                          = input.real()*_iPNOnTime;
    despread_q                          = input.imag()*_qPNOnTime;
    _lockDetectorSum                    += despread_i*despread_i + despread_q*despread_q;
  }
  else
    _lockDetectorSum                    += despread_i;

//
// Lock statistic:
//
  if(++_lockDetectCounter >= _lockDetectorCount)                // Should we use _numberChips? and walshChips?
  {
//
// calculate normalized threshold based on input amplitude:
//
    if(_rotationInsensitive)
      norm_threshold            = 2.*_lockThreshold*_inputAmplitude*_inputAmplitude;
    else
      norm_threshold            = _lockThreshold*_inputAmplitude;
    _lockDetectorSum            /= _lockDetectorCount;
    if(_lockDetectorSum >= norm_threshold)
        return YES;
    _lockDetectorSum            = 0.;
    _lockDetectCounter          = 0;
  }
  return NO;
}

#define TAN_BIG_NUMBER  1000.
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
void PNLock::findPDOutput(Complex &input)
{
  int           quad;
  Complex       filter_out;
//
// First, check for loop locked (this call should probably be located elsewhere):
//
  if(!_loopLocked  && (_lockDetectorType != NO_LOCK_DETECTOR) )
  {
    _loopLocked         = loopLocked(input);
    if(_loopLocked)
    {
      initializeLoopFilterFor(_loopFilterType, _inputAmplitude, YES);   // Set new loop bandwidth
    }
  }
      
//
// Get output of early and late arms after quadratic detector:
// For rotation insensitive, we sum the outputs, not really needed for rotation sensitive
//
  switch(_phaseDetectorType)
  {
    case COMPLEX_PHASE_DETECTOR:
    default:
        _iEarlyOut      += _iPNEarly*input.real()       + _qPNEarly*input.imag();       // input * conj(pn)
        _qEarlyOut      += _iPNEarly*input.imag()       - _qPNEarly*input.real();       // needed for rot.
        _iLateOut       += _iPNLate*input.real()        + _qPNLate*input.imag();        // input * conj(pn)
        _qLateOut       += _iPNLate*input.imag()        - _qPNLate*input.real();        // needed for rot.
        _iOnTimeOut     += _iPNOnTime*input.real()      + _qPNOnTime*input.imag();      // input * conj(pn)
        _qOnTimeOut     += _iPNOnTime*input.imag()      - _qPNOnTime*input.real();      // needed for rot.
        break;
    case REAL_PHASE_DETECTOR:
        _iEarlyOut      += _iPNEarly*input.real();                              // input * pn
        _iLateOut       += _iPNLate*input.real();                               // input * pn
        _qEarlyOut      = _qLateOut     = 0.;
        break;
  }
//
// Check for integrate and dump time:
//
  _integrateSamples++;
  if(_integrateSamples >= _integrateSize)
  {
    if(_loopLocked)                                                     // Accumulate estimators if locked
    {
        _timingAverage          += _integratedFilterOutput;
        _timingVariance         += _integratedFilterOutput*_integratedFilterOutput;
        if(_qOnTimeOut > 0.)                                            // count samples for each quadrants
        {
          if(_iOnTimeOut > 0.)
            quad        = 0;
          else
            quad        = 1;
        }
        else
        {
          if(_iOnTimeOut > 0.)
            quad        = 3;
          else
            quad        = 2;
        }
        _quadrantCount[quad]++;
        if(_iOnTimeOut != 0.)
          _qOverIAverage[quad]  += _qOnTimeOut/_iOnTimeOut;
        else
          _qOverIAverage[quad]  += _qOnTimeOut*TAN_BIG_NUMBER;          // +/- 90 degrees
        _numberTimingAvgs++;
    }
    if(_rotationInsensitive)
    {
      _pdOutputReal             = _iEarlyOut*_iEarlyOut + _qEarlyOut*_qEarlyOut;
      _pdOutputReal             -= _iLateOut*_iLateOut + _qLateOut*_qLateOut;
      _pdOutputReal             /= _integrateSize;                      // need to divide by N^2
      _phaseRotationOutput      = atan2(_qOnTimeOut, _iOnTimeOut);      // This can be deleted for fixed pt, use _qRotation
    }
    else
    {
      _pdOutputReal             = _iEarlyOut - _iLateOut;
    }
    _pdOutputReal               /= _integrateSize;
    _iEarlyOut                  = _iLateOut     = _iOnTimeOut   = 0.;   // Reset I/D accumulators
    _qEarlyOut                  = _qLateOut     = _qOnTimeOut   = 0.;
  }
//
// Save the output for display:
//
  switch(_phaseDetectorType)
  {
     case COMPLEX_PHASE_DETECTOR:
     case REAL_PHASE_DETECTOR:
     default:
       _pdOutput        = Complex(_pdOutputReal, _phaseRotationOutput);
       break;
     case 2:                                                            // These are for debug
       _pdOutput        = Complex(_iPNEarly, _qPNEarly);
       break;
     case 3:
       _pdOutput        = Complex(_iPNLate, _qPNLate);
       break;
     case 4:
       _pdOutput        = Complex(_iOnTimeOut, _qOnTimeOut);
       break;                                                           // input to loop filter
     case 5:
       _pdOutput        = Complex(_iEarlyOut, _qEarlyOut);
       break;                                                           // input to loop filter
     case 6:
       _pdOutput        = Complex(_iLateOut, _qLateOut);
       break;
  }

  return;
}

#define INSENSITIVE_GAIN        4.
#define SENSITIVE_GAIN          2.
// ############################# Private Function ###############################
// findPDVCOGain -- Find the gain of the phase detector times the VCO gain
//
// Input:       inputAmplitude: Maximum amplitude of input signal
//
// Output:                      The phase detector * VCO gain
//
// Notes:
// ############################# Private Function ###############################
float PNLock::findPDVCOGain(float inputAmplitude)
{
  float pd_gain, vco_gain, pd_vco_gain;
//
// Assume  VCO gain = fs rad/s/v
// Phase detector gain = A/PI v/rad
// -> pd_vco_gain is in 1/sec
//
  vco_gain              = _samplingFrequency;
  if(_rotationInsensitive)
  {
    switch(_phaseDetectorType)
    {
      case COMPLEX_PHASE_DETECTOR:                               // (input_i + j input_q)*(pn_i + j pn_q)
        pd_gain         = INSENSITIVE_GAIN*inputAmplitude*inputAmplitude/PI; 
        break;
      case REAL_PHASE_DETECTOR:
      default:
        pd_gain         = inputAmplitude*inputAmplitude/PI;     // (input_i)*(pn_i)
        break;
    }
  }
  else                                                          // Not rotation insensitive
  {
    switch(_phaseDetectorType)
    {
      case COMPLEX_PHASE_DETECTOR:
        pd_gain         = SENSITIVE_GAIN*inputAmplitude/PI;     // (input_i + j input_q)*(pn_i + j pn_q)
        break;
      case REAL_PHASE_DETECTOR:
      default:
        pd_gain         = inputAmplitude/PI;                    // (input_i)*(pn_i)
        break;
    }
  }
  pd_vco_gain           = pd_gain*vco_gain;
  return pd_vco_gain;
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
float PNLock::findSamplingTime()
{
  return        _integrateSize/_samplingFrequency;
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
// 4. The variable _qRotation and _iRotation accumulate the outputs of the on
//    time variables in order to find the phase rotation,  in this way, only
//    one call to atan2() is required.
// ############################# Private Function ###############################
void PNLock::processPDOutput()
{
//
// Loop filter output:
//
  if(!_openLoopOperation)
  {
    if(_integrateSamples >= _integrateSize)
    {
      _loopFilterOutput         = loopFilter(_pdOutputReal);
      _integrateSamples         = 0;
    }
//
// VCO output, mod it 2PI to prevent overflows.
// NOTE: _vcoOffset is used to pre-tune the VCO
//
    _vcoOutput                  += _vcoOffset + _loopFilterOutput;
    _unModuloVCOOutput          += _vcoOffset + _loopFilterOutput;      // Not mod'd by TWOPI
    _integratedFilterOutput     += _loopFilterOutput;
  }
  else
  {
    if(_integrateSamples >= _integrateSize)
    {
      _integrateSamples = 0;
    }
    _vcoOutput          += _vcoOffset;
    _unModuloVCOOutput  += _vcoOffset;                          // Not mod'd by TWOPI
  }

//
// The following assumes delta = 1 chip
//
  if(_vcoOutput > TWOPI)
  {
    _numberChips++;
    _vcoOutput          -= TWOPI;
  }
  _oldVCOOutput         = _vcoOutput;
//
// Get the PN generator outputs.
//
  _oscillatorOutput     = Complex(_iPNOnTime, _qPNOnTime);
  _iPNEarly             = _iPN->getEarlyLateOutput(_numberChips+1, _vcoOutput);
  _iPNLate              = _iPN->getEarlyLateOutput(_numberChips, _vcoOutput);
  _qPNEarly             = _qPN->getEarlyLateOutput(_numberChips+1, _vcoOutput);
  _qPNLate              = _qPN->getEarlyLateOutput(_numberChips, _vcoOutput);
  _iPNOnTime            = _iPN->getOnTimeOutput(_numberChips+1);
  _qPNOnTime            = _qPN->getOnTimeOutput(_numberChips+1);
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the PNLock class.
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
PNLock::PNLock(int loopFilter, int output, float input, float zeta, float fn, float sampling,
                         float tune, int pnPolyI, IS95PNGen *gen)
            :PhaseLocked(loopFilter, output, input, zeta, fn, sampling, tune)
{
  if(gen == NULL)
  {
    _iPN                = new PNFiltered(pnPolyI);
    _qPN                = new PNFiltered(IS95_Q_PN);
  }
  else
  {
    _iPN                = new PNFiltered(gen->iData(0), gen->pnLength());
    _qPN                = new PNFiltered(gen->qData(0), gen->pnLength());
  }
  setRotationInsensitive(NO);
  setSamplingAndTune(sampling, tune);
  setPNOffset(DEFAULT_OFFSET);
  setIntegrateSize(DEFAULT_INTEGRATE);
  setLockDetectorCount(ROUND(sampling/tune)*WALSH_CHIPS*2);
  initializePLL();
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the PNLock class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
PNLock::~PNLock()
{
  delete        _iPN;
  delete        _qPN;
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
void PNLock::initializePLL()
{
  int   i;
//
// Perform super's function
//
  PhaseLocked::initializePLL();
//
// Also reset the PN generator:
//
  _iPN->initGenerator();
  _qPN->initGenerator(); 
//
// reset integrators:
//
  _integrateSamples             = 0;
//
// reset estimator averagers:
//
  _timingAverage                = 0.;
  _timingVariance               = 0.;
  _numberTimingAvgs             = 0;
  for(i=0; i<NUMBER_QUADS; i++)
  {
    _quadrantCount[i]           = 0;
    _qOverIAverage[i]           = 0.;
  }
//
// reset accumulators (for rotation insensitive)
//
  _iEarlyOut            = _iLateOut     = _iOnTimeOut   = 0.;
  _qEarlyOut            = _qLateOut     = _qOnTimeOut   = 0.;
//
// Initialize VCO accumulator to get 1/2 sample time advance to synchronize
// with input waveform (a bit of a kluge), also add in _pnOffset to get
// user defined advance/retard:
//
   _vcoOutput                   = _vcoOffset/2.;
   _numberChips                 = _pnOffset;
//
// Reset lock detector variables:
//
   _samplesPerChip      = ROUND(_samplingFrequency/_tuneFrequency);
   _lockDetectCounter   = _pnOffset*_samplesPerChip;
   _lockDetectorSum     = 0.;

//
// Initialize PN outputs:
//
  _iPNEarly             = _iPNLate      = _iPNOnTime    = 0;
  _qPNEarly             = _qPNLate      = _qPNOnTime    = 0;

  return;
}

// ############################# Public Function ###############################
// setIPoly -- Sets a a new I channel PN polynomial in octal.
//
// Input:       poly:           new PN polynomial (octal representation)
//
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void PNLock::setIPoly(int poly)
{
  _iPN->setPNPolynomial(poly);
  return;
}

// ############################# Public Function ###############################
// setQPoly -- Sets a a new Q channel PN polynomial in octal.
//
// Input:       poly:           new PN polynomial (octal representation)
//
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void PNLock::setQPoly(int poly)
{
  _qPN->setPNPolynomial(poly);
  return;
}

// ############################# Public Function ###############################
// setIntegrateSize -- Sets a a new integrate and dump size.
//
// Input:       size:           new integrate size
//
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void PNLock::setIntegrateSize(int size)
{
  _integrateSize        = MAX(1, size);
  initializePLL();
  return;
}

// ############################# Public Function ###############################
// setRotationInsensitive -- Turns on or off rotation insensitivity.
//
// Input:       flag:           YES = loop insensitive to rotations in input constellation
//
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void PNLock::setRotationInsensitive(LOGICAL flag)
{
  _rotationInsensitive          = flag;
  return;
}

// ############################# Public Function ###############################
// setPNOffset -- Sets a new value for the PN offset.
//
// Input:       offset:         Offset in chips
//
// Output:                      None
//
// Notes:
// 1. If the offset is negative, then delay getting chips from the PN generator
//    until offset is achieved.  Otherwise, pre-fetch chips and throw away.
// ############################# Public Function ###############################
void PNLock::setPNOffset(int offset)
{
  _pnOffset             = MIN(MAX_OFFSET, offset);
  _pnOffset             = MAX((-MAX_OFFSET), _pnOffset);
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
void PNLock::setSamplingAndTune(float sampling, float tune)
{
//
// Perform super's functions
//
  PhaseLocked::setSamplingFrequency(sampling);
  PhaseLocked::setTuneFrequency(tune);

  return;
}

// ############################# Public Function ###############################
// timingEstimate -- Returns the timing estimate from the integrated loop filter.
//
// Input:                       None
//
// Output:                      Timing estimate in fractions of a chip
//
// Notes:
// 1. If the input signal is delayed relative to the initial PN offset of the DLL,
//    the timing estimator will be negative, otherwise, it will be positive.
// 2. The integrated loop filter output is scaled for 2*PI radians per chip (of
//    course this is arbitrary).  We need to account for this scaling below.
// ############################# Public Function ###############################
float PNLock::timingEstimate()
{
  float estimate        = 0.;
  if(_numberTimingAvgs > 0)
    estimate    = _timingAverage/TWOPI/_numberTimingAvgs;

  return estimate;
}

// ############################# Public Function ###############################
// timingVariance -- Returns the timing estimate variance from the integrated loop filter.
//
// Input:                       None
//
// Output:                      Timing estimate variance in chips^2
//
// Notes:
// 1. The integrated loop filter output is scaled for 2*PI radians per chip (of
//    course this is arbitrary).  We need to account for this scaling below.
// ############################# Public Function ###############################
float PNLock::timingVariance()
{
  float variance        = 0.;
  if(_numberTimingAvgs > 1)
    variance            = (_timingVariance-_timingAverage*_timingAverage/_numberTimingAvgs)/(_numberTimingAvgs-1);
  variance              /= TWOPI*TWOPI;
  return variance;
}

// ############################# Public Function ###############################
// rotationEstimate -- Returns the phase rotation estimate from the early and
//                     late integrators, given in radians.
//
// Input:                       None
//
// Output:                      phase rotation in radians
//
// Notes:
// ############################# Public Function ###############################
float PNLock::rotationEstimate()
{
  int   i, max_quadrant, max_count;
  float estimate        = 0.;

  max_quadrant  = 0;
  max_count     = _quadrantCount[0];
  for(i=1; i<NUMBER_QUADS; i++)
  {
    if(_quadrantCount[i] > max_count)
    {
      max_quadrant      = i;
      max_count         = _quadrantCount[i];
    }
  }
  if(max_count > 0)
    estimate    = atan(_qOverIAverage[max_quadrant]/max_count);         // result is in -PI/2->PI/2
  if(max_quadrant == 1)
    estimate    += PI;
  if(max_quadrant == 2)
    estimate    -= PI;

  return estimate;
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
Complex PNLock::pllOutputFor(Complex &input)
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
    case PN_EARLY_OUTPUT:
      _pllOutput        = Complex(_iPNEarly, _qPNEarly);
      break;
    case PN_LATE_OUTPUT:
      _pllOutput        = Complex(_iPNLate, _qPNLate);
      break;
    case PN_ON_OUTPUT:
      _pllOutput        = Complex(_iPNOnTime, _qPNOnTime);
      break;
    case DESPREAD_OUTPUT:
      _pllOutput        = Complex(_iPNOnTime, _qPNOnTime);              // For now
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
Complex *PNLock::pllOutputFor(Complex *input, int number)
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
Complex *PNLock::pllOutputFor(const float *realInput, const float *imagInput, int number)
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
