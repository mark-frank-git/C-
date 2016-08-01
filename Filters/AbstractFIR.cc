/************************************************************************
 *                                                                      *
 * This subclass of DigitalFilter adds functionality for filtering      *
 * using an FIR structure.  For calculating filter coefficients, see    *
 * the subclasses.                                                      *
 *                                                                      *
 * File:AbstractFIR.cc                                                  *
 *                                                                      *
 * The filter is stored in the forms:                                   *
 *                                                                      *
 *    H(z) =  b[n] + b[n-1]z + ... + b[0]z^n                            *
 *         =  b[0] + b[1]z^(-1) + ... + b[n]z^-(n)                      *
 *                                                                      *
 * Where the a and b polynomials are stored in the polynomial objects.  *
 * NOTE: The order of the filter is the # taps - 1.                     *
 *                                                                      *
 *                                                                      *
 ************************************************************************/

#include "AbstractFIR.h"                                        // Object prototypes
#if defined(WIN32)
#include <GNU/Complex.h>
#include <Polynomials/DoublePoly.h>
#include <C_Libraries/constants.h>
#else
#include "Complex.h"
#include "DoublePoly.h"
#include "constants.h"
#endif

#include <math.h>
#include <stdio.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))


// ############################# Private Method ###############################
// transferResponseAt -- Finds the response of the filter at a complex point using
//                       the transfer function representation of the filter.
//
// Input:       point:          Complex point in the s plane
//          
// Output:                      None
//
// Notes:
// 1. Overridden from AbstractFilter
// ############################# Private Method ###############################
Complex AbstractFIR::transferResponseAt(Complex &point)
{
  Complex num, complex_response;

//
// Evaluate numerator polynomial at the point:
//
  num                   = _bPolyObject->evaluateAtComplexPoint(point);
  complex_response      = num;
  
  return complex_response;
}

// ############################# Private Method ###############################
// outputResponseFor -- Convert Complex response to form desired by user.  The
//                      response type is determined by the instance variable,
//                      __responseType
//
// Input:       cResponse:          Filter's complex response
//              omega:              Frequency in rad/s of filter's response
//          
// Output:                          A floating point value giving the filter's
//                                  response.
//
// Notes:
// 1. This method overrides the super's method to calculate normalized phase.
// ############################# Private Method ###############################
float AbstractFIR::outputResponseFor(Complex cResponse , float omega)
{
  float output_response;

  switch(_responseType)
  {
    default:
//
// Call super's method
//
      output_response   = AbstractFilter::outputResponseFor(cResponse, omega);
      break;
    case PHASE_NORMALIZED:                      // This is overridden from super
      output_response   = atan2(imag(cResponse), real(cResponse));
      if(ABS((_thetaOld-output_response)) > PI)         // Check for 0 crossing
      {
        if(_thetaOld > 0.)
          _subAngle     -= TWOPI;
        else
          _subAngle     += TWOPI;
      }
      _thetaOld         = output_response;
      output_response   -= _subAngle;
      output_response   -= omega*(_numberTaps-1.)/2./_fs;
      if(_phaseInDegrees)
        output_response *= DEG_RAD;
      break;
  }
  return output_response;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the AbstractFIR class.
//
// Input:           pass:           Filter pass type: LOW_PASS, BAND_PASS, etc.
//                  centerFreq:     Filter center frequency in Hertz
//                  cutoffFreq:     Filter cutoff frequency in Hertz
//                  samplingFreq:   Digital sampling frequency in Hertz
//                  order           Filter order
//
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
AbstractFIR::AbstractFIR(int pass, double centerFreq, double cutoffFreq,
                                  double samplingFreq, int order)
              :DigitalFilter(pass, centerFreq, cutoffFreq, samplingFreq, order)
{
  setNumberTaps(order+1);
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the AbstractFIR class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
AbstractFIR::~AbstractFIR()
{
  return;
}

// ############################# Public Method ###############################
// setNumberTaps -- Sets a new number of taps for the filter.
// Input:       taps            New number of filter taps
//          
// Output:                      none
// Notes:
// ############################# Public Method ###############################
void AbstractFIR::setNumberTaps(int taps)
{
  _numberTaps   = MAX(taps, MIN_TAPS);
  _numberTaps   = MIN(_numberTaps, MAX_TAPS);
//
// Call super's method
//
  DigitalFilter::setFilterOrder(_numberTaps-1);
  return;
}

// ############################# Public Method ###############################
// setFilterOrder -- Sets the order of the filter (numberTaps - 1).
// Input:       taps            New number of filter taps
//          
// Output:                      none
// Notes:
// ############################# Public Method ###############################
void AbstractFIR::setFilterOrder(int order)
{
  setNumberTaps(order+1);
  return;
}

// ############################# Public Method ###############################
// filterPoles -- Returns the complex poles of the digital filter.
// Input:                       None
//          
// Output:                      An array of complex pole positions
// Notes:
//  1. To find the size of the array, call numberOfPoles().
//  2. Since we are an FIR filter, we really don't have any poles, so set poles
//     at origin.
// ############################# Public Method ###############################
const Complex *AbstractFIR::filterPoles()
{
  int       i;

  if(_bPolyObject != NULL)
  {
    _numberPoles        = _bPolyObject->getOrder();
    _filterPoles        = new Complex[_numberPoles];
    for(i=0; i<_numberPoles; i++)
    _filterPoles[i]     = Complex(0., 0.);
  }
  
  return _filterPoles;
}

// ############################# Public Method ###############################
// filterZeros -- Returns the complex zeros of the digital filter.
// Input:                       None
//          
// Output:                      An array of complex zeros positions
// Notes:
//  1. To find the size of the array, call numberOfZeros().
//  2. The root finder does not always work, especiall for high order polys
// ############################# Public Method ###############################
const Complex *AbstractFIR::filterZeros()
{
  int       i;
  Complex   *roots;

  if(_bPolyObject != NULL)
  {
    delete [] _filterZeros;
    _numberZeros        = _bPolyObject->getOrder();
    _filterZeros        = new Complex[_numberZeros];
    roots               = _bPolyObject->roots();
    for(i=0; i<_numberZeros; i++)
      _filterZeros[i]   = roots[i];
  }
  return _filterZeros;
}

// ############################# Public Method ###############################
// setPassBandGain -- Calculate and set a new _passBandGain.  Useful when loading
//                    in new filter coefficients
// Input:       usePolesZeros:  If yes, set pass band gain using poles and zeros
//                              else, use transfer function
//          
// Output:                      None
//
// ############################# Public Method ###############################
void AbstractFIR::setPassBandGain(LOGICAL usePolesZeros)
{
  Complex   arg, num, response;
  double    wo, wc;
//
// Calculate the pass band gain:
//
  switch(_filterPassType)
  {
    case LOW_PASS:
    case BAND_STOP:
    default:
      arg               = Complex(1., 0.);                              // omega = 0
      break;
    case HIGH_PASS:
      wc                = TWOPI*_fc/_fs;
      arg               = Complex(0., (wc + PI)/2.);                    // omega = (wc+PI)/2.
      arg               = exp(arg);
      break;
    case BAND_PASS:
      wo                = TWOPI*_fo/_fs;
      arg               = Complex(0., wo);                              // omega = wo;
      arg               = exp(arg);
      break;
  }
//
// Gain = numerator/denominator
//
  if(usePolesZeros)
  {
    response            = poleZeroResponseAt(arg, _filterPoles, _filterZeros);
    _passBandGain       = abs(response);
    _passBandGain       = MAX(1.e-40, _passBandGain);
  }
  else
  {
    num                 = _bPolyObject->evaluateAtComplexPoint(arg);
    _passBandGain       = abs(num);
  }
  return;
}

// ############################# Public Method ###############################
// zeroOutTaps -- checks the size of the shift register, and zeros out the
//                digital filter's taps.
// Input:                       None
//          
// Output:                      None
// ############################# Public Method ###############################
void AbstractFIR::zeroOutTaps()
{
//
// Perform super's function
//
  DigitalFilter::zeroOutTaps();
  _shiftPtr     = MAX((_shiftSize - 1), 0);     // The pointer decrements, so set to end
  return;
}

// ############################# Public Method ###############################
// filterFloatArray -- Filter the input array and overwrite it with output.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   x:              An array of floating point data, modified on output
//          numberPts:      Number of points in the input array.
//          
// Output:                  YES if went OK, No if _bPolyObject doesn't exist.
// ############################# Public Method ###############################
LOGICAL AbstractFIR::filterFloatArray(float * input, int numberPts)
{
  int       i, j, order, shift_ptr;
  double    gain, yout;
  const     double *b, *b_ptr;

   if(_bPolyObject == NULL)
     return NO;
   
   if(_passBandGain != 0.)
     gain   = 1./_passBandGain;
   else
     gain   = 1.;
//
// Now, get the coefficients
//
  b_ptr         = _bPolyObject->getCoefficients();
  order         = _bPolyObject->getOrder();
//
// Filter the data:
//
  for(j=0; j<numberPts; j++)
  {
    b           = b_ptr;
    yout        = (*b) * input[j];
    b++;
    shift_ptr   = _shiftPtr;
    for(i=0; i<order; i++)
    {
      yout      += (*b) * _doubleShift[shift_ptr];
      b++;
      shift_ptr++;
      if(shift_ptr == _shiftSize)
        shift_ptr       = 0;
    }
//
// Add the new data to the shift register
//
    _shiftPtr--;
    if(_shiftPtr < 0)
      _shiftPtr                 = MAX((_shiftSize - 1), 0);
    _doubleShift[_shiftPtr]     = input[j];
//
// Save the output:
//
    input[j]                    = gain*yout;
  }

  return YES;
}

// ############################# Public Method ###############################
// filterFloatData -- Filter an input data sample, and return the sampled value.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   input:          An input data sample
//          
// Output:                  The filtered data sample
// ############################# Public Method ###############################
double AbstractFIR::filterFloatData(float input)
{
  int           i, shift_ptr;
  double        yout, gain;
  const         double *b;

   if(_bPolyObject == NULL)
     return 0.;
   
   if(_passBandGain != 0.)
     gain       = 1./_passBandGain;
   else
     gain       = 1.;
//
// Now, get the coefficients
//
  b         = _bPolyObject->getCoefficients();
//
// Filter the data:
//
  yout          = (*b) * input;
  b++;
  shift_ptr     = _shiftPtr;
  for(i=0; i<_bPolyObject->getOrder(); i++)
  {
    yout        += (*b) * _doubleShift[shift_ptr];
    b++;
    shift_ptr++;
    if(shift_ptr == _shiftSize)
      shift_ptr = 0;
  }
  yout      *= gain;
//
// Add the new data to the shift register
//
  _shiftPtr--;
  if(_shiftPtr < 0)
    _shiftPtr                   = MAX((_shiftSize - 1), 0);
  _doubleShift[_shiftPtr]       = input;

  return yout;
}

// ############################# Public Method ###############################
// filterDoubleData -- Filter an input data sample, and return the sampled value.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   input:          An input data sample
//          
// Output:                  The filtered data sample
// ############################# Public Method ###############################
double AbstractFIR::filterDoubleData(double input)
{
  int           i, shift_ptr;
  double        yout, gain;
  const         double *b;

  if(_bPolyObject == NULL)
     return 0.;
   
  if(_passBandGain != 0.)
     gain       = 1./_passBandGain;
  else
     gain       = 1.;
//
// Now, get the coefficients
//
  b         = _bPolyObject->getCoefficients();
//
// Filter the data:
//
  yout          = (*b) * input;
  b++;
  shift_ptr     = _shiftPtr;
  for(i=0; i<_bPolyObject->getOrder(); i++)
  {
    yout        += (*b) * _doubleShift[shift_ptr];
    b++;
    shift_ptr++;
    if(shift_ptr == _shiftSize)
      shift_ptr = 0;
  }
  yout      *= gain;
//
// Add the new data to the shift register
//
  _shiftPtr--;
  if(_shiftPtr < 0)
    _shiftPtr                   = MAX((_shiftSize - 1), 0);
  _doubleShift[_shiftPtr]       = input;

  return yout;
}

// ############################# Public Method ###############################
// filterComplexArray -- Filter the input array and overwrite it with output.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   x:              An array of floating point data, modified on output
//          numberPts:      Number of points in the input array.
//          
// Output:                  YES if went OK, No if _bPolyObject doesn't exist.
// ############################# Public Method ###############################
LOGICAL AbstractFIR::filterComplexArray(Complex *x, int numberPts)
{
  int           j, i, order, shift_ptr;
  double        gain;
  const         double *b, *b_ptr;
  Complex       yout;

   if(_bPolyObject == NULL)
     return NO;
   
   if(_passBandGain != 0.)
     gain   = 1./_passBandGain;
   else
     gain   = 1.;
//
// Now, get the coefficients
//
  b_ptr         = _bPolyObject->getCoefficients();
  order         = _bPolyObject->getOrder();
//
// Filter the data:
//
  for(j=0; j<numberPts; j++)
  {
    b           = b_ptr;
    yout        = (*b) * x[j];
    b++;
    shift_ptr   = _shiftPtr;
    for(i=0; i<order; i++)
    {
      yout      += (*b) * _complexShift[shift_ptr];
      b++;
      shift_ptr++;
      if(shift_ptr == _shiftSize)
        shift_ptr       = 0;
    }
//
// Add the new data to the shift register
//
    _shiftPtr--;
    if(_shiftPtr < 0)
      _shiftPtr                 = MAX((_shiftSize - 1), 0);
    _complexShift[_shiftPtr]    = x[j];
//
// Save the output:
//
    x[j]                        = gain*yout;
  }

  return YES;
}

// ############################# Public Method ###############################
// filterComplexData -- Filter a complex data sample, and return the sampled value.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   xin:            An input data sample
//          
// Output:                  The filtered data sample
// ############################# Public Method ###############################
Complex  AbstractFIR::filterComplexData(Complex &input)
{
  int           i, shift_ptr;
  double        gain;
  const         double *b;
  Complex       yout;

  if(_bPolyObject == NULL)
     return 0.;

  if(_passBandGain != 0.)
     gain       = 1./_passBandGain;
  else
     gain       = 1.;
//
// Now, get the coefficients
//
  b         = _bPolyObject->getCoefficients();
//
// Filter the data:
//
  yout          = (*b) * input;
  b++;
  shift_ptr     = _shiftPtr;
  for(i=0; i<_bPolyObject->getOrder(); i++)
  {
    yout        += (*b) * _complexShift[shift_ptr];
    b++;
    shift_ptr++;
    if(shift_ptr == _shiftSize)
      shift_ptr = 0;
  }
  yout      *= gain;
//
// Add the new data to the shift register
//
  _shiftPtr--;
  if(_shiftPtr < 0)
    _shiftPtr                   = MAX((_shiftSize - 1), 0);
  _complexShift[_shiftPtr]      = input;

  return yout;
}


// ############################# Public Method ###############################
// complexFilterComplexArray -- Filter an array of complex data samples using
//                              complex coefficients.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:       coeff:          The set of complex coefficients
//              taps:           The number of coefficients
//              x:              An array of complex data, filtered on output
//              numberPts:      # of points in above array.
//          
// Output:                      None
// ############################# Public Method ###############################
void  AbstractFIR::complexFilterComplexArray(Complex *coeff, int taps, Complex *x, int numberPts)
{
  int           j, i, order, shift_ptr;
  Complex       *b_ptr, *b;
  Complex       yout;
//
// Filter the data:
//
  b_ptr         = coeff;
  order         = taps-1;
  for(j=0; j<numberPts; j++)
  {
    b           = b_ptr;
    yout        = (*b) * x[j];
    b++;
    shift_ptr   = _shiftPtr;
    for(i=0; i<order; i++)
    {
      yout      += (*b) * _complexShift[shift_ptr];
      b++;
      shift_ptr++;
      if(shift_ptr == _shiftSize)
        shift_ptr       = 0;
    }
//
// Add the new data to the shift register
//
    _shiftPtr--;
    if(_shiftPtr < 0)
      _shiftPtr                 = MAX((_shiftSize - 1), 0);
    _complexShift[_shiftPtr]    = x[j];
//
// Save the output:
//
    x[j]                        = yout;
  }

  return;
}

// ############################# Public Method ###############################
// filterResponseAtFrequencies -- This routine calculates the response of the filter.
// Input:       frequencies:    An array of frequencies given in rad/s
//              numberPts:      Number of frequencies
//              calc:           YES = re-calculate poles and zeros
//
// Output:                      An array of filter responses at the input frequencies
// Notes:
// ############################# Public Method ###############################
float *AbstractFIR::filterResponseAtFrequencies(float *frequencies, int numberPts, LOGICAL calc)
{
  int     n, m;
  Complex complex_response, carg, exponent, alias_response;
//
// First calculate filter coefficients and/or poles and zeros:
//
  if(calc)
    findTransferFunction();                             // This needs to be overridden in subclasses
//
// deltaOmega is used in outputResponseFor() for calculating group delay:
// _subAngle is used for PI crossings.
//
  if(numberPts>1)
    _deltaOmega = (frequencies[numberPts-1] - frequencies[0])/numberPts;
  _deltaOmega   *= _fs;
  _deltaOmega   = MAX(1e-30, _deltaOmega);
  _subAngle     = _thetaOld     = 0.;
    
/***************************
 * Get space for output    *
 * arrays:                 *
 ***************************/  
  if(numberPts > _oldNumberFrequencies)
  {
    _oldNumberFrequencies = numberPts;
    if(_filterResponse != NULL)
      delete [] _filterResponse;
    _filterResponse  = new float[numberPts];
  }
/********************************
 * Loop over frequencies:       *
 * Start the loop at less than  *
 * zero in order to initialize  *
 * phase/group delay.           *
 ********************************/
  setPassBandGain();
  for(n=-2; n<numberPts; n++)
  {
    m                           = MAX(0, n);
    if( _responseType == ALIAS_RESPONSE )
    {
      alias_response            = aliasResponseAt(frequencies[m]);
      _filterResponse[m]        = abs(alias_response)/_passBandGain;
    }
    else if (_responseType == DB_ALIAS_RESPONSE)
    {
      alias_response            = aliasResponseAt(frequencies[m]);
      _filterResponse[m]        = abs(alias_response)/_passBandGain;
      if(_filterResponse[m] > 0.)
        _filterResponse[m]      = 20.*log10(_filterResponse[m]);
    }
    else
    {
      exponent                  = Complex(0.,frequencies[m]);
      carg                      = exp(exponent);
      complex_response          = transferResponseAt(carg);
      _filterResponse[m]        = outputResponseFor(complex_response, _fs*frequencies[m]);
    }
  }
  return _filterResponse;
}

// ############################# Public Method ###############################
// filterResponseAtFrequency -- This routine calculates the response of the filter
//                               at a single input frequency
// Input:       frequency:      digital frequency given in rad/s
//              calc:           YES = re-calculate poles and zeros
//
// Output:                      Response at the input frequency
// Notes:
// ############################# Public Method ###############################
float AbstractFIR::filterResponseAtFrequency(float frequency, LOGICAL calc)
{
  float     filter_response;
  Complex   complex_response, carg, exponent, alias_response;
  
//
// First calculate filter coefficients and/or poles and zeros:
//
  if(calc)
    findTransferFunction();                             // This needs to be overridden in subclasses

/********************************
 * Calculate response at        *
 * given frequency:             *
 ********************************/
  setPassBandGain();
  if( _responseType == ALIAS_RESPONSE )
  {
      alias_response            = aliasResponseAt(frequency);
      filter_response           = abs(alias_response)/_passBandGain;
  }
  else if (_responseType == DB_ALIAS_RESPONSE)
  {
      alias_response            = aliasResponseAt(frequency);
      filter_response           = abs(alias_response)/_passBandGain;
      if(filter_response > 0.)
        filter_response         = 20.*log10(filter_response);
  }
  else
  {
    exponent            = Complex(0., frequency);
    carg                = exp(exponent);
    complex_response    = transferResponseAt(carg);
    filter_response     = outputResponseFor(complex_response, _fs*frequency);
  }

  return filter_response;
}

// ############################# Public Method ###############################
// complexResponseAtFrequency -- This routine calculates the response of the filter
//                               at a single input frequency
// Input:       frequency:      digital frequency given in rad/s
//              calc:           YES = re-calculate poles and zeros
//
// Output:                      Response at the input frequency
// Notes:
// ############################# Public Method ###############################
Complex AbstractFIR::complexResponseAtFrequency(float frequency, LOGICAL calc)
{
  Complex   complex_response, carg, exponent;
  
//
// First calculate filter coefficients and/or poles and zeros:
//
  if(calc)
    findTransferFunction();                             // This needs to be overridden in subclasses

  /********************************
 * Calculate response at freq:  *
 ********************************/
  if( _responseType == ALIAS_RESPONSE )
  {
      complex_response          = aliasResponseAt(frequency);
  }
  else
  {
    exponent                    = Complex(0., frequency);
    carg                        = exp(exponent);
    complex_response            = transferResponseAt(carg);
  }

  return complex_response;
}

// ############################# Public Method ###############################
// findTransferFunction -- Find the filter's transfer function
//
// Input:                   None
//          
// Output:                  None
//
// Notes:
//  1. This needs to be overridden in subclasses.
// ############################# Public Method ###############################
void AbstractFIR::findTransferFunction()
{
  zeroOutTaps();
  return;
}

// ############################# Public Method ###############################
// findPower -- Find and return the normalized power out of the filter.
//
// Input:               None
//          
// Output:              noise power in Watts
//
// Notes:
//  1. The transfer function must first be calculated.
//  2. The noise power is given by, P = (1/2PI)*Integral(-PI,PI, [H(jw)]**2),
//     which for an FIR filter is equal to the sum of the coefficients squared.
// ############################# Public Method ###############################
double AbstractFIR::findPower()
{
  int           i;
  double        f_power;
  const double *b_coeffs;
//
// Sum the squares of the coefficients:
//
  b_coeffs      = bCoeffs();
  f_power       = 0.;
  for(i=0; i<_numberTaps; i++)
    f_power     += b_coeffs[i]*b_coeffs[i];
  return f_power;
}

