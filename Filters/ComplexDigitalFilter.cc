/************************************************************************
 *                                                                      *
 * This subclass of ComplexFilter implements a digital filtering        *
 * object.  Note that this class has some capability for designing      *
 * its own coefficients, otherwise, they must be calculated externally, *
 * and loaded in.                                                       *
 *                                                                      *
 * File:ComplexDigitalFilter.cc                                         *
 *                                                                      *
 * The filter is stored in two forms:                                   *
 *                                                                      *
 *           b[0] + b[1]z^(-1) + ... + b[n]z^-(n)                       *
 *  H(z)   = ------------------------------------                       *
 *           1    + a[1]z^(-1) + ... + a[n]z^-(n)                       *
 *                                                                      *
 *  and:                                                                *
 *           (z-zero[0]) * (z-zero[1]) ... (z-zero[n_zero])             *
 *    H(z) = ----------------------------------------------             *
 *           (z-pole[0]) * (z-pole[1]) ... (z-pole[n_pole])             *
 *                                                                      *
 * Where the a and b polynomials are stored in the polynomial objects.  *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 10/11/03  - Started.                                             *
 ************************************************************************/

#include "ComplexDigitalFilter.h"                                       // Object prototypes

#if defined(WIN32)
#include <Polynomials/ComplexPoly.h>
#include <GNU/Complex.h>
#include <C_Libraries/constants.h>
#else
#include "ComplexPoly.h"
#include "Complex.h"
#include "constants.h"
#endif

#include <math.h>
#include <stdio.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))


// ############################# Private Method ###############################
// checkShiftSize -- Allocates the arrays for the digital shift registers.
//
// Input:           size:       New shift size needed
//          
// Output:                      None
//
// Notes:
// ############################# Private Method ###############################
void ComplexDigitalFilter::checkShiftSize(int size)
{  
//
// Check if we are big enough:
//
  if(_shiftSize < size)
  {
    delete [] _complexShift;

    _shiftSize      = size;
    _complexShift   = new Complex[_shiftSize+1];
    zeroOutTaps();                              // Clear the filter's taps
  }
  return;
}

// ############################# Private Method ###############################
// filterNextDouble -- Filter a single complex input variable.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   input:          Next complex point input
//          a:              denominator coefficients
//          b:              numerator coefficients
//          a_order:        order of denominator coefficients.
//          b_order:        order of numerator coefficients.
//          
// Output:                  YES if went OK, No if _aPolyObject doesn't exist.
// ############################# Private Method ###############################
Complex ComplexDigitalFilter::filterNextDouble(double input, const Complex *a, const Complex *b,
                                                int a_order, int b_order)
{
  int       i, i_minus_1, i_minus_2;
  Complex   yout, complex_input;

  complex_input = Complex(input, 0.);

//
// First consider case when a_order == b_order:
//
  if(a_order == b_order)
  {
//
// Check for boundary case:
//
    if( a_order == 0)
    {
       return complex_input*b[0];
    }
    yout        = b[1]*_complexShift[0];
    for(i=a_order; i>1 ;i--)
    {
       i_minus_1                = i-1;
       yout                     += b[i]*_complexShift[i_minus_1];
       complex_input            -= a[i]*_complexShift[i_minus_1];
       _complexShift[i_minus_1] = _complexShift[i-2];
    }
    complex_input               -= a[1]*_complexShift[0];
    _complexShift[0]            = complex_input;
    yout                        += b[0]*complex_input;
  }
//
// Else a_order < b_order
//
  else if(a_order < b_order)
  {
//
// Check for boundary case:
//
    if( b_order == 0)
    {
       return complex_input*b[0];
    }
//
// First consider part of delay line that doesn't include a[i]:
//
    yout        = 0.;
    for(i=b_order; i>a_order ;i--)
    {
       i_minus_1                = i-1;
       i_minus_2                = MAX(0, i-2);
       yout                     += b[i]*_complexShift[i_minus_1];
       _complexShift[i_minus_1] = _complexShift[i_minus_2];
    }
//
// Now, consider both a[i] and b[i]:
//
    for(i=a_order; i>1; i--)
    {
       i_minus_1                = i-1;
       yout                     += b[i]*_complexShift[i_minus_1];
       complex_input            -= a[i]*_complexShift[i_minus_1];
       _complexShift[i_minus_1] = _complexShift[i-2];
    }
    if(a_order > 0)
      complex_input             -= a[1]*_complexShift[0];
    _complexShift[0]            = complex_input;
    yout                        += b[0]*complex_input;
  }
//
// Finally, b_order < a_order:
//
  else
  {
//
// Check for boundary case:
//
    if( a_order == 0)
    {
       return complex_input*b[0];
    }
    if(b_order > 0)
      yout                      = b[1]*_complexShift[0];
    else
      yout                      = Complex(0., 0.);
    for(i=a_order; i>1 ;i--)
    {
      i_minus_1                = i-1;
      if(i<=b_order)
        yout                     += b[i]*_complexShift[i_minus_1];
      complex_input             -= a[i]*_complexShift[i_minus_1];
      _complexShift[i_minus_1] = _complexShift[i-2];
    }
    complex_input               -= a[1]*_complexShift[0];
    _complexShift[0]            = complex_input;
    yout                        += b[0]*complex_input;
  }

  return yout;
}
// ############################# Private Method ###############################
// filterNextComplex -- Filter a single complex input variable.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   input:          Next complex point input
//          a:              denominator coefficients
//          b:              numerator coefficients
//          a_order:        order of denominator coefficients.
//          b_order:        order of numerator coefficients.
//          
// Output:                  YES if went OK, No if _aPolyObject doesn't exist.
// ############################# Private Method ###############################
Complex ComplexDigitalFilter::filterNextComplex(Complex input, const Complex *a, const Complex *b,
                                                int a_order, int b_order)
{
  int       i, i_minus_1, i_minus_2;
  Complex   yout;

//
// First consider case when a_order == b_order:
//
  if(a_order == b_order)
  {
//
// Check for boundary case:
//
    if( a_order == 0)
    {
       return input*b[0];
    }
    yout        = b[1]*_complexShift[0];
    for(i=a_order; i>1 ;i--)
    {
       i_minus_1                = i-1;
       yout                     += b[i]*_complexShift[i_minus_1];
       input                    -= a[i]*_complexShift[i_minus_1];
       _complexShift[i_minus_1] = _complexShift[i-2];
    }
    input                       -= a[1]*_complexShift[0];
    _complexShift[0]            = input;
    yout                        += b[0]*input;
  }
//
// Else a_order < b_order
//
  else if(a_order < b_order)
  {
//
// Check for boundary case:
//
    if( b_order == 0)
    {
       return input*b[0];
    }
//
// First consider part of delay line that doesn't include a[i]:
//
    yout        = 0.;
    for(i=b_order; i>a_order ;i--)
    {
       i_minus_1                = i-1;
       i_minus_2                = MAX(0, i-2);
       yout                     += b[i]*_complexShift[i_minus_1];
       _complexShift[i_minus_1] = _complexShift[i_minus_2];
    }
//
// Now, consider both a[i] and b[i]:
//
    for(i=a_order; i>1; i--)
    {
       i_minus_1                = i-1;
       yout                     += b[i]*_complexShift[i_minus_1];
       input                    -= a[i]*_complexShift[i_minus_1];
       _complexShift[i_minus_1] = _complexShift[i-2];
    }
    if(a_order > 0)
      input                     -= a[1]*_complexShift[0];
    _complexShift[0]            = input;
    yout                        += b[0]*input;
  }
//
// Finally, b_order < a_order:
//
  else
  {
//
// Check for boundary case:
//
    if( a_order == 0)
    {
       return input*b[0];
    }
    if(b_order > 0)
      yout                      = b[1]*_complexShift[0];
    else
      yout                      = Complex(0., 0.);
    for(i=a_order; i>1 ;i--)
    {
      i_minus_1                = i-1;
      if(i<=b_order)
        yout                     += b[i]*_complexShift[i_minus_1];
      input                    -= a[i]*_complexShift[i_minus_1];
      _complexShift[i_minus_1] = _complexShift[i-2];
    }
    input                       -= a[1]*_complexShift[0];
    _complexShift[0]            = input;
    yout                        += b[0]*input;
  }

  return yout;
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
float ComplexDigitalFilter::outputResponseFor(Complex cResponse , float omega)
{
  float output_response, center_of_mass;

  switch(_responseType)
  {
    default:
//
// Call super's method
//
      output_response   = AbstractFilter::outputResponseFor(cResponse, omega);
      break;
    case PHASE_NORMALIZED:                      // This is overridden from super
      center_of_mass    = denominatorCenterOfMass();
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
      output_response   -= omega*center_of_mass/_fs;            // Only makes sense for FIR filters
      if(_phaseInDegrees)
        output_response *= RAD_DEG;
      break;
  }
  return output_response;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the ComplexDigitalFilter class.
//
// Input:           aCoeffs:    The coefficients of the denominator polynomial
//                  bCoeffs:    The coefficients of the numerator polynomial
//                  order:      The order of the polynomials
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
ComplexDigitalFilter::ComplexDigitalFilter(Complex *aCoeffs, Complex *bCoeffs, int order)
                     :ComplexFilter(aCoeffs, bCoeffs, order)
{
  _shiftSize            = -1;
  _complexShift         = NULL;
  _fs                   = 1.;  
  
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the ComplexDigitalFilter class.
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
ComplexDigitalFilter::ComplexDigitalFilter(int pass, double centerFreq, double cutoffFreq,
                                  double samplingFreq, int order)
              :ComplexFilter(pass, centerFreq, cutoffFreq, order)
{
  _shiftSize            = -1;
  _complexShift         = NULL;
  _fs                   = samplingFreq;  
  
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the ComplexDigitalFilter class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
ComplexDigitalFilter::~ComplexDigitalFilter()
{
  delete [] _complexShift;

  return;
}

// ############################# Public Method ###############################
// filterFloatData -- Filter an input data sample, and return the sampled value.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   input:          An input data sample
//          
// Output:                  The filtered data sample
// ############################# Public Method ###############################
Complex ComplexDigitalFilter::filterFloatData(float input)
{
  double        gain;
  Complex       yout;
  const         Complex *a, *b;

   if(_aPolyObject == NULL)
     return 0.;
   
   if(_passBandGain != 0.)
     gain       = 1./_passBandGain;
   else
     gain       = 1.;
//
// Now, get the coefficients
//
  a         = _aPolyObject->getCoefficients();
  b         = _bPolyObject->getCoefficients();
//
// Filter the data:
//
  yout      = filterNextDouble(input, a, b, _aPolyObject->getOrder(), _bPolyObject->getOrder());
  yout      *= gain;
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
Complex ComplexDigitalFilter::filterDoubleData(double input)
{
  double        gain;
  Complex       yout;
  const         Complex *a, *b;

   if(_aPolyObject == NULL)
     return 0.;
   
   if(_passBandGain != 0.)
     gain       = 1./_passBandGain;
   else
     gain       = 1.;
//
// Now, get the coefficients
//
  a         = _aPolyObject->getCoefficients();
  b         = _bPolyObject->getCoefficients();
//
// Filter the data:
//
  yout      = filterNextDouble(input, a, b, _aPolyObject->getOrder(), _bPolyObject->getOrder());
  yout      *= gain;
  return yout;
}

// ############################# Public Method ###############################
// filterComplexArray -- Filter the input array and overwrite it with output.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   x:              An array of floating point data, modified on output
//          numberPts:      Number of points in the input array.
//          
// Output:                  YES if went OK, No if _aPolyObject doesn't exist.
// ############################# Public Method ###############################
LOGICAL ComplexDigitalFilter::filterComplexArray(Complex *x, int numberPts)
{
  int       j;
  int       a_order, b_order;
  double    gain;
  const     Complex *a, *b;

   if(_aPolyObject == NULL)
     return NO;
   
   if(_passBandGain != 0.)
     gain   = 1./_passBandGain;
   else
     gain   = 1.;
//
// Now, get the coefficients
//
   a        = _aPolyObject->getCoefficients();
   b        = _bPolyObject->getCoefficients();
   a_order  = _aPolyObject->getOrder();
   b_order  = _bPolyObject->getOrder();
      
   for(j=0; j<numberPts; j++)
     x[j]   = gain*filterNextComplex(x[j], a, b, a_order, b_order);
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
Complex  ComplexDigitalFilter::filterComplexData(Complex &input)
{
  double    gain;
  const     Complex *a, *b;
  Complex   yout;

   if(_aPolyObject == NULL)
     return 0.;
   
   if(_passBandGain != 0.)
     gain       = 1./_passBandGain;
   else
     gain       = 1.;
//
// Now, get the coefficients
//
  a         = _aPolyObject->getCoefficients();
  b         = _bPolyObject->getCoefficients();
//
// Filter the data:
//
  yout      = filterNextComplex(input, a, b, _aPolyObject->getOrder(), _bPolyObject->getOrder());
  yout      *= gain;
  return yout;
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
void ComplexDigitalFilter::setPassBandGain(LOGICAL usePolesZeros)
{
  Complex   arg, num, den, response;
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
    den                 = _aPolyObject->evaluateAtComplexPoint(arg);
    if(abs(den) > 0.)
      _passBandGain     = abs(num/den);
    else
      _passBandGain     = 1.;
  }
  return;
}

// ############################# Public Method ###############################
// zeroOutTaps -- checks the sizer of the shift register, and zeros out the
//                digital filter's taps.
// Input:                       None
//          
// Output:                      None
// ############################# Public Method ###############################
void ComplexDigitalFilter::zeroOutTaps()
{
  int i, max_shift;
//
// Check the size of the shift register:
//
  max_shift     = MAX(_aPolyObject->getOrder(), _bPolyObject->getOrder());
//
// Check if we are big enough:
//
  if(_shiftSize < max_shift)
  {
    delete [] _complexShift;

    _shiftSize      = max_shift;
    _complexShift   = new Complex[_shiftSize+1];
  }

  for(i=0; i<_shiftSize; i++)
  {
    _complexShift[i]    = Complex(0., 0.);
  }
  return;
}

// ############################# Public Method ###############################
// setFilterACoeffs -- Sets filter denominator polynomials
// Input:   aCoeff:         coefficients of denominator polynomial
//          
// Output:                  None
// ############################# Public Method ###############################
void ComplexDigitalFilter::setFilterACoeffs(const Complex *aCoeff)
{
  int   max_shift;
//
// Perform super's function
//
  ComplexFilter::setFilterACoeffs(aCoeff);
//
// Add this
//
// Check the shift register:
//
  if((_aPolyObject!=NULL) && (_bPolyObject!=NULL) )
  {
    max_shift     = MAX(_aPolyObject->getOrder(), _bPolyObject->getOrder());
    checkShiftSize(max_shift);
    zeroOutTaps();                              // Clear the filter's taps
  }
  return;
}

// ############################# Public Method ###############################
// setFilterBCoeffs -- Sets filter numerator polynomials
// Input:   bCoeff:         coefficients of numerator polynomial
//          
// Output:                  None
// ############################# Public Method ###############################
void ComplexDigitalFilter::setFilterBCoeffs(const Complex *bCoeff)
{
  int   max_shift;
//
// Perform super's function
//
  ComplexFilter::setFilterBCoeffs(bCoeff);
//
// Add this
//
// Check the shift register:
//
  if((_aPolyObject!=NULL) && (_bPolyObject!=NULL) )
  {
    max_shift     = MAX(_aPolyObject->getOrder(), _bPolyObject->getOrder());
    checkShiftSize(max_shift);
    zeroOutTaps();                              // Clear the filter's taps
  }
  return;
}

// ############################# Public Method ###############################
// scaleCoefficientsByPassBandGain -- Scales the filter's numerator (b) coefficients
//                                    by the pass band gain
// Input:               None
//          
// Output:              None
//
// Notes
// 1. The instance variables _bPolyObject, and _passBandGain are modified.
// ############################# Public Method ###############################
void ComplexDigitalFilter::scaleCoefficientsByPassBandGain()
{
//
// Calculate the pass band gain based on filter pass type
//
  setPassBandGain();
  if(_passBandGain != 0.)
  {
    *_bPolyObject       /= Complex(_passBandGain, 0.);
    _passBandGain       = 1.;
  }
  return;
}

// ############################# Public Method ###############################
// filterPoles -- Returns the complex poles of the digital filter.
// Input:                       None
//          
// Output:                      An array of complex pole positions
// Notes:
//  1. To find the size of the array, call numberOfPoles().
// ############################# Public Method ###############################
const Complex *ComplexDigitalFilter::filterPoles()
{
  int       i;
  Complex   *roots;

  findTransferFunction();                               // calculate the poles and zeros or xfer fn
  if(_aPolyObject != NULL)
  {
        delete []       _filterPoles;                           // Find poles and zeros from xfer function
        _numberPoles    = _aPolyObject->getOrder();
        _filterPoles    = new Complex[_numberPoles];
        roots           = _aPolyObject->roots();
        for(i=0; i<_numberPoles; i++)
          _filterPoles[i]  = roots[i];
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
// ############################# Public Method ###############################
const Complex *ComplexDigitalFilter::filterZeros()
{
  int       i;
  Complex   *roots;

  findTransferFunction();                                       // calculate the poles and zeros or xfer fn
  if(_bPolyObject != NULL)
  {
        delete []       _filterZeros;
        _numberZeros    = _bPolyObject->getOrder();
        _filterZeros    = new Complex[_numberZeros];
        roots           = _bPolyObject->roots();
        for(i=0; i<_numberZeros; i++)
        {
          _filterZeros[i]       = roots[i];
        }
  }
  return _filterZeros;
}

// ############################# Public Method ###############################
// filterResponseAtFrequencies -- This routine calculates the response of the filter.
// Input:       frequencies:    An array of frequencies given in rad/s
//              numberPts:      Number of frequencies
//
// Output:                      An array of filter responses at the input frequencies
// Notes:
// ############################# Public Method ###############################
float *ComplexDigitalFilter::filterResponseAtFrequencies(float *frequencies, int numberPts)
{
  int     n, m;
  Complex complex_response, carg, exponent;
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
    _oldNumberFrequencies       = numberPts;
    if(_filterResponse != NULL)
      delete [] _filterResponse;
    _filterResponse             = new float[numberPts];
  }
/********************************
 * Loop over frequencies:       *
 * Start the loop at less than  *
 * zero in order to initialize  *
 * phase/group delay.           *
 ********************************/
  for(n=-2; n<numberPts; n++)
  {
    m                   = MAX(0, n);
    exponent            = Complex(0.,frequencies[m]);
    carg                = exp(exponent);
    if(_filterStructureType == POLE_ZERO)
      complex_response  = poleZeroResponseAt(carg, _filterPoles, _filterZeros);
    else 
      complex_response  = transferResponseAt(carg);
    
    _filterResponse[m]   = outputResponseFor(complex_response, _fs*frequencies[m]);
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
float ComplexDigitalFilter::filterResponseAtFrequency(float frequency)
{
  float     filter_response;
  Complex   complex_response, carg, exponent;

/********************************
 * Calculate response at        *
 * given frequency:             *
 ********************************/
  exponent              = Complex(0., frequency);
  carg                  = exp(exponent);
  if(_filterStructureType == POLE_ZERO)
      complex_response  = poleZeroResponseAt(carg, _filterPoles, _filterZeros);
  else 
      complex_response  = transferResponseAt(carg);
    
  filter_response       = outputResponseFor(complex_response, _fs*frequency);

  return filter_response;
}

// ############################# Public Method ###############################
// complexResponseAtFrequency -- This routine calculates the response of the filter
//                               at a single input frequency
// Input:       frequency:      digital frequency given in rad/s
//
// Output:                      Response at the input frequency
// Notes:
// ############################# Public Method ###############################
Complex ComplexDigitalFilter::complexResponseAtFrequency(float frequency)
{
  Complex   complex_response, carg, exponent;

  /********************************
 * Calculate response at freq:  *
 ********************************/
  exponent          = Complex(0., frequency);
  carg              = exp(exponent);
  if(_filterStructureType == POLE_ZERO)
      complex_response  = poleZeroResponseAt(carg, _filterPoles, _filterZeros);
  else 
      complex_response  = transferResponseAt(carg);

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
void ComplexDigitalFilter::findTransferFunction()
{
  zeroOutTaps();
  return;
}


// ############################# Public Method ###############################
// shiftCoefficientsByFrequency -- Shift the filter's coefficients to re-center
//                                 the filter to a new frequency.
//
// Input:       analogFrequency:        Analog frequency in Hertz
//          
// Output:                              None
//
// Notes:
//  1. This routine could be used to generate a complex (asymmetrical about DC)
//     bandpass filter from a low pass filter.  For example, load the (real) filter
//     coefficients into the a and b polynomials, set the sampling frequency, and
//     then call this function to get complex coefficients.
//  2. See Oppenheim and Schafer, section 3.4.3
// ############################# Public Method ###############################
void ComplexDigitalFilter::shiftCoefficientsByFrequency(float analogFrequency)
{
  int           i, a_size, b_size;
  double        omega_0;
  Complex       exponential, exponential_to_i;
  const         Complex *old_a;
  const         Complex *old_b;
  Complex       *new_a, *new_b;
// Find omega_0 in the digital domain
  if(_fs != 0.)
    omega_0     = TWOPI*analogFrequency/_fs;
  else
    omega_0     = 1.;
// Find exp(j*omeg_0);
  exponential   = Complex(0., omega_0);
  exponential   = exp(exponential);
// Get the old coefficients, and allocate new
  old_a         = aCoeffs();
  old_b         = bCoeffs();
  a_size        = aOrder()+1;
  b_size        = bOrder()+1;
  new_a         = new Complex[a_size];
  new_b         = new Complex[b_size];
//
// Scale the new coefficients by the exponential:
//
  exponential_to_i      = Complex(1., 0.);
  for(i=0; i<a_size; i++)
  {
    new_a[i]            = exponential_to_i*old_a[i];
    exponential_to_i    *= exponential;
  }
  for(i=0; i<b_size; i++)
  {
    new_b[i]            = exponential_to_i*old_b[i];
    exponential_to_i    *= exponential;
  }
//
// Set the new coefficients:
//
  setFilterACoeffs(new_a);
  setFilterBCoeffs(new_b);
  delete [] new_a;
  delete [] new_b;
  return;
}

// ############################# Public Method ###############################
// convertCoefficientsToZ -- Convert coefficients previously loaded using
//                           setFilterACoeffs() or setFilterBCoeffs() to z
//                           domain using s -> 2(z-1)/(z+1)/T (bilinear)
//
// Input:                   None
//          
// Output:                  None
//
// Notes:
//  1. The instanceVariables, _aPolyObject and _bPolyObject are modified.
//  2. The s domain poly is of the form: p(s) = b[n] + b[n-1]s + ... + b[0]s**(n)
//
// ############################# Public Method ###############################
void ComplexDigitalFilter::convertCoefficientsToZ()
{
  int           i, j, power;
  int           a_order, b_order;
  double        ts;
  const         Complex *a_coeffs, *b_coeffs;
  Complex       a0;
  Complex       tran_coeff[2];                  // 2z-2
  Complex       den_coeff[2];                   // z+1
  ComplexPoly   sum_poly, trans_poly, temp_poly;
  ComplexPoly   den_poly, mult_poly;
//
// set the transform polynomial:
//
  tran_coeff[0] = Complex(2., 0.);
  tran_coeff[1] = Complex(-2., 0.);
  den_coeff[0]  = Complex(1., 0.);
  den_coeff[1]  = Complex(1., 0.);
  trans_poly.assign(1, tran_coeff);
  den_poly.assign(1, den_coeff);
  mult_poly     = 1.;
//
// First convert a polynomial:
//
  if(_fs > 0.)
    ts          = 1./_fs;
  else
    ts          = 1.;
  a_order       = _aPolyObject->getOrder();
  a_coeffs      = _aPolyObject->getCoefficients();
  for(i=0; i<a_order; i++)
  {
    temp_poly   = trans_poly;
    temp_poly   *= a_coeffs[i];
    power       = a_order-i-1;
    for(j=0; j<power; j++)
      temp_poly *= trans_poly;
//
// Mult by [T*(z+1)]^i
//
    temp_poly   *= pow(ts, i);
    temp_poly   *= mult_poly;
    mult_poly   *= den_poly;
//
// Add term to sum:
//
    if(i==0)
      sum_poly  = temp_poly;
    else
      sum_poly  += temp_poly;
  }
//
// Add the constant term, and copy to A polynomial
//
  temp_poly     = a_coeffs[a_order]*pow(ts,a_order);
  temp_poly     *= mult_poly;
  sum_poly      += temp_poly;
  *_aPolyObject = sum_poly;
//
// Now convert b polynomial:
//
  mult_poly     = Complex(1., 0.);
  b_order       = _bPolyObject->getOrder();
  b_coeffs      = _bPolyObject->getCoefficients();
  for(i=0; i<b_order; i++)
  {
    temp_poly   = trans_poly;
    temp_poly   *= b_coeffs[i];
    power       = b_order-i-1;
    for(j=0; j<power; j++)
      temp_poly *= trans_poly;
//
// Mult by [T*(z+1)]^i
//
    temp_poly   *= pow(ts, i);
    temp_poly   *= mult_poly;
    mult_poly   *= den_poly;
//
// Add term to sum:
//
    if(i==0)
      sum_poly  = temp_poly;
    else
      sum_poly  += temp_poly;
  }
//
// Add the constant term and copy to B polynomial:
//
  temp_poly     = b_coeffs[b_order]*pow(ts,b_order);
  temp_poly     *= mult_poly;
  sum_poly      += temp_poly;
  *_bPolyObject = sum_poly;
//
// Normalize coefficients so a[0] = 1.
//
  a_coeffs      = _aPolyObject->getCoefficients();
  a0            = a_coeffs[0];
  if(abs(a0) != 0.)
  {
    *_aPolyObject       /= a0;
    *_bPolyObject       /= a0;
  }
  
  return;
}


// ############################# Public Method ###############################
// convertCoefficientsToZ -- Convert coefficients previously loaded using
//                           setFilterACoeffs() or setFilterBCoeffs() to z
//                           domain using s->(z-1)/zT
//
// Input:                   None
//          
// Output:                  None
//
// Notes:
//  1. The instanceVariables, _aPolyObject and _bPolyObject are modified.
//  2. The s domain poly is of the form: p(s) = b[n] + b[n-1]s + ... + b[0]s**(n)
//
// ############################# Public Method ###############################
void ComplexDigitalFilter::convertCoefficientsToZUsingBackwardDifference()
{
  int           i, j, power;
  int           a_order, b_order;
  double        ts;
  const         Complex *a_coeffs, *b_coeffs;
  Complex       a0;
  Complex       poly_coeff[2];                  // z-1
  ComplexPoly   sum_poly, trans_poly, temp_poly;
  ComplexPoly   den_poly, mult_poly;
//
// set the transform polynomial:
//
  poly_coeff[0] = Complex(1.,0.);
  poly_coeff[1] = Complex(-1.,0.);
  trans_poly.assign(1, poly_coeff);
//
// First convert a polynomial:
//
  if(_fs > 0.)
    ts          = 1./_fs;
  else
    ts          = 1.;
  a_order       = _aPolyObject->getOrder();
  a_coeffs      = _aPolyObject->getCoefficients();
  for(i=0; i<a_order; i++)
  {
    temp_poly   = trans_poly;
    temp_poly   *= a_coeffs[i];
    power       = a_order-i-1;
    for(j=0; j<power; j++)
      temp_poly *= trans_poly;
//
// Mult by (Tz)^i
//
    temp_poly   *= pow(ts, i);
    temp_poly   <<= i;
//
// Add term to sum:
//
    if(i==0)
      sum_poly  = temp_poly;
    else
      sum_poly  += temp_poly;
  }
//
// Add the constant term, and copy to A polynomial
//
  temp_poly     = a_coeffs[a_order]*pow(ts,a_order);
  temp_poly     <<= a_order;
  sum_poly      += temp_poly;
  *_aPolyObject = sum_poly;
//
// Now convert b polynomial:
//
  b_order       = _bPolyObject->getOrder();
  b_coeffs      = _bPolyObject->getCoefficients();
  for(i=0; i<b_order; i++)
  {
    temp_poly   = trans_poly;
    temp_poly   *= b_coeffs[i];
    power       = b_order-i-1;
    for(j=0; j<power; j++)
      temp_poly *= trans_poly;
//
// Mult by (Tz)^i
//
    temp_poly   *= pow(ts, i);
    temp_poly   <<= i;
//
// Add term to sum:
//
    if(i==0)
      sum_poly  = temp_poly;
    else
      sum_poly  += temp_poly;
  }
//
// Add the constant term and copy to B polynomial:
//
  temp_poly     = b_coeffs[b_order]*pow(ts,b_order);
  temp_poly     <<= b_order;
  sum_poly      += temp_poly;
  *_bPolyObject = sum_poly;
  return;
}

// ############################# Public Method ###############################
// warpFrequency -- warp the input frequency according to fs
//
// Input:       omega:          Analog frequency in radians/second
//              fs:             Sampling rate in Hertz
//
// Output:                      Warped frequency
//
// Notes:
// ############################# Public Method ###############################
double ComplexDigitalFilter::warpFrequency(double omega, double fs)
{
  double        omega_warp      = 1.;
  if(fs != 0.)
  {
    omega_warp  = omega/fs;
  }
  omega_warp    = 2.*fs*tan(omega_warp/2.);
  return omega_warp;
}


