/************************************************************************
 *                                                                      *
 * This subclass of RealFilter implements a digital filtering           *
 * object.  Note that this class is somewhat abstract, see the sub-     *
 * classes, FIRFilter and IIRFilter for actual implementations.         *
 * At the present time, it contains code for general purpose IIR        *
 * filtering, but can not calculate the filter coefficients.            *
 *                                                                      *
 * File:DigitalFilter.cc                                                *
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
 *  1. 12/30/94  - Objective-C -> C++                                   *
 *  2. 01/17/95  - Add filterDoubleData().                              *
 *  3. 01/27/95  - x.real() -> real(x), add ASINH_EPS.                  *
 *  4. 02/27/95  - Fixed destructor to account for ~AbstractFilter being*
 *                 called.  Set window -> NULL in constructor.          *
 *  5. 07/12/95  - Fixed new on analogPoles, analogZeros.               *
 *  6. 02/20/97  - Add correlateVectors().                              *
 *  7. 03/12/97  - Add setPassBandGain().                               *
 *  8. 03/13/97  - Add CIC filter type.                                 *
 *  9. 03/13/97  - Base filtering on polynomial orders.                 *
 * 10. 04/29/97  - Move windowing stuff to DataWindow class.            *
 * 11. 06/04/97  - Added remez routine.                                 *
 * 12. 07/05/97  - Added raised cosine stuff.                           *
 * 13. 07/14/97  - Add root raised cosine stuff.                        *
 * 14. 01/06/98  - Made superclass for IIRFilter and FIRFilter.         *
 * 15. 03/03/98  - Modify zeroOutTaps to check shift size.              *
 * 16. 12/18/00  - Add ALIAS_RESPONSE                                   *
 * 17. 03/27/02  - Add findPower.                                       *
 ************************************************************************/

#include "DigitalFilter.h"                                      // Object prototypes

#if defined(WIN32)
#include <GNU/Complex.h>
#include <C_Libraries/constants.h>
#include <Polynomials/DoublePoly.h>
#else
#include "Complex.h"
#include "constants.h"
#include "DoublePoly.h"
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
void DigitalFilter::checkShiftSize(int size)
{  
//
// Check if we are big enough:
//
  if(_shiftSize < size)
  {
    delete [] _complexShift;
    delete [] _doubleShift;

    _shiftSize      = size;
    _doubleShift    = new double[_shiftSize+1];
    _complexShift   = new Complex[_shiftSize+1];
    zeroOutTaps();                              // Clear the filter's taps
  }
  return;
}

// ############################# Private Method ###############################
// filterNextDouble -- Filter a single double input variable.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   input:          Next double point input
//          a:              denominator coefficients
//          b:              numerator coefficients
//          a_order:        order of denominator coefficients.
//          b_order:        order of numerator coefficients.
//          
// Output:                  YES if went OK, No if _aPolyObject doesn't exist.
// ############################# Private Method ###############################
double DigitalFilter::filterNextDouble(double input, const double *a, const double *b, int a_order,
                                      int b_order)
{
  int       i, i_minus_1, i_minus_2;
  double    yout;

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
    yout        = b[1]*_doubleShift[0];
    for(i=a_order; i>1 ;i--)
    {
       i_minus_1                = i-1;
       yout                     += b[i]*_doubleShift[i_minus_1];
       input                    -= a[i]*_doubleShift[i_minus_1];
       _doubleShift[i_minus_1]  = _doubleShift[i-2];
    }
    input                       -= a[1]*_doubleShift[0];
    _doubleShift[0]             = input;
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
       yout                     += b[i]*_doubleShift[i_minus_1];
       _doubleShift[i_minus_1]  = _doubleShift[i_minus_2];
    }
//
// Now, consider both a[i] and b[i]:
//
    for(i=a_order; i>1; i--)
    {
       i_minus_1                = i-1;
       yout                     += b[i]*_doubleShift[i_minus_1];
       input                    -= a[i]*_doubleShift[i_minus_1];
       _doubleShift[i_minus_1]  = _doubleShift[i-2];
    }
    if(a_order > 0)
      input                     -= a[1]*_doubleShift[0];
    _doubleShift[0]             = input;
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
      yout      = b[1]*_doubleShift[0];
    else
      yout      = 0.;
    for(i=a_order; i>1 ;i--)
    {
      i_minus_1                 = i-1;
      if(i <= b_order)
        yout                    += b[i]*_doubleShift[i_minus_1];
      input                     -= a[i]*_doubleShift[i_minus_1];
      _doubleShift[i_minus_1]   = _doubleShift[i-2];
    }
    input                       -= a[1]*_doubleShift[0];
    _doubleShift[0]             = input;
    yout                        += b[0]*input;
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
Complex DigitalFilter::filterNextComplex(Complex input, const double *a, const double *b, int a_order,
                                         int b_order)
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
      i_minus_1                 = i-1;
      if(i <= b_order)
        yout                    += b[i]*_complexShift[i_minus_1];
      input                     -= a[i]*_complexShift[i_minus_1];
      _complexShift[i_minus_1]  = _complexShift[i-2];
    }
    input                       -= a[1]*_complexShift[0];
    _complexShift[0]            = input;
    yout                += b[0]*input;
  }

  return yout;
}

// ############################# Private Method ###############################
// complexFilterNextComplex -- Filter a single complex input variable.
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
Complex DigitalFilter::complexFilterNextComplex(Complex input, const Complex *a, const Complex *b,
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
// aliasResponseAt -- Find the (complex magnitude) alias response at the input
//                    frequency value.
// Input:       frequency:      digital frequency given in rad/s
//          
// Output:                      Complex alias response
// NOTES:
// 1. This function sums the aliasing terms per (2.65a) in Crochiere and Rabiner
// ############################# Private Method ###############################
Complex DigitalFilter::aliasResponseAt(float frequency)
{
  int           j, old_type;
  float         alias_frequency;
  Complex       alias_response;

  old_type              = _responseType;
  _responseType         = MAGNITUDE;
  alias_response        = Complex(0.,0.);
  for(j=1; j<_decimateFactor; j++)
  {
    alias_frequency     = frequency - TWOPI*j/_decimateFactor;
    alias_response      += complexResponseAtFrequency(alias_frequency, NO);
  }
  _responseType         = old_type;
  return alias_response;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the DigitalFilter class.
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
DigitalFilter::DigitalFilter(int pass, double centerFreq, double cutoffFreq,
                                  double samplingFreq, int order)
              :RealFilter(pass, centerFreq, cutoffFreq, order)
{
  _shiftSize            = -1;
  _doubleShift          = NULL;
  _complexShift         = NULL;
  _fs                   = samplingFreq;
  _decimateFactor       = 1;
  
  return;
}
// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the DigitalFilter class.
//
// Input:           aCoeffs:    The coefficients of the denominator polynomial
//                  bCoeffs:    The coefficients of the numerator polynomial
//                  order:      The order of the polynomials
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
DigitalFilter::DigitalFilter(double *aCoeffs, double *bCoeffs, int order)
        :RealFilter(aCoeffs, bCoeffs, order)
{
  _shiftSize            = -1;
  _doubleShift          = NULL;
  _complexShift         = NULL;
  _fs                   = 1.;
  _decimateFactor       = 1;
  
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the DigitalFilter class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
DigitalFilter::~DigitalFilter()
{
  delete [] _doubleShift;
  delete [] _complexShift;

  return;
}

// ############################# Public Method ###############################
// filterFloatArray -- Filter the input array and overwrite it with output.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   x:              An array of floating point data, modified on output
//          numberPts:      Number of points in the input array.
//          
// Output:                  YES if went OK, No if _aPolyObject doesn't exist.
// ############################# Public Method ###############################
LOGICAL DigitalFilter::filterFloatArray(float * input, int numberPts)
{
  int       j, a_order, b_order;
  double    gain;
  const     double *a, *b;

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
     input[j]   = gain*filterNextDouble(input[j], a, b, a_order, b_order);
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
double DigitalFilter::filterFloatData(float input)
{
  double    yout, gain;
  const     double *a, *b;

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
double DigitalFilter::filterDoubleData(double input)
{
  double    yout, gain;
  const     double *a, *b;

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
LOGICAL DigitalFilter::filterComplexArray(Complex *x, int numberPts)
{
  int       j;
  int       a_order, b_order;
  double    gain;
  const     double *a, *b;

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
// complexFilterComplexArray -- Filter the input array and overwrite it with output.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   a:              The complex filter coefficients
//          b:              The complex filter numerator coefficients
//          n:              The order of the filter
//          x:              An array of floating point data, modified on output
//          numberPts:      Number of points in the input array.
//          
// Output:                  YES if went OK, No if _aPolyObject doesn't exist.
// ############################# Public Method ###############################
LOGICAL DigitalFilter::complexFilterComplexArray(Complex *a, Complex *b, int n, Complex *x, int numberPts)
{
  int j;

//
// Check the shift register:
//
  checkShiftSize(n);
//
// Now, filter:
//
  for(j=0; j<numberPts; j++)
     x[j]   = complexFilterNextComplex(x[j], a, b, n, n);

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
Complex  DigitalFilter::filterComplexData(Complex &input)
{
  double    gain;
  const     double *a, *b;
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
// setFilterParameters -- Sets new filter frequency parameters
// Input:   centerFreq:         Filter center frequency in Hertz
//          cutoff:             Filter cutoff frequency in Hertz
//          samplingFreq:       Filter sampling frequency in Hertz
//          
// Output:                      None
// ############################# Public Method ###############################
void DigitalFilter::setFilterParameters(double centerFreq, double cutoff, double samplingFreq)
{
  _fo = centerFreq;
  _fc = cutoff;
  _fs = samplingFreq;
  return;
}

// ############################# Public Method ###############################
// setFilterOrder -- Sets a new filter order
// Input:       order:          New filter order
//          
// Output:                      None
//
// Notes:
// 1. This is just a dummy for overriding in AbstractFIR
// ############################# Public Method ###############################
void DigitalFilter::setFilterOrder(int order)
{
  RealFilter::setFilterOrder(order);
  return;
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
void DigitalFilter::setPassBandGain(LOGICAL usePolesZeros)
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
      arg               = Complex(1., 0.);              // omega = 0
      break;
    case HIGH_PASS:
      wc                = TWOPI*_fc/_fs;
      arg               = Complex(0., (wc + PI)/2.);    // omega = (wc+PI)/2.
      arg               = exp(arg);
      break;
    case BAND_PASS:
      wo                = TWOPI*_fo/_fs;
      arg               = Complex(0., wo);              // omega = wo;
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
void DigitalFilter::zeroOutTaps()
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
    delete [] _doubleShift;

    _shiftSize      = max_shift;
    _doubleShift    = new double[_shiftSize+1];
    _complexShift   = new Complex[_shiftSize+1];
  }

  for(i=0; i<_shiftSize; i++)
  {
    _doubleShift[i]     = 0;
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
void DigitalFilter::setFilterACoeffs(const double *aCoeff)
{
  int   max_shift;
//
// Perform super's function
//
  RealFilter::setFilterACoeffs(aCoeff);
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
void DigitalFilter::setFilterBCoeffs(const double *bCoeff)
{
  int   max_shift;
//
// Perform super's function
//
  RealFilter::setFilterBCoeffs(bCoeff);
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
void DigitalFilter::scaleCoefficientsByPassBandGain()
{
//
// Calculate the pass band gain based on filter pass type
//
  setPassBandGain();
  if(_passBandGain != 0.)
  {
    *_bPolyObject       /= _passBandGain;
    _passBandGain       = 1.;
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
float *DigitalFilter::filterResponseAtFrequencies(float *frequencies, int numberPts, LOGICAL calc, LOGICAL normalize)
{
  int     n, m;
  Complex complex_response, carg, exponent;
  Complex alias_response;
//
//  Error check:
//
  if( (numberPts<1) || (frequencies==NULL) )
    return NULL;
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
    m                   = MAX(0, n);
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
      if(_filterStructureType == POLE_ZERO)
        complex_response        = poleZeroResponseAt(carg, _filterPoles, _filterZeros);
      else 
        complex_response        = transferResponseAt(carg);
      _filterResponse[m]        = outputResponseFor(complex_response, _fs*frequencies[m]);
    }
  }
//
// Normalize if requested:
//
  if(normalize)
  {
   switch(_responseType)
   {
    case MAGNITUDE:
      _passBandGain     = 0.;
      for(n=0; n<numberPts; n++)
      {
        if(_filterResponse[n] > _passBandGain)
          _passBandGain = _filterResponse[n];
      }
      if(_passBandGain > 0.)
        for(n=0; n<numberPts; n++)
          _filterResponse[n]    /= _passBandGain;
      break;
    case DB_MAGNITUDE:
      _passBandGain             = -1000.;
      for(n=0; n<numberPts; n++)
      {
        if(_filterResponse[n] > _passBandGain)
          _passBandGain         = _filterResponse[n];
      }
      for(n=0; n<numberPts; n++)
        _filterResponse[n]      -= _passBandGain;
      break;
    default:
      break;
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
// 1. This includes _passBandGain.
// ############################# Public Method ###############################
float DigitalFilter::filterResponseAtFrequency(float frequency, LOGICAL calc)
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
    exponent                    = Complex(0., frequency);
    carg                        = exp(exponent);
    if(_filterStructureType == POLE_ZERO)
      complex_response          = poleZeroResponseAt(carg, _filterPoles, _filterZeros);
    else 
      complex_response          = transferResponseAt(carg);
    
    filter_response             = outputResponseFor(complex_response, _fs*frequency);
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
Complex DigitalFilter::complexResponseAtFrequency(float frequency, LOGICAL calc)
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
    if(_filterStructureType == POLE_ZERO)
      complex_response          = poleZeroResponseAt(carg, _filterPoles, _filterZeros);
    else 
      complex_response          = transferResponseAt(carg);
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
void DigitalFilter::findTransferFunction()
{
  zeroOutTaps();
  return;
}

// ############################# Public Method ###############################
// findNBW -- Find and return the noise bandwidth in Hertz
//
// Input:               None
//          
// Output:              noise bandwidth in Hertz
//
// Notes:
//  1. The transfer function must first be calculated.
//  2. For bandpass and bandstop filters, the noise bandwidth is calculated
//     over [0, fs/2], otherwise over [-fs/2, fs/2], therefore the 2 sided
//     noise bandwidth of low pass filters is returned.
// ############################# Public Method ###############################
float DigitalFilter::findNBW()
{
  float         noise_bw, mag_squared_max;
  //TODO: Note the findPower() in IIRFilter.cpp doesn't seem to work
  noise_bw              = DigitalFilter::findPower();
  mag_squared_max       = findMagSquaredMax();
  if(mag_squared_max > 0.)
    noise_bw    *= _fs/mag_squared_max;
  else
    noise_bw    = 0.;
  return noise_bw;
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
//  2. The noise power is given by, P = (1/2PI)*Integral(-PI,PI, [H(jw)]**2).
// ############################# Public Method ###############################
double DigitalFilter::findPower()
{
  int           i, old_resp_type;
  double        f_power, magnitude_sq;
  float         freq_min, freq_max, delta_f;
  float         frequency;
//
//  Set response type to magnitude:
//
  old_resp_type = _responseType;
  _responseType = MAGNITUDE;
//
// First, find the frequency range to calculate NBW:
//
  if( (_filterPassType == BAND_PASS) || (_filterPassType == BAND_STOP) )
    freq_min    = 0.;
  else
    freq_min    = -PI;
  freq_max      = PI;
  delta_f       = (freq_max - freq_min)/(NUMBER_NBW_POINTS-1);
//
// Integrate the magnitude function, and find max value:
//
  f_power               = 0.;
  frequency             = freq_min;
  for(i=0; i<NUMBER_NBW_POINTS-1; i++)          // don't integrate last point
  {
    magnitude_sq        = filterResponseAtFrequency(frequency, NO);
    magnitude_sq        *= magnitude_sq;
    f_power             += magnitude_sq;
    frequency           += delta_f;
  }
//
// filterResponseAtFrequency is scaled by 1/_passBandGain:
//
  f_power               *= _passBandGain*_passBandGain*delta_f/TWOPI;
//
// Restore response type:
//
  _responseType         = old_resp_type;
  return f_power;
}

// ############################# Public Method ###############################
// findMagSquaredMax -- Find and return the maximum magnitude squared value
//                      of the filter.
//
// Input:               None
//          
// Output:              maximum magnitude squared value
//
// Notes:
//  1. The transfer function must first be calculated.
// ############################# Public Method ###############################
double DigitalFilter::findMagSquaredMax()
{
  int           i, old_resp_type;
  double        mag_squared_max, magnitude_sq;
  float         freq_min, freq_max, delta_f;
  float         frequency;
//
//  Set response type to magnitude:
//
  old_resp_type = _responseType;
  _responseType = MAGNITUDE;
//
// First, find the frequency range to calculate mag squared:
//
  if( (_filterPassType == BAND_PASS) || (_filterPassType == BAND_STOP) )
    freq_min    = 0.;
  else
    freq_min    = -PI;
  freq_max      = PI;
  delta_f       = (freq_max - freq_min)/(NUMBER_NBW_POINTS-1);
//
// Find max value:
//
  mag_squared_max       = 0.;
  frequency             = freq_min;
  for(i=0; i<NUMBER_NBW_POINTS-1; i++)          // don't integrate last point
  {
    magnitude_sq        = filterResponseAtFrequency(frequency, NO);
    magnitude_sq        *= magnitude_sq;
    if(magnitude_sq > mag_squared_max)
      mag_squared_max   = magnitude_sq;
    frequency           += delta_f;
  }
//
// filterResponseAtFrequency is scaled by 1/_passBandGain:
//
 mag_squared_max        *= _passBandGain*_passBandGain;
//
// Restore response type:
//
  _responseType         = old_resp_type;
  return mag_squared_max;
}
