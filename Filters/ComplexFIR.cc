/************************************************************************
 *                                                                      *
 * This subclass of ComplexDigitalFilter adds functionality for         *
 * filtering using an FIR structure.  For calculating filter            *
 * coefficients, see the subclasses.                                    *
 *                                                                      *
 * File:ComplexFIR.h                                                    *
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
 * Revision history:                                                    *
 *  1. 04/20/05 - Modified from AbstractFIR.                            *
 *  2. 04/25/06 - Added findPower()                                     *
 ************************************************************************/
#include "ComplexDigitalFilter.h"
#include "ComplexFIR.h"                                 // Object prototypes
#if defined(WIN32)
#include <Polynomials/ComplexPoly.h>
#include <GNU/Complex.h>
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


// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the ComplexFIR class.
//
// Input:       bCoeffs:                Filter coefficients
//              numberCoeffs:           number of coefficients (order + 1)
//
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
ComplexFIR::ComplexFIR(Complex *bCoeffs, int numberCoeffs)
           :ComplexDigitalFilter(NULL, bCoeffs, numberCoeffs)
{
  _numberTaps   = numberCoeffs;
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the ComplexFIR class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
ComplexFIR::~ComplexFIR()
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
void ComplexFIR::setNumberTaps(int taps)
{
  _numberTaps   = MAX(taps, MIN_TAPS);
  _numberTaps   = MIN(_numberTaps, MAX_TAPS);
//
// Call super's method
//
  ComplexDigitalFilter::setFilterOrder(_numberTaps-1);
  return;
}

// ############################# Public Method ###############################
// setFilterOrder -- Sets the order of the filter (numberTaps - 1).
// Input:       taps            New number of filter taps
//          
// Output:                      none
// Notes:
// ############################# Public Method ###############################
void ComplexFIR::setFilterOrder(int order)
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
const Complex *ComplexFIR::filterPoles()
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
const Complex *ComplexFIR::filterZeros()
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
void ComplexFIR::setPassBandGain(LOGICAL usePolesZeros)
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
void ComplexFIR::zeroOutTaps()
{
//
// Perform super's function
//
  ComplexDigitalFilter::zeroOutTaps();
  _shiftPtr     = MAX((_shiftSize - 1), 0);     // The pointer decrements, so set to end
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
Complex ComplexFIR::filterFloatData(float input)
{
  int           i, shift_ptr;
  double        gain;
  Complex       yout;
  const         Complex *b;

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
// filterDoubleData -- Filter an input data sample, and return the sampled value.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   input:          An input data sample
//          
// Output:                  The filtered data sample
// ############################# Public Method ###############################
Complex ComplexFIR::filterDoubleData(double input)
{
  int           i, shift_ptr;
  double        gain;
  Complex       yout;
  const         Complex *b;

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
// filterComplexArray -- Filter the input array and overwrite it with output.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   x:              An array of floating point data, modified on output
//          numberPts:      Number of points in the input array.
//          
// Output:                  YES if went OK, No if _bPolyObject doesn't exist.
// ############################# Public Method ###############################
LOGICAL ComplexFIR::filterComplexArray(Complex *x, int numberPts)
{
  int           j, i, order, shift_ptr;
  double        gain;
  const         Complex *b, *b_ptr;
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
Complex  ComplexFIR::filterComplexData(const Complex &input)
{
  int           i, shift_ptr;
  double        gain;
  const         Complex *b;
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
double ComplexFIR::findPower()
{
  int           i;
  double        f_power;
  const Complex *b_coeffs;
//
// Sum the squares of the coefficients:
//
  b_coeffs      = bCoeffs();
  f_power       = 0.;
  for(i=0; i<_numberTaps; i++)
    f_power     += norm(b_coeffs[i]);
  return f_power;
}

