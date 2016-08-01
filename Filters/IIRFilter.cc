/************************************************************************
 *                                                                      *
 * This subclass of DigitalFilter adds functionality for calculating    *
 * the transfer function and poles and zeros of certain types of IIR    *
 * filters.  The actual filtering functions are defined in the super-   *
 * class.                                                               *
 *                                                                      *
 * File:IIRFilter.h                                                     *
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
 * NOTE: The order of the filter is the # taps - 1.                     *
 *                                                                      *
 ************************************************************************/

#include "IIRFilter.h"                                  // Object prototypes
#include "AnalogFilter.h"

#if defined(WIN32)
#include <Polynomials/DoublePoly.h>
#include <Polynomials/ComplexPoly.h>
#include <GNU/Complex.h>
#include <C_Libraries/constants.h>
#define ADD_RESIDUE     1
#else
#include "DoublePoly.h"
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
// findDigitalPolesZerosFromAnalog --Find the filter's poles and zeros from
//                                   analog prototype
//
// Input:                   None
//          
// Output:                  None
//
// Notes:
//  1. The instance variables, _filterPoles, _filterZeros are modified.
// ############################# Private Method ###############################
void IIRFilter::findDigitalPolesZerosFromAnalog()
{
  int           i;
  double        wc, wo;
  double        wc_warp, wo_warp, w1_warp, w2_warp;
  Complex       carg, s_plane_w;
  const         Complex *analog_poles;

/***************************
 * Conversions:            *
 ***************************/
  wc            = TWOPI*_fc;
  wo            = TWOPI*_fo;
  _fs           = MAX(1.e-5, _fs);
  wo_warp       = warpFrequency(wo, _fs);               // pre-warp center freq
  wc_warp       = warpFrequency(wc, _fs);               // pre-warp cutoff freq
  w1_warp       = (wo-wc);                              // pre-warp bandwidth
  if(w1_warp>0.)
    w1_warp     = warpFrequency(w1_warp, _fs);
  w2_warp       = (wo+wc);                          // pre-warp bandwidth
  if(w2_warp>0.)
    w2_warp     = warpFrequency(w2_warp, _fs);

/***************************
 * Find Analog poles and   *
 * zeros, then convert to  *
 * z domain using bilinear *
 * tranformation.          *
 ***************************/
  if(_analogFilter==NULL)
    _analogFilter       = new AnalogFilter(_filterPassType, wo_warp/TWOPI, wc_warp/TWOPI, _filterOrder);
  _analogFilter->setPassType(_filterPassType);
  _analogFilter->setFilterStructureType(_filterStructureType);
  _analogFilter->setAnalogType(_analogType);
  _analogFilter->setFilterOrder(_filterOrder);
  _analogFilter->setFilterFrequencies(wo_warp/TWOPI, wc_warp/TWOPI);
  _analogFilter->setPassbandRipple(_passbandRipple);

  analog_poles          = _analogFilter->filterPoles();
  _numberPoles          = _analogFilter->numberOfPoles();
  for(i=0; i<_numberPoles; i++)
    _filterPoles[i]    = sToZ(analog_poles[i]);

/*********************************
 * Now find zeros:               *
 *    LPF = (z+1)**m             *
 *    HPF = (z-1)**m             *
 *    BPF = [(z-1)(z+1)]**m      *
 *    BSF = [(s-jwo)(s+jwo)]**m  *
 *********************************/
  _numberZeros = _numberPoles;
  switch(_filterPassType)
  {
    case LOW_PASS:
    default:
      for(i=0; i<_numberPoles; i++)
        _filterZeros[i] = Complex(-1., 0.);
      break;
    case HIGH_PASS:
      for(i=0; i<_numberPoles; i++)
        _filterZeros[i] = Complex(1., 0.);
      break;
    case BAND_PASS:
      for(i=0; i<_numberPoles; i++)
      {
        if(i%2)
          _filterZeros[i] = Complex(1., 0.);
        else
          _filterZeros[i] = Complex(-1., 0.);
      }
      break;
    case BAND_STOP:
      s_plane_w = Complex(0., wo_warp);
      carg      = sToZ( s_plane_w );
      for(i=0; i<_numberPoles; i++)
      {
        if(i<_numberPoles/2)
          _filterZeros[i] = carg;
        else
          _filterZeros[i] = conj(carg);
      }
      break;
  }
//
// Set pass band gain for unity:
//
  setPassBandGain(YES);

  return;                     
}

// ############################# Private Method ###############################
// transferForIIRPowerLawFilter -- This routine finds the coefficients of the IIR
//                              filter for implementing a 1/(f^alpha) power law
//                              response.
//
// Input:                   None
//          
// Output:                  Return YES if OK, return NO if polynomials don't exist
//
// Notes:
//  1. The instanceVariables, _aPolyObject and _bPolyObject are modified.
//  2. The number of taps is equal to order + 1.
//  3. The instance variables, _powerLawAlpha  affect
//     the output, in addition to the normal filter parameters.
//  4. Based on N.J. Kasdin, "Discrete simulation of ...", Proc. IEEE, vol. 83, no.5
//  5. A constant, BETA, is used to move the poles inside the unit circle.
// ############################# Private Method ###############################
LOGICAL IIRFilter::transferForIIRPowerLawFilter()
{
  int       n, taps, old_type;
  double    *a, *b;
  
  if(_aPolyObject == NULL)
     return NO;
   
  taps          = _filterOrder+1;

//
// Set up the filter coefficients according to equation (116) in the reference
// Note the typo in the equation.
//
  if(taps)
  {
    a           = new double[taps];
    b           = new double[1];
    a[0]        = 1.;
    for(n=1; n<taps; n++)
      a[n]      = POWER_LAW_BETA*(n - 1. - _powerLawAlpha/2.)*a[n-1]/n;
//
// Store coefficients in polynomials:
//
    b[0] = 1.;
    _aPolyObject->assign(_filterOrder, a);
    _bPolyObject->assign(0, b);
    delete [] a;
    delete [] b;
//
// Set unity pass band gain:
//
    old_type            = _filterPassType;
    _filterPassType     = LOW_PASS;
    setPassBandGain();
    _filterPassType     = old_type;
  }

  return YES;
}

#define MAX_CIC_STAGES      10
// ############################# Private Method ###############################
// transferForCICFilter -- This routine finds the coefficients of the low pass
//                         CIC filter as described in Hogenauer, "An economical
//                         class ...", IEEE ASSP, Vol. ASSP-29, no. 2, April 1981.
//
// Input:                   None
//          
// Output:                  Return YES if OK, return NO if polynomials don't exist
//
// Notes:
//  1. The instanceVariables, _aPolyObject and _bPolyObject are modified.
//  2. The _filterOrder is used to set the number of stages, _decimateFactor is also
//     used.
// ############################# Private Method ###############################
LOGICAL IIRFilter::transferForCICFilter()
{
  int       i, n;
  double    *c, *d;
  DoublePoly c_poly, d_poly;

  if(_aPolyObject == NULL)
     return NO;
//
// Allocate multiplication polynomials:
//
  c             = new double[_decimateFactor + 1];
  d             = new double[2];
  for(i=1; i<_decimateFactor; i++)
    c[i]        = 0.;
  d[0]          = c[0]                  = 1.;
  d[1]          = c[_decimateFactor]    = -1.;
//
// Initialize the coefficients in the polynomials:
//
  _aPolyObject->assign(1, d);
  _bPolyObject->assign(_decimateFactor, c);
  c_poly.assign(_decimateFactor, c);
  d_poly.assign(1, d);

//
// Now, find the polynomials by multiplying
//
  n     = MIN(_filterOrder, MAX_CIC_STAGES);
  for(i=1; i<n; i++)
  {
    *_aPolyObject       *= d_poly;
    *_bPolyObject       *= c_poly;
  }
//
// Divide a out of b, there should be no remainder
//
  *_bPolyObject         /= *_aPolyObject;
  _aPolyObject->assign(0, d);

//
// Delete allocated arrays:
//
  delete [] c;
  delete [] d;
//
// Set unity pass band gain:
//
  setPassBandGain();
  return YES;
}

// ############################# Private Method ###############################
// initInstanceVariables -- This routine initializes the object's instance variables
//
// Input:                       None
//          
// Output:                      None
//
// Notes:
// ############################# Private Method ###############################
void IIRFilter::initInstanceVariables()
{
  _analogFilter         = NULL;
//
// Add these:
//
  setIIRFilterType(ANALOG_PROTOTYPE);
  setPowerLawAlpha(DEFAULT_POWER_ALPHA);
  setDecimateFactor(DEFAULT_CIC_DECIMATION);
  setPassbandRipple(DEFAULT_RIPPLE);

  
  return;
}

// ############################# Private Method ###############################
// residueAtPoint -- This routine finds the residue at a complex point
//
// Input:       numerator       Numerator polynomial b(z)%[a(z)*a(z^-1)]
//              denominator     Denominator polynomial a(z)
//              denomRev        Denominator reverse a(z^-1)
//              point           Point to find residue
//          
// Output:                      Residue
//
// Notes:
// ############################# Private Method ###############################
Complex IIRFilter::residueAtPoint(DoublePoly &numerator, DoublePoly &denominator, DoublePoly &denomRev, Complex &point)
{
  Complex       residue_coeff[2], denominator_value, reciprocal_point;
  Complex       numerator_value = 0.;
#if(ADD_RESIDUE)
  ComplexPoly   complex_num, complex_den, residue_poly, complex_den_rev;

//
// Copy polys to complex polys
//
  complex_num           = numerator;
  complex_den           = denominator;
  complex_den_rev       = denomRev;
//
// Divide the denominator by the residue.
//
  residue_coeff[1]      = Complex(1., 0.);
  residue_coeff[0]      = -1.*point;
  residue_poly.assign(1, residue_coeff);
  complex_den_rev       /= residue_poly;
  if(norm(point) != 0.)
    reciprocal_point    = Complex(1.,0.)/point;
  else
    reciprocal_point    = point;
//
// Now, evaluate numerator and denominator at point/pole
//
  numerator_value       = complex_num.evaluateAtComplexPoint(reciprocal_point);
  denominator_value     = complex_den.evaluateAtComplexPoint(reciprocal_point);
  denominator_value     *= complex_den_rev.evaluateAtComplexPoint(reciprocal_point);
//
// Divide numerator by denominator, and return result:
//
  if(norm(denominator_value) != 0.)
    numerator_value     /= denominator_value;
#endif
  
  return numerator_value;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the IIRFilter class.
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
IIRFilter::IIRFilter(int pass, double centerFreq, double cutoffFreq,
                                  double samplingFreq, int order)
              :DigitalFilter(pass, centerFreq, cutoffFreq, samplingFreq, order)
{
  initInstanceVariables();
  
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the IIRFilter class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
IIRFilter::~IIRFilter()
{
  delete _analogFilter;
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
const Complex *IIRFilter::filterPoles()
{
  int       i;
  Complex   *roots;

  findTransferFunction();                               // calculate the poles and zeros or xfer fn
  switch(_iirFilterType)
  {
    case ANALOG_PROTOTYPE:
      findDigitalPolesZerosFromAnalog();
      break;
    default:
      if(_aPolyObject != NULL)
      {
        delete []       _filterPoles;                   // Find poles and zeros from xfer function
        _numberPoles    = _aPolyObject->getOrder();
        _filterPoles    = new Complex[_numberPoles];
        roots           = _aPolyObject->roots();
        for(i=0; i<_numberPoles; i++)
          _filterPoles[i]  = roots[i];
       }
       break;
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
const Complex *IIRFilter::filterZeros()
{
  int       i;
  Complex   *roots;

  findTransferFunction();                               // calculate the poles and zeros or xfer fn
  switch(_iirFilterType)
  {
    case ANALOG_PROTOTYPE:
      findDigitalPolesZerosFromAnalog();
      break;
    default:
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
       break;
  }
  return _filterZeros;
}

// ############################# Public Method ###############################
// findTransferFunction -- Find the filter's transfer function
//
// Input:                   None
//          
// Output:                  None
//
// Notes:
//  1. The instanceVariables, _aPolyObject and _bPolyObject are modified.
// ############################# Public Method ###############################
void IIRFilter::findTransferFunction()
{
  switch(_iirFilterType)
  {
    case IIR_POWER_LAW:
      transferForIIRPowerLawFilter();
      break;
    case CIC:
      transferForCICFilter();
      break;
    case ANALOG_PROTOTYPE:
    default:
      findDigitalPolesZerosFromAnalog();
      transferFromPoles(_filterPoles, _filterZeros, _numberPoles, _numberZeros, _analogFilter->realAxisPoles());
      break;
    case USER_DEFINED_FILTER:
      break;                                                            // Do nothing
  }
  zeroOutTaps();
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
double IIRFilter::warpFrequency(double omega, double fs)
{
  double        omega_warp      = 0.;
  if(fs != 0.)
  {
    omega_warp  = omega/fs;
  }
  omega_warp    = 2.*fs*tan(omega_warp/2.);
  return omega_warp;
}

// ############################# Public Method ###############################
// sToZ -- Convert s plane pole to z plane pole using Bilinear
//         transformation: s -> 2(z-1)/(z+1)/T
//
// Input:       sPoint:     Complex point in the s plane
//
// Output:                  Complex point in the z plane
//
// Notes:
// ############################# Public Method ###############################
Complex IIRFilter::sToZ(const Complex &sPoint)
{
 double  t;
 Complex p_t_2, cnum, cden;

 if(_fs>0.)
   t   = 1/_fs;
 else
   t   = 1.;
 p_t_2 = sPoint * t/2.;

 cnum  = 1. + p_t_2;
 cden  = 1. - p_t_2;
 if(abs(cden) > 0.)
   p_t_2   = cnum/cden;

 return p_t_2;
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
void IIRFilter::convertCoefficientsToZ()
{
  int           i, j, power;
  int           a_order, b_order;
  const         double *a_coeffs, *b_coeffs;
  double        ts, a0;
  double        tran_coeff[2]   = {2., -2.};            // 2z-2
  double        den_coeff[2]    = {1., 1.};             // z+1
  DoublePoly    sum_poly, trans_poly, temp_poly;
  DoublePoly    den_poly, mult_poly;
//
// set the transform polynomial:
//
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
  mult_poly     = 1.;
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
  if(a0 != 0.)
  {
    *_aPolyObject       /= a0;
    *_bPolyObject       /= a0;
  }
  
  return;
}

// ############################# Public Method ###############################
// convertCoefficientsToZUsingBackwardDifference -- Convert coefficients previously
//                           loaded using setFilterACoeffs() or setFilterBCoeffs()
//                           to z domain using s->(z-1)/zT
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
void IIRFilter::convertCoefficientsToZUsingBackwardDifference()
{
  int           i, j, power;
  int           a_order, b_order;
  const         double *a_coeffs, *b_coeffs;
  double        ts;
  double        poly_coeff[2]   = {1., -1.};            // z-1
  DoublePoly    sum_poly, trans_poly, temp_poly;
//
// set the transform polynomial:
//
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
// findPower -- Find and return the normalized power out of the filter.
//
// Input:               None
//          
// Output:              noise power in Watts
//
// Notes:
//  1. The transfer function must first be calculated.
//  2. The noise power is given by, P = (1/2PI)*Integral(-PI,PI, [H(jw)]**2).
//     From Oppenheim and Schafer, Eqn (A.66), Power = Sum(k=1,N; Ak), where
//     from Eqn (A.65) Ak = H(z)H(z^-1)(z-zp) - evaluated at z = zp.  That
//     is, the residue of H(z)H(z^-1).
// ############################# Public Method ###############################
double IIRFilter::findPower()
{
  int           i, n_roots;
  double        f_power;
  const double  *poly_coeffs;
  DoublePoly    a, b, a_rev, b_rev, b_div_a;
  Complex       *a_roots, residue_sum;
  Complex       residue;
//
//  Copy the polynomials:
//
  a             = *_aPolyObject;
  b             = *_bPolyObject;
  a_rev         = a;
  b_rev         = b;
  a_rev.reverse();
  b_rev.reverse();
//
// Get the roots of the denominator, a, and find the residues at the roots:
//
  n_roots       = a.getOrder();
  a_roots       = a.roots();
//
// Find H(z)*H(z^-1)
//
  a             *= a_rev;
  b             *= b_rev;
//
// If order of b >= a, then divide it out.
//
  if(b.getOrder() >= a.getOrder())
  {
    b_div_a     = b/a;
    b           %= a;
    poly_coeffs = b_div_a.getCoefficients();
    residue_sum = Complex(poly_coeffs[0], 0.);
  }
  else
  {
    residue_sum = Complex(0., 0.);
  }
//
// Sum the residues, first restore a to ease accuracy
// requirements for finding residue.
//
  a             = a_rev;
  a.reverse();
  for(i=0; i<n_roots; i++)
  {
    residue     = residueAtPoint(b, a, a_rev, a_roots[i]);
    residue_sum += residueAtPoint(b, a, a_rev, a_roots[i]);
  }
//
// Find the power, it should be real:
//
  f_power       = real(residue_sum);
//
// Restore response type:
//
  return f_power;
}
