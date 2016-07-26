/************************************************************************
 *                                                                      *
 * This subclass of object implements an abstract filter class.  For    *
 * actual implementations, see the subclasses, DigitalFilter and        *
 * AbstractFilter.                                                      *
 *                                                                      *
 * File:AbstractFilter.cc                                               *
 *                                                                      *
 * The filter is stored in two forms:                                   *
 *                                                                      *
 *          b[n] + b[n-1]s + ... + b[0]s**(n)                           *
 *    H(s) = ------------------------------------                       *
 *          1    + a[n-1]s + ... + a[0]s**(n)                           *
 *                                                                      *
 *           (s-zero[0]) * (s-zero[1]) ... (s-zero[n_zero])             *
 *    H(s) = ----------------------------------------------             *
 *           (s-pole[0]) * (s-pole[1]) ... (s-pole[n_pole])             *
 *                                                                      *
 * or equivalently, for digital filters:                                *
 *                                                                      *
 * NOTE: The order of the filter is n, but the number of coefficients   *
 *       is equal to n+1.                                               *
 * Where the a and b polynomials are stored in the polynomial objects.  *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 01/05/04  - Abstracted from AbstractFilter                       *
 ************************************************************************/

#include "AbstractFilter.h"                                     // Object prototypes

#if defined(WIN32)
#include <GNU/Complex.h>
#include <C_Libraries/constants.h>
#else
#include "Complex.h"
#include "constants.h"
#endif

#include <math.h>
#include <stdio.h>


#define ABS(a)          ((a) >= 0 ? (a) : (-a))
#define SIGN(a,b)       ( ((b)>=0 ) ? ABS((a)) : -ABS((a)) )
#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define SGN01(a)        ( ((a)>=0.) ? 1 :  0 )
#define MAX_ORDER   10                          /* max prototype order */


// ############################# Private Method ###############################
// initInstanceVariables -- This routine initializes the object's instance variables
//
// Input:                       None
//          
// Output:                      None
//
// Notes:
// ############################# Private Method ###############################
void AbstractFilter::initInstanceVariables()
{

// init instance variables:
  setPassType(LOW_PASS);
  setFilterStructureType(TRANSFER_FUNCTION);
  setAnalogType(BUTTERWORTH);
  setResponseType(MAGNITUDE);
  setPhaseInDegrees(NO);


  _filterPoles          = _filterZeros  = NULL;
  _filterResponse       = NULL;

  _passBandGain         = 1.;
  _oldNumberFrequencies = 0;
  _thetaOld             = 0.;
  _subAngle             = 0.;
  
  return;
}

// ############################# Private Method ###############################
// initPoleZeroArrays -- Allocates the arrays for the Complex poles and zeros.
//
// Input:                       None
//          
// Output:                      None
//
// Notes:
// ############################# Private Method ###############################
void AbstractFilter::initPoleZeroArrays()
{
  int   n;

  delete [] _filterPoles;
  delete [] _filterZeros;

  n             = MAX(_numberPoles, _numberZeros);
  _filterPoles  = new Complex[2*n];
  _filterZeros  = new Complex[2*n];
  return;
}

// ############################# Private Method ###############################
// poleZeroResponseAt -- Finds the response of the filter at a complex point using
//                       the pole zero representation of the filter.
//
// Input:       point:          Complex point in the s plane
//              poles:          Complex array of pole locations
//              zeros:          Complex array of zero locations
//          
// Output:                      None
//
// Notes:
//  1. Since this routine can be also be used for a digital filter, we pass
//     in the poles and zeros arrays.
// ############################# Private Method ###############################
Complex AbstractFilter::poleZeroResponseAt(Complex &point, Complex *poles, Complex *zeros)
{
  int     i;
  Complex num, den, complex_response;

/****************************
 * Find numerator and de-   *
 * nominator for pole zero: *
 ****************************/
  num   = Complex(1., 0.);
  for(i=0; i<_numberZeros; i++)
    num *= point - zeros[i];

  den   = Complex(1., 0.);
  for(i=0; i<_numberPoles; i++)
    den *= point - poles[i];
  if(norm(den) > 0.)
    complex_response = num/den;
  else
    complex_response = Complex(0.,0.);
  
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
// ############################# Private Method ###############################
float AbstractFilter::outputResponseFor(Complex cResponse , float omega)
{
  float output_response;
  double temp, theta_new;

  if(_passBandGain == 0.)
    _passBandGain       = 1.;                   // Check for error condition
  switch(_responseType)
  {
    case MAGNITUDE:
    default:
      output_response   = abs(cResponse);
      output_response   /= _passBandGain;
      break;
    case DB_MAGNITUDE:
      temp              = norm(cResponse);
      temp              /= _passBandGain*_passBandGain;
      if(temp > 0.)
        output_response = 10.*log10(temp);
      else
        output_response = 0.;
      break;
    case PHASE:
    case PHASE_NORMALIZED:                      // This is handled in FIR filter
      output_response   = arg(cResponse);
      if(ABS((_thetaOld-output_response)) > PI)         // Check for 0 crossing
      {
        if(_thetaOld > 0.)
          _subAngle     -= TWOPI;
        else
          _subAngle     += TWOPI;
      }
      _thetaOld         = output_response;
      output_response   -= _subAngle;
      if(_phaseInDegrees)
        output_response *= DEG_RAD;
      break;
    case PHASE_DELAY:
      output_response   = atan2(imag(cResponse), real(cResponse));
      if((_thetaOld<=0.) && (output_response>0.))        /* Pi crossing*/
        _subAngle       += TWOPI;
      _thetaOld         = output_response;
      output_response   = output_response - _subAngle;
      if(omega > 0.)
        output_response /= -omega;
      break;
    case GROUP_DELAY:
      theta_new         = atan2(imag(cResponse), real(cResponse));
      if(ABS((_thetaOld-theta_new)) < PI)               // Check for 0 crossing
        output_response = (_thetaOld - theta_new)/_deltaOmega;
      else
      {
        if(_thetaOld > 0.)
          output_response       = (_thetaOld - theta_new - TWOPI)/_deltaOmega;
        else
          output_response       = (_thetaOld - theta_new + TWOPI)/_deltaOmega;
      }
      _thetaOld         = theta_new;
      break;
  }
  return output_response;
}
  
// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the AbstractFilter class.
//
// Input:           type:       Band pass, low pass, etc.
//                  centerFreq: Filter center frequency in Hertz
//                  cutoffFreq: Filter cutoff frequency in Hertz
//                  order:      The order of the filter
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
AbstractFilter::AbstractFilter(int type, double centerFreq, double cutoffFreq, int order)
{
  initInstanceVariables();
  setFilterOrder(order);
  setPassType(type);
  setFilterFrequencies(centerFreq, cutoffFreq);

  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the AbstractFilter class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
AbstractFilter::~AbstractFilter()
{
  delete [] _filterPoles;
  delete [] _filterZeros;
  delete [] _filterResponse;
  return;
}

// ############################# Public Method ###############################
// setPassType -- Sets a new filter pass type
// Input:   type:               LOW_PASS, BAND_PASS, etc..
//          
// Output:                      None
// ############################# Public Method ###############################
void AbstractFilter::setPassType(int type)
{
  _filterPassType       = type;
  return;
}

// ############################# Public Method ###############################
// setFilterStructureType -- Sets a new filter structure type
// Input:   type:               TRANSFER_FUNCTION, POLE_ZERO.
//          
// Output:                      None
// ############################# Public Method ###############################
void AbstractFilter::setFilterStructureType(int type)
{
  _filterStructureType  = type;
  return;
}

// ############################# Public Method ###############################
// setAnalogType -- Sets the type of analog (prototype) filter.
// Input:   type:               CHEBYSHEV, BUTTERWORTH, etc.
//          
// Output:                      None
// ############################# Public Method ###############################
void AbstractFilter::setAnalogType(int type)
{
  _analogType   = type;
  return;
}

// ############################# Public Method ###############################
// setResponseType -- Sets a filter response type
// Input:   type:               MAGNITUDE, DB_MAGNITUDE, etc.
//          
// Output:                      None
// ############################# Public Method ###############################
void AbstractFilter::setResponseType(int type)
{
  _responseType = type;
  return;
}


// ############################# Public Method ###############################
// setPhaseInDegrees -- Sets whether or not to calculate phase in degrees
// Input:   flag:               Yes = calculate in degrees.
//          
// Output:                      None
// ############################# Public Method ###############################
void AbstractFilter:: setPhaseInDegrees (LOGICAL flag)
{
  _phaseInDegrees       = flag;
  return;
}


// ############################# Public Method ###############################
// setAbstract_filterOrder -- Sets the order of the filter
// Input:   order:              New filter order
//          
// Output:                      None
// ############################# Public Method ###############################
void AbstractFilter::setFilterOrder(int order)
{
  _filterOrder       = MAX(0, order);
  _numberPoles       = _numberZeros = _filterOrder;
  initPoleZeroArrays();
  return;
}

// ############################# Public Method ###############################
// setFilterFrequencies -- Sets filter center and cutoff frequencies
// Input:   centerFreq:         filter's center frequency in Hertz
//          cutoff:             filter's cutoff frequency in Hertz
//          
// Output:                      None
// ############################# Public Method ###############################
void AbstractFilter::setFilterFrequencies(double centerFreq, double cutoff)
{
  _fo = centerFreq;
  _fc = cutoff;
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
// NOTES:
// 1. This needs to be overridden in subclasses.
// ############################# Public Method ###############################
void AbstractFilter::setPassBandGain(LOGICAL usePolesZeros)
{
  return;
}
