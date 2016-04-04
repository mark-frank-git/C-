/************************************************************************
 *                                                                      *
 * This subclass of object implements an analog filter object           *
 *                                                                      *
 * File:AnalogFilter.cc                                                 *
 *                                                                      *
 * The filter is stored in two forms:                                   *
 *                                                                      *
 *          b[0]s**(n) + b[1]s**(n-1) + ... +b[n-1]                     *
 *    H(s) = ------------------------------------                       *
 *          s**(n)    +  a[1]s**(n-1) + ... +a[n-1]                     *
 *                                                                      *
 *           (s-zero[0]) * (s-zero[1]) ... (s-zero[n_zero])             *
 *    H(s) = ----------------------------------------------             *
 *           (s-pole[0]) * (s-pole[1]) ... (s-pole[n_pole])             *
 *                                                                      *
 * Where the a and b polynomials are stored in the polynomial objects.  *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 12/30/94  - Objective-C -> C++                                   *
 *  2. 01/27/95  - x.real() -> real(x), add ASINH_EPS.                  *
 *  3. 02/27/95  - Fixed bug in destructor.                             *
 *  4. 07/12/95  - Fixed new on _filterPoles, _filterZeros.             *
 *  5. 02/20/97  - Reversed order of denominator, numerator polys.      *
 *  6. 02/25/97  - Add getACoeffs(), getBCoeffs().                      *
 *  7. 01/05/98  - make subclass of AbstractFilter.                     *
 *  8. 10/08/99  - make subclass of RealFilter.                         *
 *  9. 11/13/00  - Add _passbandRipple (used to be equal to 3 dB)       *
 ************************************************************************/

#include "AnalogFilter.h"                                       // Object prototypes

#if defined(WIN32)
#include <Polynomials/DoublePoly.h>
#include <GNU/Complex.h>
#include <C_Libraries/constants.h>
#include <Specfuns/specfuns.h>
#else
#include "DoublePoly.h"
#include "Complex.h"
#include "constants.h"
#include "specfuns.h"
#endif

#include <math.h>

#define ABS(a)      ((a) >= 0 ? (a) : (-a))
#define SIGN(a,b)   ( ((b)>=0 ) ? ABS((a)) : -ABS((a)) )
#define MAX(a,b)    ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)    ( ((a)>(b)) ? (b) : (a) )
#define SGN01(a)  ( ((a)>=0.) ? 1 :  0 )

// ############################# Private Method ###############################
// findPolesZeros -- Find poles and zeros based on cutoff freq, filter type, etc.
// Input:                           None
//          
// Output:                          None
//
// Notes:
//  1. The instance variables, _filterPoles and _filterZeros are modified.
//  2. The function, findAnalogPolesZerosFor() does all the work.
// ############################# Private Method ###############################
void AnalogFilter::findPolesZeros()
{
  double  wc, wo;
  wc = TWOPI*_fc;
  wo = TWOPI*_fo;
  
  findAnalogPolesZerosFor(wc, wo, (wc+wc));

  return;
}

// ############################# Private Method ###############################
// findAnalogPolesZerosFor -- Find poles and zeros based on cutoff freq, filter type, etc.
// Input:           wc:         Filter cutoff frequency in rad/s
//                  wo:         Filter center frequency in rad/s
//                  bw:         Filter two sided bandwidth in rad/s
//          
// Output:                      None
//
// Notes:
//  1. The instance variables, _filterPoles and _filterZeros are modified.
// ############################# Private Method ###############################
void AnalogFilter::findAnalogPolesZerosFor(double wc, double wo, double bw)
{
  int i, n;
  Complex carg, response;

  if(bw <= 0.)
    bw                  = 1.;
  _realAxisPoles        = NO;
  _filterOrder          = MIN(MAX_ORDER, _filterOrder);
  n                     = _filterOrder;
  if(n<=0)
  {
    _filterPoles[0] = 0.;
    _filterZeros[0] = 0.;
  }
  else if( (wc==0.)  )              /* degenerate case     */
    for(i=0; i<n; i++)
    {
      _filterPoles[0] = 0.;
      _filterZeros[0] = 0.;
    }
  else
  {
    switch(_analogType)              /* LPF wc = 1 */
    {
      case BUTTERWORTH:
      default:
        butterPrototype();      
        break;
      case CHEBYSHEV:
        chebyPrototype();
        break;
      case BESSEL:
        besselPrototype();
        break;
    }
    switch(_filterPassType)          /* convert to wc filter */
    {
      case LOW_PASS:
      default:
        lowPassToLowPass(wc);
        carg    = 0.;
        break;
      case HIGH_PASS:
        lowPassToHighPass(wc);
        carg    = Complex(0., 4.*wc);
        break;
      case BAND_PASS:
        if(wo == 0.)
          lowPassToLowPass(wc);
        else
          lowPassToBandPass(wo, bw);
          carg  = Complex(0., wo);
          break;
      case BAND_STOP:
        if(wo==0.)
          lowPassToHighPass(wc);
        else
          lowPassToBandStop(wo, bw);
        carg    = Complex(0., 0.);
        break;
    }
    response            = poleZeroResponseAt(carg, _filterPoles, _filterZeros);
    _passBandGain       = abs(response);
    _passBandGain       = MAX(1.e-40, _passBandGain);
  }
  return;
}

// ############################# Private Method ###############################
// rippleToEps -- Converts passband ripple given in dB to epsilon.
//
// Input:       ripple:         Passband ripple given in dB
//          
// Output:                      epsilon
//
// Notes:
// ############################# Private Method ###############################
double AnalogFilter::rippleToEps(double ripple)
{
  double epsilon, arg;

  arg   = pow(10., ripple/10.);
  if(arg > 1.)
    epsilon     = sqrt(arg - 1.);
  else
    epsilon     = 1.;
  return epsilon;
}

// ############################# Private Method ###############################
// butterPrototype -- This routine calculates the poles for a Butterworth filter
//                      with cut-off frequency of wc = 1 rad/s and a pass band
//                      ripple of 3 dB.  See Budak, Passive and Active Network
//                      Analysis Synthesis, p. 506  
//
// Input:                       None
//          
// Output:                      None
//
// Notes:
//  1. The instance variables, _filterPoles and _filterZeros are modified.
// ############################# Private Method ###############################
void AnalogFilter::butterPrototype()
{
  int   i, k, n;
  double arg, eps_to_n;
//
// Convert passband ripple to epsilon:
//
  n                     = _filterOrder;
  eps_to_n              = rippleToEps(_passbandRipple);
  if(n > 0)
    eps_to_n            = 1./pow(eps_to_n, 1./n);
  else
    eps_to_n            = 1.;

/***************************
 * Find poles for wc = 1   *
 ***************************/
  for(k=0; k<n/2; k++)
  {                         
    arg                 = PI*(1.+k+k)/2./n;
    _filterPoles[k]     = Complex(-eps_to_n*sin(arg), eps_to_n*cos(arg));
  }
/*********************************
 * If n is odd, add pole at s = -1*
 *********************************/
  if( n%2 )
  {
    _filterPoles[k]     = Complex(-eps_to_n, 0.);
    k++;
  }
/********************************
 * Copy upper half to lower *
 ********************************/
  i = 0;
  for( ; k<n; k++)
  {
    _filterPoles[k]     = conj(_filterPoles[i]);
    i++;
  }
  return;
}

// ############################# Private Method ###############################
// chebyPrototype -- This routine calculates the poles for a Chebyshev filter
//                      with cut-off frequency of wc = 1 rad/s and a pass band
//                      ripple of 3 dB.  See Budak, Passive and Active Network
//                      Analysis Synthesis, p. 515  
//
// Input:                       None
//          
// Output:                      None
//
// Notes:
//  1. The instance variables, _filterPoles and _filterZeros are modified.
// ############################# Private Method ###############################
void AnalogFilter::chebyPrototype()
{
  int   i, k, n;
  double arg, sinh_value, cosh_value;
  double asinh_eps, eps;
//
// Convert passband ripple to epsilon:
//
  n                     = _filterOrder;
  eps                   = rippleToEps(_passbandRipple);
  asinh_eps             = asinh(1./eps);

/***************************
 * Find poles for wc = 1   *
 ***************************/
  sinh_value    = sinh(asinh_eps/n);
  cosh_value    = cosh(asinh_eps/n);
  for(k=0; k<n/2; k++)
  {                         
    arg         = PI*(1.+k+k)/2./n;
    _filterPoles[k]     = Complex(-sin(arg)*sinh_value, cos(arg)*cosh_value);
  }
/*********************************
 * If n is odd, add real axis pole*
 *********************************/
  if( n%2 )
  {
    _filterPoles[k]     = Complex( -sinh_value, 0.);
    k++;
  }

/********************************
 * Copy upper half to lower *
 ********************************/
  i = 0;
  for( ; k<n; k++)
  {
    _filterPoles[k]     = conj(_filterPoles[i]);
    i++;
  }
  return;
}

// ############################# Private Method ###############################
// besselPrototype -- This routine calculates the poles for a Bessel filter
//                      with cut-off frequency of wc = 1 rad/s See Budak, Passive
//                      and Active Network Analysis Synthesis, p. 525   
//
// Input:                       None
//          
// Output:                      None
//
// Notes:
//  1. The instance variables, _filterPoles and _filterZeros are modified.
// ############################# Private Method ###############################
void AnalogFilter::besselPrototype()
{
  int           i, j, n;
  double        *a, omega_c;
  Complex       *roots, ctemp, response;
  DoublePoly    *polynomial;

/***************************
 * Find poles for wc = 1   *
 ***************************/
  n     = _filterOrder;
  
  if(n>0)
  {
    a   = new double[(n+1)];
    a[0]        = a[1] = 1.;
    for(i=0; i<n; i++)
      a[i+1]    = a[i]*2.*(n-i)/(n+n-i)/(i+1.);
    polynomial  = new DoublePoly(n, a);
    if( (roots = polynomial->roots()) != NULL)
    {
        for(i=0; i<n; i++)
          _filterPoles[i]       = roots[i];
    }
    delete polynomial;
    delete [] a;
  }
//
// We need to normalize frequency
// so that response passes through 3 dB point
// at w == 1
//
  j             = _numberZeros;
  _numberZeros  = 0;
  ctemp         = Complex(0., 0.);
  response      = poleZeroResponseAt(ctemp, _filterPoles, _filterZeros);
  _passBandGain = abs(response);
  _passBandGain = MAX(1.e-40, _passBandGain);
  omega_c       = zeroIn(0.5, 5.);
  omega_c       = MAX(omega_c, 1.e-5);
  omega_c       = 1./omega_c;
  _numberZeros  = j;
  for(i=0; i<_numberPoles; i++)
    _filterPoles[i] *= omega_c;

  return;
}

// ############################# Private Method ###############################
// lowPassToLowPass -- This routine moves the poles of a filter having a cut off
//                     frequency of  1 rad/s to wc rad/s    
//
// Input:               wc:     New cutoff frequency in rad/s
//          
// Output:                      None
//
// Notes:
//  1. The instance variables, _filterPoles are modified.
//  2. Should we be setting _filterZeros to something???????????
// ############################# Private Method ###############################
void AnalogFilter::lowPassToLowPass(double wc)
{
  int n;

  _numberZeros  = 0;
  n             = _numberPoles = _filterOrder;
/***************************
 * Multiply poles by wc    *
 ***************************/
  while(n--)
    _filterPoles[n]     *= wc;

  return;
}

// ############################# Private Method ###############################
// lowPassToHighPass -- This routine moves the poles of a filter having a cut off
//                     frequency of  1 rad/s to wc rad/s, and transforms them to
//                     high pass poles.
//
// Input:               wc:     New cutoff frequency in rad/s
//          
// Output:                      None
//
// Notes:
//  1. The instance variables, _filterPoles are modified.
// ############################# Private Method ###############################
void AnalogFilter::lowPassToHighPass(double wc)
{
  int       n;
  Complex   ctemp;

  n             = _numberZeros = _numberPoles = _filterOrder;
/***************************
 * pole -> wc/pole         *
 * zeros: s**n             *
 ***************************/
  ctemp                 = Complex(wc, 0.);
  while(n--)
  {
    if(norm(_filterPoles[n])>0.)
      _filterPoles[n]   = ctemp/_filterPoles[n];
    _filterZeros[n]     = Complex(0., 0.);
  }

  return;
}

// ############################# Private Method ###############################
// lowPassToBandPass -- This routine moves the poles of a filter having a cut off
//                      frequency of  1 rad/s to wc rad/s, and transforms them to
//                      band pass poles.
//
// Input:               wo:     New center frequency in rad/s
//                      bw:     Two-sided bandwidth in rad/s
//          
// Output:                      None
//
// Notes:
//  1. The instance variables, _filterPoles and _filterZeros are modified.
// ############################# Private Method ###############################
void AnalogFilter::lowPassToBandPass(double wo, double bw)
{
  int     i, j, k, n;
  double  wo_lp, q_lp, delta, delta_squared, q_bp, wo_bp[2];
  double  real_part, imag_part, temp1, temp2;
  Complex *new_poles;
 
/***************************
 * Find # of conjugate low *
 * pass poles:             *
 ***************************/
  _numberZeros  = n = _filterOrder;
  k             = _filterOrder/2;
  _numberPoles  = n+n;
  real_part     = imag_part = temp1 = temp2 = 0.;
  
/***************************
 * Note: Each low pass pole*
 * produces two bp poles   *
 ***************************/
  new_poles     = new Complex[n];

/***************************
 * Find bandpass poles     *
 * as in p. 580, Budak     *
 ***************************/
  i = 0;
  while(k--)
  {
    wo_lp           = abs(_filterPoles[k]);
    real_part       = real(_filterPoles[i]);
    if(real_part != 0.)
      q_lp          = -wo_lp/2./real_part;
    else
    {
      wo_lp         = 1.;                                   // Error condition
      q_lp          = 1.;
    }
    delta           = bw*wo_lp/wo;
    delta_squared   = delta*delta;
    q_bp            = sqrt((1+4./delta_squared)*(1+4./delta_squared) -4./delta_squared/q_lp/q_lp);
    q_bp            = q_lp*sqrt(1+4./delta_squared + q_bp)/SQRT2;
    wo_bp[0]        = 0.5*wo*(delta*q_bp/q_lp - sqrt(delta_squared*q_bp*q_bp/q_lp/q_lp -4.) );
    wo_bp[1]        = 0.5*wo*(delta*q_bp/q_lp + sqrt(delta_squared*q_bp*q_bp/q_lp/q_lp -4.) );
    for(j=0; j<2; j++)
    {
      real_part     = -wo_bp[j]/2./q_bp;
      imag_part     = sqrt(wo_bp[j]*wo_bp[j] - real_part*real_part);
      new_poles[i]  = Complex(real_part, imag_part);
      i++;
    }
  }
/**********************
 * Real axis pole:    *
 **********************/
  if(n%2)
  {
    real_part           = real(bw*_filterPoles[n/2]);
    temp1               = 4.*wo*wo;
    temp2               = real_part*real_part;
    if(temp1>=temp2)
    {
      imag_part         = sqrt(temp1 - temp2)/2.;
      new_poles[i]      = Complex(real_part/2., imag_part);
    }
    else
    {
      _realAxisPoles    = YES;
      new_poles[i]      = Complex((real_part+sqrt(temp2-temp1))/2., 0.);
    }
  }
/*********************
 * Copy over poles   *
 * zeros: s**(n/2)   *
 *********************/
  for(i=0; i<n; i++)
  {
    _filterPoles[i]     = new_poles[i];
    _filterZeros[i]     = Complex(0., 0.);
  }
  for( ; i<_numberPoles; i++)                                // conjugate poles
    _filterPoles[i]     = conj(_filterPoles[i-n]);
  if(_realAxisPoles)
    _filterPoles[_numberPoles-1]        = Complex( (real_part - sqrt(temp2-temp1))/2., 0.);

  delete [] new_poles;
  return;
}

// ############################# Private Method ###############################
// lowPassToBandStop -- This routine moves the poles of a filter having a cut off
//                      frequency of  1 rad/s to wc rad/s, and transforms them to
//                      band stop poles.
//
// Input:               wo:     New center frequency in rad/s
//                      bw:     Two-sided bandwidth in rad/s
//          
// Output:                      None
//
// Notes:
//  1. The instance variables, _filterPoles and _filterZeros are modified.
// ############################# Private Method ###############################
void AnalogFilter::lowPassToBandStop(double wo, double bw)
{
  int     i, j, k, n;
  double  wo_lp, q_lp, delta, delta_squared, q_bs, wo_bs[2];
  double  real_part, imag_part, temp1, temp2;
  Complex *new_poles;
 
/***************************
 * Find # of conjugate low *
 * pass poles:             *
 ***************************/
  n             = _filterOrder;
  k             = _filterOrder/2;
  _numberZeros  = _numberPoles = n+n;
  real_part     = imag_part = temp1 = temp2 = 0.;
  
/***************************
 * Note: Each low pass pole*
 * produces two bs poles   *
 ***************************/
  new_poles = new Complex[n];

/***************************
 * Find bandstop poles     *
 * as in p. 630, Budak     *
 ***************************/
  i = 0;
  while(k--)
  {
    wo_lp           = abs(_filterPoles[k]);
    real_part       = real(_filterPoles[i]);
    if(real_part != 0.)
      q_lp          = -wo_lp/2./real_part;
    else
    {
      wo_lp         = 1.;                                   // Error condition
      q_lp          = 1.;
    }
    delta           = bw/wo_lp/wo;
    delta_squared   = delta*delta;
    q_bs            = sqrt((1+4./delta_squared)*(1+4./delta_squared) -4./delta_squared/q_lp/q_lp);
    q_bs            = q_lp*sqrt(1+4./delta_squared + q_bs)/SQRT2;
    wo_bs[0]        = 0.5*wo*(delta*q_bs/q_lp - sqrt(delta_squared*q_bs*q_bs/q_lp/q_lp -4.) );
    wo_bs[1]        = 0.5*wo*(delta*q_bs/q_lp + sqrt(delta_squared*q_bs*q_bs/q_lp/q_lp -4.) );
    for(j=0; j<2; j++)
    {
      real_part     = -wo_bs[j]/2./q_bs;
      imag_part     = sqrt(wo_bs[j]*wo_bs[j] - real_part*real_part);
      new_poles[i]  = Complex(real_part, imag_part);
      i++;
    }
  }
/**********************
 * Real axis pole:    *
 **********************/
  if(n%2)
  {
    if(real(_filterPoles[n/2]) != 0.)
    {
      real_part = real(bw/_filterPoles[n/2]);
      temp1     = 4.*wo*wo;
      temp2     = real_part*real_part;
      if(temp1>=temp2)
      {
        imag_part       = sqrt(temp1 - temp2)/2.;
        new_poles[i]    = Complex(real_part/2., imag_part);
      }
      else
      {
        _realAxisPoles  = YES;
        new_poles[i]    = Complex((real_part+sqrt(temp2-temp1))/2., 0.);
      }
    }
  }
/**********************
 * Copy over poles,   *
 * zeros:             *
 * (s-jwo)(s+jwo)     *
 **********************/
  for(i=0; i<n; i++)
  {
    _filterPoles[i]     = new_poles[i];
    _filterZeros[i]     = Complex(0., -wo);
  }
  for( ; i<_numberPoles; i++)        /* conjugate poles */
  {
    _filterPoles[i]     = conj(_filterPoles[i-n]);
    _filterZeros[i]     = Complex(0., wo);
  }
  if(_realAxisPoles)
    _filterPoles[_numberPoles-1] = Complex( (real_part - sqrt(temp2-temp1))/2., 0.);
  
  delete [] new_poles;
  return;
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
void AnalogFilter::initInstanceVariables()
{
  double tol1;
  
  _realAxisPoles        = NO;
  setPassbandRipple(DEFAULT_RIPPLE);
//
//  Also, find machine eps
//
  _machineEps   = 1.;
  tol1          = 1.1;
  while(tol1 > 1.)
  {
    _machineEps /= 2.;
    tol1        = 1. + _machineEps;
  }
  
  return;
}
// ############################# Private Method ###############################
// zeroIn -- Find a zero in [ax, bx] of a function supplied by the sender
// Input:           tol:        tolerance on zero finder
//                  ax:         lower limit of search region
//                  bx:         upper limit of search region
//          
// Output:                      Value of x such that f(x) = 0.
// ############################# Private Method ###############################
double AnalogFilter::zeroIn(double ax, double bx, double tol)
{
  int    begin_step, converged, bisection;
  double a, b, c, d, e, fa, fb, fd, tol1;
  double xm, p, q, r, s;

/********************
 * Initialization   *
 ********************/
  c = d = e = fd = 0.;
  a  = ax;
  b  = bx;
  fa = functionToFindRoot(a);
  fb = functionToFindRoot(b);
  
/*******************
 * Begin step:     *
 *******************/
  converged  = 0;
  begin_step = 1;
  while(!converged)
  {
    if(begin_step)
    {
      c  = a;
      fd = fa;
      d  = b - a;
      e  = d;
    }
    if( ABS(fd) < ABS(fb) )
    {
      a  = b;
      b  = c;
      c  = a;
      fa = fb;
      fb = fd;
      fd = fa;
    }
  
/*******************
 * Convergence test*
 *******************/
  tol1 = 2.*_machineEps*ABS(b) + 0.5*tol;
  xm   = 0.5*(c-b);
  if( (ABS(xm)<=tol1) || (fb == 0.) )
  {
    converged = 1;
    break;
  }

/********************
 * Bisection        *
 * necessary?       *
 ********************/
  bisection = 1;
  if( (ABS(e)>=tol1) && (ABS(fa)>ABS(fb)) )
  {
  
/********************
 * quadratic interp *
 * possible?        *
 ********************/
     if( a == c)
     {
/********************
 * linear interp.   *
 ********************/
        s = fb/fa;
        p = 2.*xm*s;
        q = 1. - s;
     }
     else
     {
/*******************
 * Inverse quadrat.*
 * interp.         *
 *******************/
       q = fa/fd;
       r = fb/fd;
       s = fb/fa;
       p = s*(2.*xm*q*(q-r) - (b-a)*(r-1.));
       q = (q-1.)*(r-1.)*(s-1.);
     }
/*******************
 * Adjust signs:   *
 *******************/
     q = (p > 0.) ? -q : q;
     p = ABS(p);
/*******************
 * Is interpolation*
 * acceptable?     *
 *******************/
     if(  ( 2.*p < (2.*xm*q - ABS(tol1*q)) )  &&
          ( p < ABS(0.5*e*q) )                 )
     {
       e = d;
       d = p/q;
       bisection = 0;
     }
     else
       bisection = 1;
   }
/*******************
 * Bisection       *
 *******************/
   if(bisection)
   {
     d = xm;
     e = d;
   }
/*******************
 * Complete step   *
 *******************/
     a  = b;
     fa = fb;
     if( ABS(d) > tol1)
       b += d;
     if( ABS(d) <= tol1)
       b += SIGN(tol1, xm);
     fb = functionToFindRoot(b);
     if(fb*fd/ABS(fd) > 0.)
       begin_step = 1;
     else
       begin_step = 0;
   }
   return b;
}

// ############################# Private Method ###############################
// functionToFindRoot --  Used in besselPrototype to find 3 dB point
//
// Input:       x:          Frequency to evaluate filter response
//          
// Output:                  Filter response at the given frequency.
//
// Notes:
// ############################# Private Method ###############################
double AnalogFilter::functionToFindRoot(double x)
{
  Complex ctemp, response;
  double  output;
  
  ctemp     = Complex(0., x);
  response  = poleZeroResponseAt(ctemp, _filterPoles, _filterZeros);
  output    = norm(response);
  output    /= _passBandGain*_passBandGain;
  output    -= 0.5;
 
  return output;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the AnalogFilter class.
//
// Input:           aCoeffs:    The coefficients of the denominator polynomial
//                  bCoeffs:    The coefficients of the numerator polynomial
//                  order:      The order of the polynomials
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
AnalogFilter::AnalogFilter(double *aCoeffs, double *bCoeffs, int order)
             :RealFilter(aCoeffs, bCoeffs, order)

{
  initInstanceVariables();
  return;
}
  
// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the AnalogFilter class.
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
AnalogFilter::AnalogFilter(int type, double centerFreq, double cutoffFreq, int order)
             :RealFilter(type, centerFreq, cutoffFreq, order)
{
  initInstanceVariables();
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the AnalogFilter class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// 1. The base class destructor will also be called.
// ############################# Class Destructor ###############################
AnalogFilter::~AnalogFilter()
{
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
void AnalogFilter::setPassBandGain(LOGICAL usePolesZeros)
{
  Complex   arg, num, den, response;
  double    wo, omega;
//
// Calculate the pass band gain:
//
  switch(_filterPassType)
  {
    case LOW_PASS:
    case BAND_STOP:
    default:
      arg               = Complex(0., 0.);                      // omega = 0
      break;
    case HIGH_PASS:
      omega             = TWOPI*2.*_fc;
      arg               = Complex(0., omega);
      break;
    case BAND_PASS:
      wo                = TWOPI*_fo;
      arg               = Complex(0., wo);                      // omega = wo;
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
// filterPoles -- Finds and returns the complex poles of the filter
// Input:                   None
//          
// Output:                  The complex poles
// ############################# Public Method ###############################
const Complex *AnalogFilter::filterPoles()
{
  findPolesZeros();
  return _filterPoles;
}

// ############################# Public Method ###############################
// filterZeros -- Finds and returns the complex zeros of the filter
// Input:                   None
//          
// Output:                  The complex poles
// ############################# Public Method ###############################
const Complex *AnalogFilter::filterZeros()
{
  findPolesZeros();
  return _filterZeros;
}

// ############################# Public Method ###############################
// filterResponseAtFrequencies -- Calculates the filter's response at input frequencies
//
// Input:   frequencies:    array of input frequencies given in rad/s
//          numberPoints    size of the above array
//          
// Output:                  _filterResponse
// ############################# Public Method ###############################
float  *AnalogFilter::filterResponseAtFrequencies(float *frequencies, int numberPts, LOGICAL calc)
{
  int     n, m;
  Complex complex_response, carg;
//
//  Error check:
//
 if( (numberPts<1) || (frequencies==NULL) )
   return NULL;

  if(calc)
  {
    findPolesZeros();
    if(_filterStructureType == TRANSFER_FUNCTION)
      transferFromPoles(_filterPoles, _filterZeros, _numberPoles, _numberZeros, _realAxisPoles);
  }
//
// _deltaOmega is used in outputResponseFor() for calculating group delay:
// _subAngle is used for PI crossings.
//
  if(numberPts>1)
    _deltaOmega  = (frequencies[numberPts-1] - frequencies[0])/numberPts;
  _deltaOmega    = MAX(1e-30, _deltaOmega);
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
  for(n=-2; n<numberPts; n++)
  {
    m       = MAX(0, n);
    carg    = Complex(0.,frequencies[m]);
    if(_filterStructureType == POLE_ZERO)
      complex_response  = poleZeroResponseAt(carg, _filterPoles, _filterZeros);
    else 
      complex_response  = transferResponseAt(carg);
    
    _filterResponse[m]   = outputResponseFor(complex_response, frequencies[m]);
  }

  return _filterResponse;
}

