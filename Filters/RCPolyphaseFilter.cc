/************************************************************************
 *                                                                      *
 * This subclass of ComplexDigitalFilter adds functionality for         *
 * calculating the transfer function of RC polyphase filters.           *
 * The actual filtering functions are defined in the superclass.        *
 *                                                                      *
 * File:RCPolyphaseFilter.h                                             *
 *                                                                      *
 * The filter is stored in the following form:                          *
 *                                                                      *
 *           b[0] + b[1]z^(-1) + ... + b[n]z^-(n)                       *
 *  H(z)   = ------------------------------------                       *
 *           1    + a[1]z^(-1) + ... + a[n]z^-(n)                       *
 * Where the a and b polynomials are stored in the polynomial objects.  *
 * NOTE: The order of the filter is the # taps - 1.                     *
 *                                                                      *
 * SEE: "RC Sequence Asymmetric Polyphase Networks fro RF Integrated    *
 *       Transceivers," by Galal, Ragaie, Tawfik.                       *
 *                                                                      *
 * Revision history:                                                    *
 *  1.  03/19/02  - Started.                                            *
 ************************************************************************/

#include        "ComplexDigitalFilter.h"
#include        "RCPolyphaseFilter.h"                                   // Object prototypes

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


// ############################# Public Method ###############################
// findSTransferFunction() -- Finds the S domain transfer function based on
//                            number of stages.
// Input:                       None
//          
// Output:                      None
// Notes:
//  1. The instance variables _aPolyObject, and _bPolyObject are modified.
//  2. The s domain poly is of the form: p(s) = b[n] + b[n-1]s + ... + b[0]s**(n)
//  3. The notch frequency is set at a negative frequency instead of a positive
//     frequency as in the reference.
// ############################# Public Method ###############################
void RCPolyphaseFilter::findSTransferFunction()
{
  double        tau1, tau2, tau3, tau4, tau5, tau6;
  Complex       a_coeffs[4], b_coeffs[4];
//
// First calculate the time constants, and warp the frequencies according to the sampling
// frequency
//
  tau1                  = _r1*_c1;
  tau2                  = _r2*_c2;
  tau3                  = _r3*_c3;
  tau4                  = _r1*_c2;
  tau5                  = _r2*_c3;
  tau6                  = _r1*_c3;
  if(tau1 != 0.)
    tau1                = 1./warpFrequency(1./tau1, _fs);
  if(tau2 != 0.)
    tau2                = 1./warpFrequency(1./tau2, _fs);
  if(tau3 != 0.)
    tau3                = 1./warpFrequency(1./tau3, _fs);
  if(tau4 != 0.)
    tau4                = 1./warpFrequency(1./tau4, _fs);
  if(tau5 != 0.)
    tau5                = 1./warpFrequency(1./tau5, _fs);
  if(tau6 != 0.)
    tau6                = 1./warpFrequency(1./tau6, _fs);
//
// Now, find the a and b coeffs for the s domain filter:
//
  switch(_numberStages)
  {
    default:
    case 1:                     // Single stage
      setFilterOrder(1);
      b_coeffs[1]               = Complex(1., 0.);
      b_coeffs[0]               = Complex(0., -tau1);                   // negative notch freq
      a_coeffs[1]               = Complex(1., 0.);
      a_coeffs[0]               = Complex(tau1, 0.);
      break;  
    case 2:                     // 2 stage
      setFilterOrder(2);
      b_coeffs[2]               = Complex(1.,0.);
      b_coeffs[1]               = Complex(0., -(tau1+tau2));            // negative notch freq
      b_coeffs[0]               = Complex((-tau1*tau2), 0.);
      a_coeffs[2]               = Complex(1.,0.);
      a_coeffs[1]               = Complex((tau1+tau2+2.*tau4),0.);
      a_coeffs[0]               = Complex((tau1*tau2), 0.);
      break;  
    case 3:                     // 3 stage
      setFilterOrder(3);
      b_coeffs[3]               = Complex(1.,0.);
      b_coeffs[2]               = Complex(0., -(tau1+tau2+tau3));
      b_coeffs[1]               = Complex((-tau1*tau2-tau1*tau3-tau2*tau3), 0.);
      b_coeffs[0]               = Complex(0., tau1*tau2*tau3);
//
      a_coeffs[3]               = Complex(1., 0.);
      a_coeffs[2]               = Complex(tau1+tau2+tau3+2*tau4+2*tau5+2*tau6,0.);
      a_coeffs[1]               = Complex(tau2*tau3+tau1*tau2+tau1*tau3+2*tau4*tau3+2*tau6*tau2+2*tau1*tau5,0.);
      a_coeffs[0]               = Complex((tau1*tau2*tau3), 0.);
      break;  
  }
  setFilterACoeffs(a_coeffs);
  setFilterBCoeffs(b_coeffs);
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the RCPolyphaseFilter class.
//
// Input:           stages:             number of RC polyphase stages
//                  samplingFreq:       Digital sampling frequency in Hertz
//
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
RCPolyphaseFilter::RCPolyphaseFilter(int numberStages)
              :ComplexDigitalFilter()
{
  setNumberStages(numberStages);
  setSamplingFrequency(DEFAULT_SAMPLING_FREQ);
  
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the RCPolyphaseFilter class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
RCPolyphaseFilter::~RCPolyphaseFilter()
{
  return;
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
void RCPolyphaseFilter::findTransferFunction()
{
  findSTransferFunction();
  convertCoefficientsToZ();

  return;
}
