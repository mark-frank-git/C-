/************************************************************************
 *                                                                      *
 * This class implements a filter which gives the classical Doppler     *
 * spectrum with peaks at +/- the Doppler frequency.                    *
 *                                                                      *
 * File:RicianDopplerFilter.cc                                          *
 *                                                                      *
 * The filter is stored in the following form:                          *
 *                                                                      *
 *           b[n] + b[n-1]z + ... + b[0]z^n                             *
 *    H(z) = ------------------------------------                       *
 *           1    + a[n-1]z + ... + a[0]z^n                             *
 *                                                                      *
 * See: D.I. Laurenson, D.G.M. Cruickshank, and G.J.R. Povey, "A        *
 * computationally efficient multipath channel simulator for the COST   *
 * 207 models," IEE Colloquium on Computer Modeling of Communication    *
 * Systems, 1994, pp. 8/1-8/6.                                          *
 *                                                                      *
 ************************************************************************/

#include "RicianDopplerFilter.h"                                // Object prototypes

#if defined(WIN32)
#include <GNU/Complex.h>
#include <C_Libraries/constants.h>
#else
#include "Complex.h"
#include "constants.h"
#endif

#include <math.h>
#include <stdio.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))
#define SGN(a)    ( ((a)>=0.) ? 1 : -1 )

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the RicianDopplerFilter class.
//
// Input:       doppler:        Doppler shift frequency in Hertz
//              sampling:       Sampling frequency in Hertz
//
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
RicianDopplerFilter::RicianDopplerFilter(float doppler, float sampling)
                       :AbstractDopplerFilter()
{
  setDopplerFrequency(doppler);
  setSamplingFrequency(sampling);

  initialize();

  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the RicianDopplerFilter class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
RicianDopplerFilter::~RicianDopplerFilter()
{
  return;
}

// ############################# Public Method ###############################
// filterComplexData -- Filter an input data sample, and return the sampled value.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   input:          An input data sample
//          
// Output:                  The filtered data sample
//
// NOTES:
// 1. For now, we just return a complex (1., 0.), and perform no filtering
// ############################# Public Method ###############################
Complex RicianDopplerFilter::filterComplexData (Complex &input)
{
  Complex       yout;
//
// Filter the data:
//
  yout      = Complex(1., 0.);
  return yout;
}
