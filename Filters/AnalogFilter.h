#ifndef _ANALOGFILTER_H
#define _ANALOGFILTER_H 1
/************************************************************************
 *                                                                      *
 * This subclass of RealFilter implements an analog filter object.      *
 *                                                                      *
 * File:AnalogFilter.h                                                  *
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
 *  2. 07/12/95  - setFilterOrder -> setAnalogFilterOrder.              *
 *  3. 02/20/97  - Reversed order of denominator, numerator polys.      *
 *  4. 02/25/97  - Add getACoeffs(), getBCoeffs().                      *
 *  5. 12/05/97  - make setFilterACoeff() virtual.                      *
 *  6. 01/05/98  - make subclass of AbstractFilter.                     *
 *  7. 10/08/99  - make subclass of RealFilter.                         *
 *  8. 11/13/00  - Add _passbandRipple (used to be equal to 3 dB)       *
 ************************************************************************/
#include "RealFilter.h"
#define MAX_ORDER               15                      // max prototype order
#define DEFAULT_RIPPLE          3.                      // Ripple in dB

class   Complex;                                        // Class prototypes

class AnalogFilter: public RealFilter
{
protected:
  LOGICAL       _realAxisPoles;                         // real poles from bp, bs?
  double        _machineEps;                            // Machine epsilon
  double        _passbandRipple;                        // Passband ripple in dB

//
// Private methods:
//
  double        rippleToEps(double ripple);             // Converts ripple in dB to epsilon
  void          findPolesZeros();                       // Finds poles and zeros of filter
  void          findAnalogPolesZerosFor(double wc, double wo, double bw);
                                                        // Finds poles and zeros for the given inputs
  void          butterPrototype();                      // Finds a Butterworth prototype filter
  void          chebyPrototype();                       // Finds a Chebyshev prototype filter
  void          besselPrototype();                      // Finds a Bessel prototype filter
  void          lowPassToLowPass(double wc);            // Low pass to low pass transformation
  void          lowPassToHighPass(double wc);           // Low pass to high pass transformation
  void          lowPassToBandPass(double wo, double bw); // Low pass to band pass transformation
  void          lowPassToBandStop(double wo, double bw); // Low pass to band stop transformation
  void          initInstanceVariables();                // Initialize instance variables
  double        zeroIn(double ax, double bx, double tol=1.e-5);
  double        functionToFindRoot(double x);           // Used for finding root of polynomial
//
// Public methods:
//
public:
/********************************
 * Constructors, destructor     *
 *******************************/
  AnalogFilter(double *aCoeffs, double *bCoeffs, int order);
  AnalogFilter(int type, double centerFreq, double cutoffFreq, int order);
  ~AnalogFilter();

/**************************
 * Setting parameters:    *
 **************************/
  void          setPassBandGain(LOGICAL usePolesZeros);
  void          setPassbandRipple(double ripple)                {_passbandRipple = ripple; return;}

/**********************
 * Get parameters:    *
 **********************/
  LOGICAL       realAxisPoles()         {return _realAxisPoles;}
  const         Complex *filterPoles();
  const         Complex *filterZeros();

/**********************
 * Getting Outputs:   *
 **********************/
  float *filterResponseAtFrequencies(float *frequencies, int numberPoints, LOGICAL calc=YES);

};

#endif
